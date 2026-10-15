#!/usr/bin/env python3
"""
Add gnomAD allele frequencies to a MAF, fetched from Genome Nexus (the service cBioPortal
uses for its gnomAD column).

Usage: annotate_gnomad.py <maf> --out <out_file> [--url https://www.genomenexus.org]
                          [--batch-size 200] [--retries 3]

Adds (or fills empty cells of) two columns:
  gnomADe_AF  gnomAD v2.1.1 exomes, all populations
  gnomADg_AF  gnomAD v2.1.1 genomes, all populations
A variant gnomAD has never seen gets empty cells.

Only the variant itself is sent (chromosome, start, end, reference and alternate allele, as
written in the MAF), never sample IDs or other columns. Each distinct variant is sent once,
in batches. If Genome Nexus cannot be reached after the retries, the script exits with an
error instead of writing an unannotated file, so the gnomAD filter is never skipped silently.
"""
import argparse
import json
import sys
import time
import urllib.error
import urllib.request

NEW_COLUMNS = ('gnomADe_AF', 'gnomADg_AF')
NEEDED = ('Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2')


def genome_nexus_chrom(chrom):
    chrom = chrom[3:] if chrom.lower().startswith('chr') else chrom
    return 'MT' if chrom == 'M' else chrom


def post(url, variants, retries):
    """POST one batch to /annotation/genomic; return the parsed JSON list."""
    body = json.dumps([
        {'chromosome': c, 'start': int(s), 'end': int(e), 'referenceAllele': r, 'variantAllele': a}
        for c, s, e, r, a in variants
    ]).encode()
    request = urllib.request.Request(
        f"{url.rstrip('/')}/annotation/genomic?fields=my_variant_info",
        data=body,
        headers={'Content-Type': 'application/json', 'Accept': 'application/json'},
    )
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                return json.load(response)
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as err:
            if attempt == retries:
                sys.exit(f"ERROR: Genome Nexus request to {url} failed after {retries} attempts: {err}")
            print(f"WARNING: Genome Nexus request failed ({err}); retrying", file=sys.stderr)
            time.sleep(5 * attempt)


def allele_frequency(annotation, source):
    """annotation.my_variant_info.annotation.<source>.alleleFrequency.af, or ''."""
    info = (annotation.get('my_variant_info') or {}).get('annotation') or {}
    af = ((info.get(source) or {}).get('alleleFrequency') or {}).get('af')
    return '' if af is None else str(af)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('maf')
    parser.add_argument('--out', required=True)
    parser.add_argument('--url', default='https://www.genomenexus.org')
    parser.add_argument('--batch-size', type=int, default=200)
    parser.add_argument('--retries', type=int, default=3)
    args = parser.parse_args()

    with open(args.maf) as fh:
        lines = [line.rstrip('\r\n') for line in fh if line.strip()]
    if not lines:
        sys.exit(f"ERROR: {args.maf} is empty")

    header = lines[0].split('\t')
    missing = [name for name in NEEDED if name not in header]
    if missing:
        sys.exit(f"ERROR: MAF lacks {', '.join(missing)}")
    for name in NEW_COLUMNS:
        if name not in header:
            header.append(name)
    col = {name: i for i, name in enumerate(header)}

    rows = [line.split('\t') for line in lines[1:]]
    keys = []
    for row in rows:
        row.extend([''] * (len(header) - len(row)))
        start, end = row[col['Start_Position']], row[col['End_Position']]
        if not (start.isdigit() and end.isdigit()):
            keys.append(None)
            continue
        keys.append((genome_nexus_chrom(row[col['Chromosome']]), start, end,
                     row[col['Reference_Allele']], row[col['Tumor_Seq_Allele2']]))

    # Genome Nexus echoes each query back as "chrom,start,end,ref,alt"
    variants = sorted({k for k in keys if k})
    frequencies = {}
    for i in range(0, len(variants), args.batch_size):
        for annotation in post(args.url, variants[i:i + args.batch_size], args.retries):
            query = annotation.get('originalVariantQuery')
            if query:
                frequencies[query] = (allele_frequency(annotation, 'gnomadExome'),
                                      allele_frequency(annotation, 'gnomadGenome'))

    found = 0
    for row, key in zip(rows, keys):
        if not key:
            continue
        exome, genome = frequencies.get(','.join(key), ('', ''))
        if exome or genome:
            found += 1
        for name, value in zip(NEW_COLUMNS, (exome, genome)):
            if row[col[name]] in ('', '.'):
                row[col[name]] = value

    with open(args.out, 'w') as out:
        out.write('\t'.join(header) + '\n')
        for row in rows:
            out.write('\t'.join(row) + '\n')

    print(f"annotate_gnomad: {len(variants)} distinct variants sent to {args.url}; "
          f"{found}/{len(rows)} MAF rows have a gnomAD frequency", file=sys.stderr)


if __name__ == '__main__':
    main()
