#!/usr/bin/env python3
"""
Add gnomAD allele frequencies to a MAF, fetched from Genome Nexus (the service cBioPortal
uses for its gnomAD column).

Usage: annotate_gnomad.py <maf> --out <out_file> [--url https://www.genomenexus.org]
                          [--batch-size 200] [--retries 3]
                          [--empty-retries 2] [--retry-wait 10]

Writes (overwriting any existing values) three columns:
  gnomADe_AF  gnomAD v2.1.1 exomes, all populations
  gnomADg_AF  gnomAD v2.1.1 genomes, all populations
  gnomAD.AF   exomes + genomes combined, (AC_e + AC_g) / (AN_e + AN_g): the "Total"
              cBioPortal shows; the gnomAD namespace makes it a visible MAF column
A frequency missing from the response is computed from AC / AN, as cBioPortal does.
A variant gnomAD has never seen gets empty cells.

Genome Nexus sometimes answers without its gnomAD block, with no error. A variant without
gnomAD data is therefore asked again in up to --empty-retries extra passes; only one that
is still empty after all of them is written as unknown.

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

NEW_COLUMNS = ('gnomADe_AF', 'gnomADg_AF', 'gnomAD.AF')
NEEDED = ('Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2')
RETRY_BATCH_SIZE = 50


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


def query(url, variants, batch_size, retries):
    """Annotate variants in batches; return {"chrom,start,end,ref,alt": my_variant_info.annotation}."""
    found = {}
    for i in range(0, len(variants), batch_size):
        for annotation in post(url, variants[i:i + batch_size], retries):
            # Genome Nexus echoes each query back as "chrom,start,end,ref,alt"
            key = annotation.get('originalVariantQuery')
            if key:
                found[key] = (annotation.get('my_variant_info') or {}).get('annotation') or {}
    return found


def has_gnomad(info):
    return bool(info) and ('gnomadExome' in info or 'gnomadGenome' in info)


def source_frequency(info, source):
    """(af, ac, an) of info.<source>; af falls back to ac / an."""
    data = info.get(source) or {}
    af = (data.get('alleleFrequency') or {}).get('af')
    ac = (data.get('alleleCount') or {}).get('ac')
    an = (data.get('alleleNumber') or {}).get('an')
    if af is None and ac is not None and an:
        af = ac / an
    return af, ac, an


def frequencies(info):
    """(exome AF, genome AF, combined AF) as strings, '' when unknown."""
    exome, genome = source_frequency(info, 'gnomadExome'), source_frequency(info, 'gnomadGenome')
    counted = [(ac, an) for _, ac, an in (exome, genome) if ac is not None and an]
    if counted:
        combined = sum(ac for ac, _ in counted) / sum(an for _, an in counted)
    else:
        afs = [af for af, _, _ in (exome, genome) if af is not None]
        combined = max(afs) if afs else None
    return tuple('' if af is None else format(af, '.6g') for af in (exome[0], genome[0], combined))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('maf')
    parser.add_argument('--out', required=True)
    parser.add_argument('--url', default='https://www.genomenexus.org')
    parser.add_argument('--batch-size', type=int, default=200)
    parser.add_argument('--retries', type=int, default=3)
    parser.add_argument('--empty-retries', type=int, default=2)
    parser.add_argument('--retry-wait', type=float, default=10)
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
        keys.append(','.join((genome_nexus_chrom(row[col['Chromosome']]), start, end,
                              row[col['Reference_Allele']], row[col['Tumor_Seq_Allele2']])))

    variants = sorted({k for k in keys if k})
    as_tuple = {k: tuple(k.split(',')) for k in variants}
    info = query(args.url, [as_tuple[k] for k in variants], args.batch_size, args.retries)

    # A missing gnomAD block can be a dropped answer: ask again before calling it unknown
    empty = [k for k in variants if not has_gnomad(info.get(k))]
    requeried, recovered = len(empty), 0
    for attempt in range(1, args.empty_retries + 1):
        if not empty:
            break
        print(f"WARNING: {len(empty)} variant(s) came back without gnomAD data; "
              f"asking again (pass {attempt}/{args.empty_retries})", file=sys.stderr)
        time.sleep(args.retry_wait * attempt)
        again = query(args.url, [as_tuple[k] for k in empty], RETRY_BATCH_SIZE, args.retries)
        for k in empty:
            if has_gnomad(again.get(k)):
                info[k] = again[k]
                recovered += 1
        empty = [k for k in empty if not has_gnomad(info.get(k))]

    # A dbSNP variant absent from gnomAD is possible but worth a manual look
    with_rsid = [k for k in empty if isinstance((info.get(k) or {}).get('dbsnp'), dict)
                 and info[k]['dbsnp'].get('rsid')]
    if with_rsid:
        print(f"WARNING: {len(with_rsid)} variant(s) with a dbSNP rsID still have no gnomAD data: "
              f"{', '.join(with_rsid[:10])}{' ...' if len(with_rsid) > 10 else ''}", file=sys.stderr)

    found = 0
    for row, key in zip(rows, keys):
        values = frequencies(info.get(key) or {}) if key else ('', '', '')
        if any(values):
            found += 1
        for name, value in zip(NEW_COLUMNS, values):
            row[col[name]] = value

    with open(args.out, 'w') as out:
        out.write('\t'.join(header) + '\n')
        for row in rows:
            out.write('\t'.join(row) + '\n')

    print(f"annotate_gnomad: {len(variants)} distinct variants sent to {args.url}; "
          f"{requeried} asked again, {recovered} recovered, {len(empty)} not in gnomAD; "
          f"{found}/{len(rows)} MAF rows have a gnomAD frequency", file=sys.stderr)


if __name__ == '__main__':
    main()
