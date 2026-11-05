#!/usr/bin/env python3
"""
Add gnomAD allele frequencies and the ClinVar classification to a MAF, fetched from Genome
Nexus (the service cBioPortal uses for its gnomAD and ClinVar columns).

Usage: annotate_gnomad.py <maf> --out <out_file> [--url https://www.genomenexus.org]
                          [--batch-size 200] [--retries 3]
                          [--empty-retries 2] [--retry-wait 10]

Writes (overwriting any existing values) these columns:
  gnomADe_AF          gnomAD v2.1.1 exomes, all populations
  gnomADg_AF          gnomAD v2.1.1 genomes, all populations
  gnomAD.AF           exomes + genomes combined, (AC_e + AC_g) / (AN_e + AN_g): the "Total"
                      cBioPortal shows; the gnomAD namespace makes it a visible MAF column
  gnomAD.AF_max       highest of the total and the 8 ancestry AFs, exomes and genomes apart
  gnomAD.AF_max_pop   where that maximum comes from, e.g. exome_afr or genome_total
  ClinVar.CLNSIG      ClinVar clinicalSignificance, e.g. Pathogenic/Likely_pathogenic
  ClinVar.CLNSIGCONF  the breakdown of a conflicting call, e.g. Pathogenic(3)|Benign(1)
A frequency missing from the response is computed from AC / AN, as cBioPortal does.
A variant gnomAD (or ClinVar) has never seen gets empty cells.

Genome Nexus sometimes answers without its gnomAD block, with no error. A variant without
gnomAD or ClinVar data is therefore asked again in up to --empty-retries extra passes; only
one that is still empty after all of them is written as unknown.

Only the variant itself is sent (chromosome, start, end, reference and alternate allele, as
written in the MAF), never sample IDs or other columns. Each distinct variant is sent once,
in batches. If Genome Nexus cannot be reached after the retries, the script exits with an
error instead of writing an unannotated file, so the filters are never skipped silently.
"""
import argparse
import json
import sys
import time
import urllib.error
import urllib.request

NEW_COLUMNS = ('gnomADe_AF', 'gnomADg_AF', 'gnomAD.AF', 'gnomAD.AF_max', 'gnomAD.AF_max_pop',
               'ClinVar.CLNSIG', 'ClinVar.CLNSIGCONF')
POPULATIONS = ('afr', 'amr', 'asj', 'eas', 'fin', 'nfe', 'oth', 'sas')
NEEDED = ('Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2')
RETRY_BATCH_SIZE = 50
NO_ANSWER = ({}, {})


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
        f"{url.rstrip('/')}/annotation/genomic?fields=my_variant_info,clinvar",
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
    """Annotate variants in batches; return {"chrom,start,end,ref,alt": (my_variant_info, clinvar)}."""
    found = {}
    for i in range(0, len(variants), batch_size):
        for annotation in post(url, variants[i:i + batch_size], retries):
            # Genome Nexus echoes each query back as "chrom,start,end,ref,alt"
            key = annotation.get('originalVariantQuery')
            if key:
                found[key] = ((annotation.get('my_variant_info') or {}).get('annotation') or {},
                              (annotation.get('clinvar') or {}).get('annotation') or {})
    return found


def has_gnomad(info):
    return bool(info) and ('gnomadExome' in info or 'gnomadGenome' in info)


def complete(answer):
    """True when the answer has both a gnomAD and a ClinVar block; else it may be a dropped answer."""
    return bool(answer) and has_gnomad(answer[0]) and bool(answer[1])


def source_frequency(info, source, pop=''):
    """(af, ac, an) of info.<source> for one population ('' = all); af falls back to ac / an."""
    data = info.get(source) or {}
    suffix = f'_{pop}' if pop else ''
    af = (data.get('alleleFrequency') or {}).get(f'af{suffix}')
    ac = (data.get('alleleCount') or {}).get(f'ac{suffix}')
    an = (data.get('alleleNumber') or {}).get(f'an{suffix}')
    if af is None and ac is not None and an:
        af = ac / an
    return af, ac, an


def max_frequency(info):
    """(highest AF, 'exome_afr'-style label) over the total and every population of both sources."""
    best = (None, '')
    for source, label in (('gnomadExome', 'exome'), ('gnomadGenome', 'genome')):
        for pop in ('',) + POPULATIONS:
            af = source_frequency(info, source, pop)[0]
            if af is not None and (best[0] is None or af > best[0]):
                best = (af, f"{label}_{pop or 'total'}")
    return best


def frequencies(info):
    """(exome AF, genome AF, combined AF, max AF, max AF label) as strings, '' when unknown."""
    exome, genome = source_frequency(info, 'gnomadExome'), source_frequency(info, 'gnomadGenome')
    counted = [(ac, an) for _, ac, an in (exome, genome) if ac is not None and an]
    if counted:
        combined = sum(ac for ac, _ in counted) / sum(an for _, an in counted)
    else:
        afs = [af for af, _, _ in (exome, genome) if af is not None]
        combined = max(afs) if afs else None
    top, top_pop = max_frequency(info)
    return tuple('' if af is None else format(af, '.6g') for af in (exome[0], genome[0], combined, top)) + (top_pop,)


def clinvar_values(clinvar):
    """(clinicalSignificance, conflictingClinicalSignificance), '' when absent."""
    return (clinvar.get('clinicalSignificance') or '', clinvar.get('conflictingClinicalSignificance') or '')


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

    # A missing gnomAD or ClinVar block can be a dropped answer: ask again before calling it unknown
    empty = [k for k in variants if not complete(info.get(k))]
    requeried, recovered = len(empty), 0
    for attempt in range(1, args.empty_retries + 1):
        if not empty:
            break
        print(f"WARNING: {len(empty)} variant(s) came back without gnomAD or ClinVar data; "
              f"asking again (pass {attempt}/{args.empty_retries})", file=sys.stderr)
        time.sleep(args.retry_wait * attempt)
        again = query(args.url, [as_tuple[k] for k in empty], RETRY_BATCH_SIZE, args.retries)
        for k in empty:
            old, new = info.get(k) or NO_ANSWER, again.get(k) or NO_ANSWER
            # Keep whichever half each answer has, so a later partial answer never loses data
            info[k] = (new[0] if has_gnomad(new[0]) else old[0], new[1] or old[1])
            if complete(info[k]):
                recovered += 1
        empty = [k for k in empty if not complete(info.get(k))]

    no_gnomad = [k for k in variants if not has_gnomad((info.get(k) or NO_ANSWER)[0])]
    no_clinvar = sum(1 for k in variants if not (info.get(k) or NO_ANSWER)[1])

    # A dbSNP variant absent from gnomAD is possible but worth a manual look
    with_rsid = [k for k in no_gnomad if isinstance((info.get(k) or NO_ANSWER)[0].get('dbsnp'), dict)
                 and info[k][0]['dbsnp'].get('rsid')]
    if with_rsid:
        print(f"WARNING: {len(with_rsid)} variant(s) with a dbSNP rsID still have no gnomAD data: "
              f"{', '.join(with_rsid[:10])}{' ...' if len(with_rsid) > 10 else ''}", file=sys.stderr)

    found = 0
    for row, key in zip(rows, keys):
        gnomad, clinvar = (info.get(key) or NO_ANSWER) if key else NO_ANSWER
        values = frequencies(gnomad) + clinvar_values(clinvar)
        if any(values[:4]):
            found += 1
        for name, value in zip(NEW_COLUMNS, values):
            row[col[name]] = value

    with open(args.out, 'w') as out:
        out.write('\t'.join(header) + '\n')
        for row in rows:
            out.write('\t'.join(row) + '\n')

    print(f"annotate_gnomad: {len(variants)} distinct variants sent to {args.url}; "
          f"{requeried} asked again, {recovered} recovered, {len(no_gnomad)} not in gnomAD, "
          f"{no_clinvar} not in ClinVar; {found}/{len(rows)} MAF rows have a gnomAD frequency",
          file=sys.stderr)


if __name__ == '__main__':
    main()
