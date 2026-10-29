#!/usr/bin/env python3
"""
Filter a MAF by gnomAD population allele frequency, tumour VAF and/or ClinVar classification.
Usage: filter_mutations_af.py <maf> --out <out_file>
                              [--max-gnomad-af X] [--gnomad-columns a,b,c]
                              [--max-vaf Y]
                              [--clinvar-keep Pathogenic,Likely_pathogenic,Uncertain_significance]

gnomAD:  a row is removed when the highest value among the listed gnomAD
         columns present in the MAF is > X. The default gnomAD.AF_max is the
         highest ancestry AF, so one population above X is enough.
VAF:     a row is removed when t_alt_count / t_depth > Y.
ClinVar: a row is kept only when ClinVar.CLNSIG holds one of the listed terms
         (any part of a combined call such as Pathogenic/Likely_pathogenic, and
         Pathogenic,_low_penetrance-style variants). A conflicting call is kept when
         ClinVar.CLNSIGCONF lists one of them. A variant without ClinVar data is removed,
         and so is every row when the ClinVar.CLNSIG column is absent.

Missing or non-numeric values (blank, '.') never remove a row for gnomAD/VAF. If a
gnomAD/VAF column is absent from the header, that criterion is skipped with a warning.
"""
import argparse
import re
import sys

CLINVAR_COLUMN = 'ClinVar.CLNSIG'
CLINVAR_CONFLICT_COLUMN = 'ClinVar.CLNSIGCONF'


def to_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def clinvar_terms(value):
    """'Pathogenic(3)|Likely pathogenic(1)' -> {'pathogenic', 'likely_pathogenic'}."""
    terms = set()
    for part in re.split(r'[/|;]', value or ''):
        part = re.sub(r'\(\d+\)', '', part).strip().lower().replace(' ', '_')
        if part:
            terms.add(part)
    return terms


def clinvar_kept(clnsig, conflict, keep):
    """True when the call, or for a conflicting call its breakdown, has a kept term."""
    def matches(terms):
        # Prefix match keeps qualified calls such as pathogenic,_low_penetrance
        return any(t == k or t.startswith(k + ',') for t in terms for k in keep)

    terms = clinvar_terms(clnsig)
    if matches(terms):
        return True
    return any(t.startswith('conflicting') for t in terms) and matches(clinvar_terms(conflict))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('maf')
    parser.add_argument('--out', required=True)
    parser.add_argument('--max-gnomad-af', type=float)
    parser.add_argument('--gnomad-columns', default='gnomAD.AF_max')
    parser.add_argument('--max-vaf', type=float)
    parser.add_argument('--clinvar-keep')
    args = parser.parse_args()

    with open(args.maf) as fh:
        lines = fh.read().splitlines()

    if not lines:
        open(args.out, 'w').close()
        return

    header = lines[0].split('\t')
    col = {name: i for i, name in enumerate(header)}

    gnomad_idx = []
    if args.max_gnomad_af is not None:
        wanted = [c.strip() for c in args.gnomad_columns.split(',') if c.strip()]
        gnomad_idx = [col[c] for c in wanted if c in col]
        if not gnomad_idx:
            print(f"WARNING: none of the gnomAD columns {wanted} found in MAF; "
                  "gnomAD filter skipped", file=sys.stderr)

    vaf_idx = None
    if args.max_vaf is not None:
        if 't_depth' in col and 't_alt_count' in col:
            vaf_idx = (col['t_alt_count'], col['t_depth'])
        else:
            print("WARNING: t_depth/t_alt_count not found in MAF; VAF filter skipped",
                  file=sys.stderr)

    clinvar_keep = None
    if args.clinvar_keep:
        clinvar_keep = clinvar_terms(args.clinvar_keep.replace(',', '|'))
        if CLINVAR_COLUMN not in col:
            # No ClinVar data for any row: every row is dropped, like a variant not in ClinVar
            print(f"WARNING: {CLINVAR_COLUMN} not found in MAF; every row dropped by the ClinVar filter",
                  file=sys.stderr)

    removed_gnomad = 0
    removed_vaf = 0
    removed_clinvar = 0
    kept = [lines[0]]

    for line in lines[1:]:
        if not line.strip():
            continue
        fields = line.split('\t')

        if gnomad_idx:
            afs = [to_float(fields[i]) for i in gnomad_idx if i < len(fields)]
            afs = [af for af in afs if af is not None]
            if afs and max(afs) > args.max_gnomad_af:
                removed_gnomad += 1
                continue

        if vaf_idx:
            alt_i, depth_i = vaf_idx
            alt = to_float(fields[alt_i]) if alt_i < len(fields) else None
            depth = to_float(fields[depth_i]) if depth_i < len(fields) else None
            if alt is not None and depth and alt / depth > args.max_vaf:
                removed_vaf += 1
                continue

        if clinvar_keep:
            def field(name):
                i = col.get(name)
                return fields[i] if i is not None and i < len(fields) else ''
            if not clinvar_kept(field(CLINVAR_COLUMN), field(CLINVAR_CONFLICT_COLUMN), clinvar_keep):
                removed_clinvar += 1
                continue

        kept.append(line)

    with open(args.out, 'w') as out:
        out.write('\n'.join(kept) + '\n')

    total = len(lines) - 1
    print(f"filter_mutations_af: {total} rows in, {len(kept) - 1} kept, "
          f"{removed_gnomad} removed by gnomAD AF > {args.max_gnomad_af}, "
          f"{removed_vaf} removed by VAF > {args.max_vaf}, "
          f"{removed_clinvar} removed by ClinVar (not {args.clinvar_keep})", file=sys.stderr)


if __name__ == '__main__':
    main()
