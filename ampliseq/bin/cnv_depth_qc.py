#!/usr/bin/env python3
"""
Decide whether a sample's CNV calls are trustworthy, based on sequencing depth
of the variants in its Pisces VCF.
Usage: cnv_depth_qc.py <pisces.vcf[.gz]> <sample_id>
                       [--min-depth 250] [--min-fraction 0.95] [--pass-only]

A sample PASSes when at least --min-fraction of its variants have
depth >= --min-depth. Depth is read from FORMAT/DP, falling back to INFO/DP.
A sample with no countable variants FAILs.

Output (stdout, one line, tab-separated):
  sample_id  n_variants  n_ge_min_depth  fraction  PASS|FAIL
"""
import argparse
import gzip
import sys


def read_depth(cols, format_idx):
    if format_idx is not None and len(cols) > format_idx + 1:
        keys = cols[format_idx].split(':')
        values = cols[format_idx + 1].split(':')
        if 'DP' in keys:
            i = keys.index('DP')
            if i < len(values):
                try:
                    return int(values[i])
                except ValueError:
                    pass
    for field in cols[7].split(';'):
        if field.startswith('DP='):
            try:
                return int(field.split('=', 1)[1])
            except ValueError:
                return None
    return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('vcf')
    parser.add_argument('sample_id')
    parser.add_argument('--min-depth', type=int, default=250)
    parser.add_argument('--min-fraction', type=float, default=0.95)
    parser.add_argument('--pass-only', action='store_true')
    args = parser.parse_args()

    opener = gzip.open if args.vcf.endswith('.gz') else open

    n_variants = 0
    n_deep = 0
    format_idx = None
    with opener(args.vcf, 'rt') as fh:
        for line in fh:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                header = line.rstrip('\n').split('\t')
                format_idx = header.index('FORMAT') if 'FORMAT' in header else None
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 8:
                continue
            if args.pass_only and cols[6] != 'PASS':
                continue
            depth = read_depth(cols, format_idx)
            if depth is None:
                print(f"WARNING: no DP at {cols[0]}:{cols[1]} in {args.sample_id}, skipping",
                      file=sys.stderr)
                continue
            n_variants += 1
            if depth >= args.min_depth:
                n_deep += 1

    if n_variants == 0:
        print(f"WARNING: no countable variants for {args.sample_id}; marking FAIL",
              file=sys.stderr)
        fraction = 0.0
    else:
        fraction = n_deep / n_variants

    status = 'PASS' if n_variants and fraction >= args.min_fraction else 'FAIL'
    print(f"{args.sample_id}\t{n_variants}\t{n_deep}\t{fraction:.4f}\t{status}")


if __name__ == '__main__':
    main()
