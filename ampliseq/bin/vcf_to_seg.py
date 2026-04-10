#!/usr/bin/env python3
"""
Convert a CNV VCF file to cBioPortal .seg format.
Usage: vcf_to_seg.py <vcf_file> <sample_id>

Output: data_seg.txt (tab-separated, appended)
Columns: ID, chrom, loc.start, loc.end, num.mark, seg.mean

Only PASS records are included.
seg.mean = log2(CN/2); CN=0 yields -3.0 (homozygous deletion sentinel).
"""
import sys
import math
import os


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <vcf_file> <sample_id>", file=sys.stderr)
        sys.exit(1)

    vcf_file = sys.argv[1]
    sample_id = sys.argv[2]
    out_file = "data_seg.txt"

    rows = []
    with open(vcf_file) as fh:
        format_idx = None
        sample_idx = None
        for line in fh:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                cols = line.strip().split('\t')
                format_idx = cols.index('FORMAT')
                sample_idx = cols.index('SAMPLE')
                continue

            if format_idx is None:
                continue

            cols = line.rstrip('\n').split('\t')
            if len(cols) <= sample_idx:
                continue

            chrom = cols[0]
            pos = int(cols[1])
            filter_val = cols[6]
            info = cols[7]
            fmt = cols[format_idx]
            sample = cols[sample_idx]

            if filter_val != 'PASS':
                continue

            # remove 'chr' prefix
            chrom = chrom[3:] if chrom.startswith('chr') else chrom

            # Parse END from INFO; fall back to POS for point variants
            end = pos
            if info != '.':
                for field in info.split(';'):
                    if field.startswith('END='):
                        end = int(field.split('=', 1)[1])
                        break

            # Parse CN value from FORMAT/sample columns
            try:
                cn = int(sample)
            except (IndexError, ValueError):
                print(f"WARNING: could not parse CN at {chrom}:{pos}, skipping", file=sys.stderr)
                continue

            ## Usually seg mean is this : 
            # seg.mean = log2(CN/2); CN=0 → homozygous deletion sentinel
            # but we will just put the copy number for ampliseq data
            seg_mean = cn
            if seg_mean != 2:
                rows.append((sample_id, chrom, pos, end, 1, seg_mean))

    write_header = not os.path.exists(out_file)
    with open(out_file, 'a') as out:
        if write_header:
            out.write('ID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean\n')
        for row in rows:
            out.write('\t'.join(str(x) for x in row) + '\n')


if __name__ == '__main__':
    main()
