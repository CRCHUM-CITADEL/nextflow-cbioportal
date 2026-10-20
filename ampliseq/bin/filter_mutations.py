#!/usr/bin/env python3
"""
Keep only the MAF mutations that overlap a good-quality mutation row of the ampliseq
analysis export TSV (used when --filter_tsv_variants is on).

Usage: filter_mutations.py <maf> <analysis_export.tsv> <out_file> <min_depth> <min_vaf>

A TSV row counts as a region when:
  - Start and End are whole numbers,
  - Depth and VAF are numbers (CNV / fusion rows have N/A there and never count),
  - Depth >= min_depth and VAF > min_vaf.

A MAF row is kept when its [Start_Position, End_Position] overlaps a region (inclusive) on
the same chromosome. Chromosome names are compared without a 'chr' prefix. Columns are found
by header name in both files; MAF lines are copied through unchanged.
"""
import re
import sys

NUMBER = re.compile(r'^[0-9]+([.][0-9]*)?([eE][-+]?[0-9]+)?$')


def no_chr(chrom):
    return chrom[3:] if chrom.startswith('chr') else chrom


def read_regions(tsv_file, min_depth, min_vaf):
    """Return {chromosome: [(start, end), ...]} for TSV rows that pass the QC thresholds."""
    regions = {}
    with open(tsv_file) as fh:
        header = [h.strip() for h in fh.readline().rstrip('\r\n').split('\t')]
        needed = ['Chr', 'Start', 'End', 'Depth', 'VAF']
        if not all(name in header for name in needed):
            sys.exit(f"ERROR: TSV {tsv_file} lacks one of Chr/Start/End/Depth/VAF "
                     "(needed by --filter_tsv_variants)")
        col = {name: header.index(name) for name in needed}

        for line in fh:
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) < len(header):
                fields += [''] * (len(header) - len(fields))
            start, end = fields[col['Start']], fields[col['End']]
            depth, vaf = fields[col['Depth']], fields[col['VAF']]

            if not (start.isdigit() and end.isdigit()):
                continue  # malformed row
            if not (NUMBER.match(depth) and NUMBER.match(vaf)):
                continue  # N/A: CNV / fusion row
            if float(depth) < min_depth or float(vaf) <= min_vaf:
                continue  # fails QC

            chrom = no_chr(fields[col['Chr']])
            regions.setdefault(chrom, []).append((int(start), int(end)))
    return regions


def main():
    if len(sys.argv) != 6:
        print(__doc__, file=sys.stderr)
        sys.exit(1)

    maf_file, tsv_file, out_file = sys.argv[1:4]
    min_depth, min_vaf = float(sys.argv[4]), float(sys.argv[5])

    regions = read_regions(tsv_file, min_depth, min_vaf)

    with open(maf_file) as fh, open(out_file, 'w') as out:
        header_line = fh.readline()
        header = header_line.rstrip('\r\n').split('\t')
        if not all(name in header for name in ('Chromosome', 'Start_Position', 'End_Position')):
            sys.exit("ERROR: MAF lacks Chromosome/Start_Position/End_Position")
        chrom_i = header.index('Chromosome')
        start_i = header.index('Start_Position')
        end_i = header.index('End_Position')
        out.write(header_line)

        total = kept = 0
        for line in fh:
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) <= max(chrom_i, start_i, end_i):
                continue
            total += 1
            try:
                start, end = int(fields[start_i]), int(fields[end_i])
            except ValueError:
                continue
            chrom_regions = regions.get(no_chr(fields[chrom_i]), [])
            if any(start <= r_end and end >= r_start for r_start, r_end in chrom_regions):
                out.write(line)
                kept += 1

    print(f"filter_mutations: {kept}/{total} MAF rows overlap a TSV mutation row "
          f"with Depth >= {sys.argv[4]} and VAF > {sys.argv[5]}", file=sys.stderr)


if __name__ == '__main__':
    main()
