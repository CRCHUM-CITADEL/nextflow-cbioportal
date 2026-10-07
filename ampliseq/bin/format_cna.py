#!/usr/bin/env python3
"""
Format DUPLICATION/DELETION rows of analysis_<number>_export.tsv into cBioPortal discrete long
CNA format (Hugo_Symbol, Sample_Id, Value).
Usage: python format_cna.py <input.tsv> <Sample_Id>

Copy number -> discrete value, after rounding to the nearest integer copy (half-up):
0 -> -2 (homozygous deletion), 1 -> -1, 2 -> neutral (dropped), 3 -> 1, >=4 -> 2.
A `Genes` cell listing several genes (comma/semicolon separated) gives one row per gene, and a
gene seen twice in a sample keeps its most extreme value (cBioPortal rejects duplicate rows).
"""
import math
import os
import re
import sys

import pandas as pd


def copy_number_to_value(cn):
    cn = math.floor(cn + 0.5)
    if cn <= 0:
        return -2
    if cn == 1:
        return -1
    if cn == 3:
        return 1
    if cn >= 4:
        return 2
    return None


def main():
    if len(sys.argv) != 3:
        print("Usage: python format_cna.py <input.tsv> <Sample_Id>")
        sys.exit(1)

    input_file, sample_id = sys.argv[1], sys.argv[2]

    df = pd.read_csv(input_file, sep='\t', dtype=str, keep_default_na=False)
    filtered = df[df['Variant Subtype'].isin(['DUPLICATION', 'DELETION'])]

    best = {}  # gene -> value, most extreme wins
    for _, row in filtered.iterrows():
        try:
            cn = float(row['Copy Number'])
        except (ValueError, TypeError):
            print(f"WARNING: unparsable Copy Number {row['Copy Number']!r} for {row['Genes']!r}, skipping",
                  file=sys.stderr)
            continue
        value = copy_number_to_value(cn)
        if value is None:
            continue
        for gene in filter(None, (g.strip() for g in re.split(r'[,;]', row['Genes']))):
            if gene not in best or abs(value) > abs(best[gene]):
                best[gene] = value

    out = pd.DataFrame([{'Hugo_Symbol': g, 'Sample_Id': sample_id, 'Value': v} for g, v in best.items()],
                       columns=['Hugo_Symbol', 'Sample_Id', 'Value'])
    write_header = not os.path.isfile('data_cna.txt') or os.path.getsize('data_cna.txt') == 0
    out.to_csv('data_cna.txt', sep='\t', index=False, mode='a', header=write_header)


if __name__ == '__main__':
    main()
