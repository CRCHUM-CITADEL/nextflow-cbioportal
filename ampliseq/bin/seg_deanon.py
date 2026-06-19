#!/usr/bin/env python3
"""
Deanonymize the ID column in a cBioPortal .seg file using a linking file.
Usage: seg_deanon.py <data_seg.txt> <linking_file.txt>

Edits the file in-place. Unmatched IDs are left unchanged with a warning.
"""
import sys
import pandas as pd


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <data_seg.txt> <linking_file.txt>", file=sys.stderr)
        sys.exit(1)

    seg_file, linking_file = sys.argv[1], sys.argv[2]

    linking = pd.read_csv(linking_file, sep='\t', header=0, usecols=[0, 1])
    linking.columns = ['Anon_Id', 'Real_Id']
    # Uppercase keys for case-insensitive matching
    id_map = {k.upper(): v for k, v in zip(linking['Anon_Id'], linking['Real_Id'])}

    df = pd.read_csv(seg_file, sep='\t')

    unmatched = set(df['ID'].str.upper()) - set(id_map)
    for uid in sorted(unmatched):
        print(f"WARNING: no linking entry for ID '{uid}', leaving unchanged", file=sys.stderr)

    df['ID'] = df['ID'].str.upper().map(id_map).fillna(df['ID'])
    df.to_csv(seg_file, sep='\t', index=False)


if __name__ == '__main__':
    main()
