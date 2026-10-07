#!/usr/bin/env python3
import pandas as pd
import sys

def main():
    if len(sys.argv) != 3:
        print("Usage: python format_cna_deanon.py <data_cna.txt> <linking_file.txt>")
        sys.exit(1)

    cna_file, linking_file = sys.argv[1], sys.argv[2]

    linking = pd.read_csv(linking_file, sep='\t', header=0, usecols=[0, 1], dtype=str, keep_default_na=False)
    linking.columns = ['Anon_Id', 'Real_Id']
    # Uppercase keys for case-insensitive matching
    id_map = {k.upper(): v for k, v in zip(linking['Anon_Id'], linking['Real_Id'])}

    df = pd.read_csv(cna_file, sep='\t', dtype=str, keep_default_na=False)

    unmatched = sorted(set(df['Sample_Id'].str.upper()) - set(id_map))
    for uid in unmatched:
        print(f"WARNING: no linking entry for Sample_Id '{uid}', leaving unchanged", file=sys.stderr)

    df['Sample_Id'] = df['Sample_Id'].str.upper().map(id_map).fillna(df['Sample_Id'])
    df.to_csv(cna_file, sep='\t', index=False)

if __name__ == '__main__':
    main()
