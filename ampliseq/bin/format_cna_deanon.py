#!/usr/bin/env python3
import pandas as pd
import sys

def main():
    if len(sys.argv) != 3:
        print("Usage: python format_cna_deanon.py <data_cna.txt> <linking_file.txt>")
        sys.exit(1)

    cna_file, linking_file = sys.argv[1], sys.argv[2]

    linking = pd.read_csv(linking_file, sep='\t', header=0, usecols=[0, 1])
    linking.columns = ['Anon_Id', 'Real_Id']
    id_map = dict(zip(linking['Anon_Id'], linking['Real_Id']))

    df = pd.read_csv(cna_file, sep='\t')
    df['Sample_Id'] = df['Sample_Id'].map(id_map).fillna(df['Sample_Id'])
    df.to_csv(cna_file, sep='\t', index=False)

if __name__ == '__main__':
    main()