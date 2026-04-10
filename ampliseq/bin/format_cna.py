#!/usr/bin/env python3
import pandas as pd
import os
import sys

def copy_number_to_value(cn):
    if cn == 0:
        return -2
    elif cn == 1:
        return -1
    elif cn == 3:
        return 1
    elif cn >= 4:
        return 2
    return None

def main():
    if len(sys.argv) != 3:
        print("Usage: python format_cna.py <input.tsv> <Sample_Id>")
        sys.exit(1)

    input_file, sample_id = sys.argv[1], sys.argv[2]

    df = pd.read_csv(input_file, sep='\t')

    filtered = df[df['Variant Subtype'].isin(['DUPLICATION', 'DELETION'])]

    rows = []
    for _, row in filtered.iterrows():
        try:
            cn = float(row['Copy Number'])
        except (ValueError, TypeError):
            continue
        value = copy_number_to_value(cn)
        if value is not None:
            rows.append({'Hugo_Symbol': row['Genes'], 'Sample_Id': sample_id, 'Value': value})

    out = pd.DataFrame(rows, columns=['Hugo_Symbol', 'Sample_Id', 'Value'])
    write_header = not os.path.isfile('data_cna.txt') or os.path.getsize('data_cna.txt') == 0
    out.to_csv('data_cna.txt', sep='\t', index=False, mode='a', header=write_header)

if __name__ == '__main__':
    main()
