#!/usr/bin/env python3
"""
Replace anonymised sample IDs in Sample_Id column of a TSV
using a linking file mapping sample_id -> deanon_sample_id.
Usage: python3 deanon_sample_ids.py <input_tsv> <linking_file>
"""

import sys
import os
import pandas as pd

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input_tsv> <linking_file>", file=sys.stderr)
        sys.exit(1)

    input_tsv, linking_file = sys.argv[1], sys.argv[2]

    df      = pd.read_csv(input_tsv, sep="\t", dtype=str)
    linking = pd.read_csv(linking_file, sep="\t", dtype=str)

    mapping = linking.set_index("sample_id")["deanon_sample_id"].to_dict()

    # Warn on any unmatched IDs
    unmatched = df[~df["Sample_Id"].isin(mapping)]["Sample_Id"].unique()
    if len(unmatched):
        print(f"WARNING: {len(unmatched)} unmatched Sample_Id value(s) left as-is:",
              file=sys.stderr)
        for uid in unmatched:
            print(f"  {uid}", file=sys.stderr)

    df["Sample_Id"] = df["Sample_Id"].map(mapping).fillna(
        df["Sample_Id"]
    )

    out_path = os.path.join(os.getcwd(), os.path.basename(input_tsv))
    df.to_csv(out_path, sep="\t", index=False)
    print(f"Written: {out_path}")

if __name__ == "__main__":
    main()