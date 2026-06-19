#!/usr/bin/env python3
"""
Build an anonymised linking file from the original linking file.

Replaces deanon_sample_id and deanon_patient_id with the lowercase of
sample_id, so the existing deanon scripts remap genomic IDs to anonymised
form instead of to real (deanonymised) IDs.

Usage: build_anon_linking.py <linking_file>

Writes anon_linking.txt to the current directory.
"""

import sys
import os
import pandas as pd


def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <linking_file>", file=sys.stderr)
        sys.exit(1)

    linking = pd.read_csv(sys.argv[1], sep="\t", dtype=str)

    linking["deanon_sample_id"] = linking["sample_id"].str.lower()
    linking["deanon_patient_id"] = linking["sample_id"].str.lower()

    out_path = os.path.join(os.getcwd(), "anon_linking.txt")
    linking.to_csv(out_path, sep="\t", index=False)
    print(f"Written: {out_path}")


if __name__ == "__main__":
    main()
