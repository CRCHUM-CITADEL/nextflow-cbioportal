#!/usr/bin/env python3
"""
Transform deanonymised clinical sample data to cBioPortal sample input format.
Usage: python3 format_clinical_sample.py <input_file>
"""

import sys
import os
import pandas as pd

HEADER_LINES = [
    "#Patient Identifier\tSample Identifier\tCancer Type\tCancer Type Detailed\tSample Type\tPrimary site",
    "#Patient identifier\tSample Identifier\tCancer type\tSub-type of the specified\tType of cancer sample\tPrimary site of cancer sample",
    "#STRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING",
    "#1\t1\t1\t1\t1\t1",
    "PATIENT_ID\tSAMPLE_ID\tCANCER_TYPE\tCANCER_TYPE_DETAILED\tSAMPLE_TYPE\tPRIMARY_SITE",
]

def transform_value(val):
    """Return NA for Unknown/empty, else the original value."""
    if pd.isna(val) or str(val).strip().lower() == "unknown":
        return "NA"
    return val

def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <input_file>", file=sys.stderr)
        sys.exit(1)

    df = pd.read_csv(sys.argv[1], sep="\t", dtype=str, usecols=range(8))

    # Drop duplicate sample_id column (columns 0 and 1 are both sample_id)
    df.columns = ["sample_id_drop", "sample_id", "patient_id", "cancer_type",
                  "cancer_type_detailed", "sample_type", "tumor_site", "tumor_purity"]
    df = df.drop(columns=["sample_id_drop", "tumor_purity"])

    out = pd.DataFrame({
        "PATIENT_ID":            df["patient_id"],
        "SAMPLE_ID":             df["sample_id"],
        "CANCER_TYPE":           df["cancer_type"].apply(transform_value),
        "CANCER_TYPE_DETAILED":  df["cancer_type_detailed"].apply(transform_value),
        "SAMPLE_TYPE":           df["sample_type"].apply(transform_value),
        "PRIMARY_SITE":          df["tumor_site"].apply(transform_value),
    })

    out_path = os.path.join(os.getcwd(), "data_clinical_sample.txt")
    with open(out_path, "w") as f:
        f.write("\n".join(HEADER_LINES) + "\n")
        out.to_csv(f, sep="\t", index=False, header=False)

    print(f"Written: {out_path} ({len(out)} samples)")

if __name__ == "__main__":
    main()