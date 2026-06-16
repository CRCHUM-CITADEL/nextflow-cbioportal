#!/usr/bin/env python3
"""
Replace deanonymised (real) IDs with anonymised IDs in a cBioPortal clinical file.

Uses the original linking file (real IDs) and the anonymised linking file to build
reverse mappings:
  real sample_id  -> anonymised sample_id
  real patient_id -> anonymised patient_id

Usage: anonymize_clinical.py <clinical_file> <orig_linking> <anon_linking>

The script preserves the cBioPortal header lines (starting with '#') and
rewrites the file in-place.
"""

import sys
import pandas as pd


def main():
    if len(sys.argv) != 4:
        print(
            f"Usage: {sys.argv[0]} <clinical_file> <orig_linking> <anon_linking>",
            file=sys.stderr,
        )
        sys.exit(1)

    clinical_file = sys.argv[1]
    orig_linking = pd.read_csv(sys.argv[2], sep="\t", dtype=str)
    anon_linking = pd.read_csv(sys.argv[3], sep="\t", dtype=str)

    # Map real sample IDs -> anonymised sample IDs
    sample_map = dict(zip(orig_linking["deanon_sample_id"], anon_linking["deanon_sample_id"]))
    # Map real patient IDs -> anonymised patient IDs
    patient_map = dict(zip(orig_linking["deanon_patient_id"], anon_linking["deanon_patient_id"]))

    # Read header lines (starting with '#') separately
    header_lines = []
    with open(clinical_file) as f:
        for line in f:
            if line.startswith("#"):
                header_lines.append(line)
            else:
                break

    df = pd.read_csv(clinical_file, sep="\t", comment="#", dtype=str)

    if "SAMPLE_ID" in df.columns:
        df["SAMPLE_ID"] = df["SAMPLE_ID"].map(sample_map).fillna(df["SAMPLE_ID"])

    if "PATIENT_ID" in df.columns:
        df["PATIENT_ID"] = df["PATIENT_ID"].map(patient_map).fillna(df["PATIENT_ID"])

    with open(clinical_file, "w") as f:
        for line in header_lines:
            f.write(line)
        df.to_csv(f, sep="\t", index=False, header=True)


if __name__ == "__main__":
    main()
