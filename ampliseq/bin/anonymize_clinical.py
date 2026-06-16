#!/usr/bin/env python3
"""
Replace deanonymised (real) IDs with anonymised IDs in a cBioPortal clinical file.

Uses the linking file to build reverse mappings:
  deanon_sample_id  -> sample_id     (for SAMPLE_ID column)
  deanon_patient_id -> subject_id    (for PATIENT_ID column, from patient_map)

Usage: anonymize_clinical.py <clinical_file> <linking_file> <patient_map_file>

patient_map_file is a two-column TSV (no header): subject_id<TAB>sample_id
used to derive the anonymised patient IDs from the samplesheet.

The script preserves the cBioPortal header lines (starting with '#') and
rewrites the file in-place.
"""

import sys
import pandas as pd


def main():
    if len(sys.argv) != 4:
        print(
            f"Usage: {sys.argv[0]} <clinical_file> <linking_file> <patient_map_file>",
            file=sys.stderr,
        )
        sys.exit(1)

    clinical_file = sys.argv[1]
    linking_file = sys.argv[2]
    patient_map_file = sys.argv[3]

    linking = pd.read_csv(linking_file, sep="\t", dtype=str)
    patient_map_df = pd.read_csv(
        patient_map_file, sep="\t", dtype=str, header=None, names=["subject_id", "sample_id"]
    )

    # sample mapping: deanon_sample_id -> sample_id (anonymised)
    sample_map = dict(zip(linking["deanon_sample_id"], linking["sample_id"]))

    # patient mapping: deanon_patient_id -> subject_id (anonymised)
    # Join linking with patient_map via sample_id to get subject_id for each deanon_patient_id
    merged = linking.merge(patient_map_df, on="sample_id", how="left")
    patient_map = dict(zip(merged["deanon_patient_id"], merged["subject_id"]))

    # Read header lines (starting with '#') and data separately
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
