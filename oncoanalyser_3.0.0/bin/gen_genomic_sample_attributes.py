#!/usr/bin/env python3
"""Write one sample's cBioPortal row of genomic sample attributes.

TMB/TML come from PURPLE's purity.tsv, HRD from CHORD's prediction.tsv. Columns are
found by name. A missing file gives NA for its columns; a present file without the
expected columns is an error, so a format change cannot silently blank the study.

Output: the 4 cBioPortal '#' header lines, the attribute-ID line and one data row,
so per-sample files merge with collectFile(keepHeader: true, skip: 5).
"""

import argparse
import csv
import sys

# (cBioPortal ID, display name, description, datatype)
COLUMNS = [
    ("PATIENT_ID", "Patient Identifier", "Identifier to uniquely specify a patient.", "STRING"),
    ("SAMPLE_ID", "Sample Identifier", "A unique sample identifier.", "STRING"),
    ("TMB_PER_MB", "TMB per Mb (PURPLE)",
     "Tumour mutational burden: passing variants per Mb (PURPLE tmbPerMb).", "NUMBER"),
    ("TMB_STATUS", "TMB Status (PURPLE)", "PURPLE TMB status: HIGH (> 10 per Mb), LOW or UNKNOWN.", "STRING"),
    ("TML", "TML (PURPLE)", "Tumour mutational load: missense variants (PURPLE tml).", "NUMBER"),
    ("TML_STATUS", "TML Status (PURPLE)", "PURPLE TML status: HIGH, LOW or UNKNOWN.", "STRING"),
    ("HRD_SCORE", "HRD Score (CHORD)", "CHORD probability of homologous recombination deficiency.", "NUMBER"),
    ("HR_STATUS", "HR Status (CHORD)", "CHORD homologous recombination status.", "STRING"),
    ("HRD_TYPE", "HRD Type (CHORD)", "CHORD HRD type: BRCA1_type, BRCA2_type or none.", "STRING"),
]

# cBioPortal ID -> candidate source columns, first match wins
PURITY_FIELDS = {"TMB_PER_MB": ["tmbPerMb"], "TMB_STATUS": ["tmbStatus"],
                 "TML": ["tml"], "TML_STATUS": ["tmlStatus"]}
CHORD_FIELDS = {"HRD_SCORE": ["hrd", "p_hrd"], "HR_STATUS": ["hr_status"], "HRD_TYPE": ["hrd_type"]}


def read_first_row(path, fields, label):
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        sys.exit(f"ERROR: {label} file {path} has no data row")
    row = rows[0]
    out = {}
    for cbio_id, candidates in fields.items():
        source = next((c for c in candidates if c in row), None)
        if source is None:
            sys.exit(f"ERROR: {label} file {path} has none of the columns {', '.join(candidates)}")
        value = (row[source] or "").strip()
        out[cbio_id] = value if value else "NA"
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--patient", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--purity", help="PURPLE <tumor>.purple.purity.tsv (optional)")
    p.add_argument("--chord", help="CHORD <tumor>.chord.prediction.tsv (optional)")
    p.add_argument("--output", required=True)
    args = p.parse_args()

    values = {cbio_id: "NA" for cbio_id, *_ in COLUMNS}
    values["PATIENT_ID"] = args.patient
    values["SAMPLE_ID"] = args.sample
    if args.purity:
        values.update(read_first_row(args.purity, PURITY_FIELDS, "PURPLE purity"))
    if args.chord:
        values.update(read_first_row(args.chord, CHORD_FIELDS, "CHORD prediction"))

    with open(args.output, "w") as out:
        out.write("#" + "\t".join(c[1] for c in COLUMNS) + "\n")
        out.write("#" + "\t".join(c[2] for c in COLUMNS) + "\n")
        out.write("#" + "\t".join(c[3] for c in COLUMNS) + "\n")
        out.write("#" + "\t".join("1" for _ in COLUMNS) + "\n")
        out.write("\t".join(c[0] for c in COLUMNS) + "\n")
        out.write("\t".join(values[c[0]] for c in COLUMNS) + "\n")


if __name__ == "__main__":
    main()
