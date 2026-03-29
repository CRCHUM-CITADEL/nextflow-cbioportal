#!/usr/bin/env python3
"""
Transform deanonymised clinical data to cBioPortal patient input format.
Usage: python3 format_clinical.py <input_file>
"""

import sys
import os
import pandas as pd

HEADER_LINES = [
    "#Patient Identifier\tAge\tSex\tOverall Survival Status\tOverall Survival (Months)\tSmoking History",
    "#Patient identifier\tAge at which the condition was first diagnosed, in years\tSex\tOverall Survival Status\tOverall survival in months since initial diagnosis\tSmoking history",
    "#STRING\tNUMBER\tSTRING\tSTRING\tNUMBER\tSTRING",
    "#1\t1\t1\t1\t1\t1",
    "PATIENT_ID\tAGE\tSEX\tOS_STATUS\tOS_MONTHS\tSMOKING_HISTORY",
]

OS_STATUS_MAP = {"0": "0:LIVING", "1": "1:DECEASED"}

def transform_value(val):
    """Return NA for Unknown/empty, else the original value."""
    if pd.isna(val) or str(val).strip().lower() == "unknown":
        return "NA"
    return val

def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <input_file>", file=sys.stderr)
        sys.exit(1)

    df = pd.read_csv(sys.argv[1], sep="\t", dtype=str)

    out = pd.DataFrame({
        "PATIENT_ID":       df["patient_id"],
        "AGE":              df["age"].apply(transform_value),
        "SEX":              df["sex"].apply(transform_value),
        "OS_STATUS":        df["os_status"].apply(
                                lambda v: OS_STATUS_MAP.get(str(v).strip(), "NA")
                                if str(v).strip().lower() != "unknown" else "NA"
                            ),
        "OS_MONTHS":        df["os_months"].apply(
                                lambda v: "NA" if transform_value(v) == "NA"
                                else str(v).replace(",", ".")
                            ),
        "SMOKING_HISTORY":  df["smoking_history"].apply(transform_value),
    })

    out_path = os.path.join(os.getcwd(), "data_clinical_patient.txt")
    with open(out_path, "w") as f:
        f.write("\n".join(HEADER_LINES) + "\n")
        out.to_csv(f, sep="\t", index=False, header=False)

    print(f"Written: {out_path} ({len(out)} patients)")

if __name__ == "__main__":
    main()