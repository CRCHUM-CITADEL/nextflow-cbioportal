#!/usr/bin/env python3
"""Generate a template cBioPortal data_clinical_sample.txt for a genomic-only study.

Only values the pipeline actually knows are written:
  PATIENT_ID     subject_id from the genomic samplesheet
  SAMPLE_ID      sample_id from the genomic samplesheet
  ONCOTREE_CODE  --cancer-type (the cancer_type param), upper-cased

Every other column is written empty. Nothing is inferred: the empty columns
must be filled in by hand before the study is imported into cBioPortal.

With --existing (the data_clinical_sample.txt already in the output folder, e.g.
filled in by hand or written by --mode clinical), that file is updated instead of
replaced: its columns and the values of every sample still in the samplesheet are
kept, samples new to the samplesheet get a row with only the known values, and rows
for samples no longer in the samplesheet are dropped (the samplesheet defines the
cohort). A rerun therefore never wipes values someone filled in.

Format: https://docs.cbioportal.org/file-formats/#clinical-data
"""

import argparse
import csv
import re
import sys

# (attribute name, display name, description, datatype, priority)
COLUMNS = [
    ("PATIENT_ID", "Patient Identifier", "Identifier to uniquely specify a patient.", "STRING", "1"),
    ("SAMPLE_ID", "Sample Identifier", "A unique sample identifier.", "STRING", "1"),
    ("ONCOTREE_CODE", "Oncotree Code", "Oncotree Code", "STRING", "1"),
    ("CANCER_TYPE", "Cancer Type", "Cancer Type", "STRING", "1"),
    ("CANCER_TYPE_DETAILED", "Cancer Type Detailed", "Cancer Type Detailed", "STRING", "1"),
    ("SAMPLE_TYPE", "Sample Type", "The type of sample (i.e., normal, primary, met, recurrence).", "STRING", "1"),
]

# cBioPortal: IDs may only contain numbers, letters, points, underscores and hyphens
INVALID_ID_CHARACTERS = re.compile(r"[^A-Za-z0-9._-]")


def read_linking(path):
    """Return (subject_id, sample_id) pairs from the linking TSV, sorted by sample_id."""
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = {"subject_id", "sample_id"} - set(reader.fieldnames or [])
        if missing:
            sys.exit(f"ERROR: {path} is missing column(s): {', '.join(sorted(missing))}")
        pairs = [((row["subject_id"] or "").strip(), (row["sample_id"] or "").strip()) for row in reader]

    if not pairs:
        sys.exit(f"ERROR: no samples found in {path}")

    seen = set()
    for subject, sample in pairs:
        if not subject or not sample:
            sys.exit(f"ERROR: blank subject_id or sample_id in {path}: '{subject}' / '{sample}'")
        if sample in seen:
            sys.exit(f"ERROR: sample_id defined twice in {path}: {sample}")
        seen.add(sample)
        for value in (subject, sample):
            if INVALID_ID_CHARACTERS.search(value):
                print(
                    f"WARNING: '{value}' contains characters cBioPortal does not allow in IDs "
                    "(only letters, numbers, points, underscores and hyphens)",
                    file=sys.stderr,
                )

    return sorted(pairs, key=lambda pair: pair[1])


def read_existing(path):
    """Return (comment_rows, header, rows) of an existing clinical sample file.

    comment_rows are the leading '#' lines, kept verbatim; rows are dicts keyed by header.
    """
    with open(path, newline="") as fh:
        lines = [line.rstrip("\r\n") for line in fh]

    comment_rows = []
    while lines and lines[0].startswith("#"):
        comment_rows.append(lines.pop(0))
    if not lines:
        sys.exit(f"ERROR: {path} has no header row")

    header = lines.pop(0).split("\t")
    missing = {"PATIENT_ID", "SAMPLE_ID"} - set(header)
    if missing:
        sys.exit(f"ERROR: existing {path} is missing column(s): {', '.join(sorted(missing))}")

    rows = []
    for number, line in enumerate(lines, start=len(comment_rows) + 2):
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) > len(header):
            sys.exit(f"ERROR: existing {path} line {number} has more fields than the header")
        rows.append(dict(zip(header, fields + [""] * (len(header) - len(fields)))))
    return comment_rows, header, rows


def write_table(path, comment_rows, header, rows):
    with open(path, "w", newline="") as fh:
        for line in comment_rows:
            fh.write(line + "\n")
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join(row.get(col, "") for col in header) + "\n")


def main():
    parser = argparse.ArgumentParser(description="Create a template cBioPortal data_clinical_sample.txt.")
    parser.add_argument("--linking", required=True, help="TSV with subject_id and sample_id columns")
    parser.add_argument("--cancer-type", default="", help="OncoTree code from the cancer_type param (may be empty)")
    parser.add_argument("--existing", default=None, help="Existing data_clinical_sample.txt to update instead of replace (optional)")
    parser.add_argument("-o", "--output", default="data_clinical_sample.txt", help="Output path (default: data_clinical_sample.txt)")
    args = parser.parse_args()

    pairs = read_linking(args.linking)
    oncotree_code = args.cancer_type.strip().upper()

    def template_row(subject, sample):
        return {"PATIENT_ID": subject, "SAMPLE_ID": sample, "ONCOTREE_CODE": oncotree_code}

    if args.existing:
        comment_rows, header, old_rows = read_existing(args.existing)
        by_sample = {row["SAMPLE_ID"]: row for row in old_rows}
        rows = [by_sample.get(sample) or template_row(subject, sample) for subject, sample in pairs]
        kept = sum(1 for _subject, sample in pairs if sample in by_sample)
        dropped = sorted(set(by_sample) - {sample for _subject, sample in pairs})
        write_table(args.output, comment_rows, header, rows)
        print(f"Updated {args.existing}: kept {kept} existing row(s), added {len(pairs) - kept} new sample(s)", file=sys.stderr)
        if dropped:
            print(f"WARNING: dropped {len(dropped)} row(s) for samples not in the samplesheet: {', '.join(dropped)}", file=sys.stderr)
        return

    # Four '#' rows describing the attributes, then the attribute names
    comment_rows = ["#" + "\t".join(col[i] for col in COLUMNS) for i in (1, 2, 3, 4)]
    header = [col[0] for col in COLUMNS]
    write_table(args.output, comment_rows, header, [template_row(subject, sample) for subject, sample in pairs])

    empty = [col[0] for col in COLUMNS if col[0] not in ("PATIENT_ID", "SAMPLE_ID")
             and not (col[0] == "ONCOTREE_CODE" and oncotree_code)]
    print(f"Wrote {len(pairs)} samples to {args.output}", file=sys.stderr)
    print(f"TEMPLATE: columns left empty, to be filled in by hand: {', '.join(empty)}", file=sys.stderr)


if __name__ == "__main__":
    main()
