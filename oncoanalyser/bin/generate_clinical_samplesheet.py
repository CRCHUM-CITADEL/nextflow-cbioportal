#!/usr/bin/env python3
"""Generate a clinical samplesheet CSV from a directory of ARGO/ICGC-ARGO clinical files."""

import argparse
import csv
import datetime
import os
import sys

KNOWN_FILETYPES = [
    "donors",
    "primary_diagnoses",
    "treatments",
    "surgeries",
    "systemic_therapies",
    "specimens",
    "radiations",
    "follow_ups",
    "sample_registrations",
    "biomarkers",
]


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Create a clinical samplesheet CSV from a directory of ARGO/ICGC-ARGO clinical CSVs.\n"
            "Files are matched by name stem (e.g. donors.csv → donors). "
            "Unrecognised CSV files are ignored with a warning."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Recognised filenames:\n"
            + "\n".join(f"  {ft}.csv" for ft in KNOWN_FILETYPES)
        ),
    )
    parser.add_argument(
        "input_dir",
        help="Directory containing clinical CSV files",
    )
    parser.add_argument(
        "-o", "--output",
        default="clinical_samplesheet.csv",
        help="Output CSV path (default: clinical_samplesheet.csv)",
    )
    parser.add_argument(
        "-d", "--date",
        default=None,
        help="Extraction date in YYYY-MM-DD format (default: today)",
    )
    parser.add_argument(
        "--study_id",
        default="",
        help="Study/group identifier written into the group_id column (default: empty)",
    )
    parser.add_argument(
        "--relative-to",
        default=None,
        metavar="BASE_DIR",
        help=(
            "If provided, file paths in the output are written relative to BASE_DIR "
            "instead of as absolute paths. Useful when the samplesheet will be used "
            "with a pipeline that resolves relative paths from its project directory."
        ),
    )
    args = parser.parse_args()

    extraction_date = args.date or datetime.date.today().isoformat()
    try:
        datetime.date.fromisoformat(extraction_date)
    except ValueError:
        sys.exit(f"ERROR: --date must be in YYYY-MM-DD format, got '{extraction_date}'")

    input_dir = os.path.abspath(args.input_dir)
    if not os.path.isdir(input_dir):
        sys.exit(f"ERROR: input directory does not exist: {input_dir}")

    base_dir = os.path.abspath(args.relative_to) if args.relative_to else None

    # Scan directory: match filenames to known filetypes
    found = {}
    for name in sorted(os.listdir(input_dir)):
        stem, ext = os.path.splitext(name)
        if ext.lower() != ".csv":
            continue
        if stem in KNOWN_FILETYPES:
            abs_path = os.path.join(input_dir, name)
            if base_dir:
                try:
                    path = os.path.relpath(abs_path, base_dir)
                except ValueError:
                    # On Windows, relpath can fail across drives — fall back to absolute
                    path = abs_path
            else:
                path = abs_path
            found[stem] = path
        else:
            print(f"WARNING: ignoring unrecognised file '{name}'", file=sys.stderr)

    if not found:
        sys.exit(f"ERROR: no recognised clinical CSV files found in {input_dir}")

    rows_written = 0
    with open(args.output, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["group_id", "filetype", "filepath", "extraction_date", "info"])
        for filetype in KNOWN_FILETYPES:
            if filetype in found:
                writer.writerow([args.study_id, filetype, found[filetype], extraction_date, ""])
                rows_written += 1

    print(f"Wrote {rows_written} row(s) to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
