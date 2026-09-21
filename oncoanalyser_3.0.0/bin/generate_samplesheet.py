#!/usr/bin/env python3
"""Generate an nf-core CSV samplesheet from a directory of sample folders."""

import argparse
import csv
import os
import sys


def parse_folder_name(name):
    """Split folder name (minus trailing _WGTS) into group and subject_id."""
    stem = name.removesuffix("_WGTS")
    parts = stem.rsplit("-", 1)
    if len(parts) != 2:
        sys.exit(f"ERROR: cannot parse folder name '{name}' into group-number format")
    group, _ = parts
    subject_id = stem
    return group, subject_id


def main():
    parser = argparse.ArgumentParser(description="Create nf-core samplesheet CSV.")
    parser.add_argument("input_dir", help="Directory containing *_WGTS/ sample folders")
    parser.add_argument("-o", "--output", default="samplesheet.csv", help="Output CSV path (default: samplesheet.csv)")
    args = parser.parse_args()

    input_dir = os.path.abspath(args.input_dir)
    folders = sorted(
        d for d in os.listdir(input_dir)
        if os.path.isdir(os.path.join(input_dir, d)) and d.endswith("_WGTS")
    )

    if not folders:
        sys.exit(f"ERROR: no *_WGTS/ folders found in {input_dir}")

    with open(args.output, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["group", "subject_id", "sample_id", "folder"])
        for folder in folders:
            group, subject_id = parse_folder_name(folder)
            sample_id = f"{subject_id}-T"
            writer.writerow([group, subject_id, sample_id, os.path.join(input_dir, folder)])

    print(f"Wrote {len(folders)} samples to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
