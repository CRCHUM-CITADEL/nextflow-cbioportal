#!/usr/bin/env python3
"""
Generate a samplesheet.csv for the ampliseq-cbioportal pipeline.

Usage:
    python3 generate_samplesheet.py <input_dir> [--output samplesheet.csv]

Arguments:
    input_dir   Root directory containing one or more cohort subdirectories.
                Each cohort subdirectory contains per-sample folders.

The cohort subdirectory name becomes the `group` column.
The subject_id is derived from the sample folder name (part after the first `_`).
The sample_id is the prefix before `-basespace-` in the pisces VCF filename.

Required files per sample folder:
    - *-basespace-pisces.final.vcf.gz
    - analysis_*_export.tsv
"""

import argparse
import csv
import re
import sys
from pathlib import Path


PISCES_PATTERN = re.compile(r"^(.+)-basespace-pisces\.final\.vcf\.gz$", re.IGNORECASE)
TSV_PATTERN = re.compile(r"^analysis_.+_export\.tsv$", re.IGNORECASE)


def find_pisces_vcf(sample_dir: Path) -> tuple[str, Path] | None:
    """Return (sample_id, vcf_path) from *-basespace-pisces.final.vcf.gz, or None."""
    for f in sample_dir.iterdir():
        m = PISCES_PATTERN.match(f.name)
        if m:
            return m.group(1), f
    return None


def find_tsv(sample_dir: Path) -> Path | None:
    """Return path to analysis_*_export.tsv, or None."""
    for f in sample_dir.iterdir():
        if TSV_PATTERN.match(f.name):
            return f
    return None


def subject_id_from_sample_id(sample_id: str) -> str:
    """Extract subject ID from sample_id (part before the first `_`).

    e.g. `25bm002356_5924` -> `25bm002356`
    """
    return sample_id.split("_", 1)[0]


def is_sample_dir(path: Path) -> bool:
    """Return True if the directory looks like a sample folder (contains a pisces VCF)."""
    return path.is_dir() and any(PISCES_PATTERN.match(f.name) for f in path.iterdir())


def process_cohort(cohort_dir: Path, group: str, writer: csv.DictWriter) -> int:
    """Process all sample subdirectories within a cohort directory. Returns number of rows written."""
    rows_written = 0
    sample_dirs = sorted(d for d in cohort_dir.iterdir() if is_sample_dir(d))

    if not sample_dirs:
        print(f"  WARNING: no sample directories found in '{cohort_dir.name}'", file=sys.stderr)
        return 0

    for sample_dir in sample_dirs:
        vcf_result = find_pisces_vcf(sample_dir)
        tsv_path = find_tsv(sample_dir)
        missing = []

        if vcf_result is None:
            missing.append("*-basespace-pisces.final.vcf.gz")
        if tsv_path is None:
            missing.append("analysis_*_export.tsv")

        if missing:
            print(
                f"  SKIP {sample_dir.name}: missing {', '.join(missing)}",
                file=sys.stderr,
            )
            continue

        sample_id, _ = vcf_result
        subject_id = subject_id_from_sample_id(sample_id)

        writer.writerow(
            {
                "group": group,
                "subject_id": subject_id,
                "sample_id": sample_id,
                "folder_location": str(sample_dir),
            }
        )
        rows_written += 1

    return rows_written


def main():
    parser = argparse.ArgumentParser(
        description="Generate samplesheet.csv for the ampliseq-cbioportal pipeline."
    )
    parser.add_argument("input_dir", help="Root directory containing cohort subdirectories")
    parser.add_argument(
        "--output", "-o", default="samplesheet.csv", help="Output CSV file (default: samplesheet.csv)"
    )
    args = parser.parse_args()

    input_dir = Path(args.input_dir).resolve()
    if not input_dir.is_dir():
        print(f"ERROR: '{input_dir}' is not a directory", file=sys.stderr)
        sys.exit(1)

    output_path = Path(args.output)
    fieldnames = ["group", "subject_id", "sample_id", "folder_location"]
    total_rows = 0

    with output_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()

        # Each immediate subdirectory that is NOT itself a sample dir is treated as a cohort dir.
        cohort_dirs = sorted(
            d for d in input_dir.iterdir() if d.is_dir() and not is_sample_dir(d)
        )
        # If the input_dir itself contains sample dirs directly, treat it as the cohort.
        direct_samples = [d for d in input_dir.iterdir() if is_sample_dir(d)]

        if cohort_dirs:
            for cohort_dir in cohort_dirs:
                group = cohort_dir.name
                print(f"Processing cohort: {group}")
                n = process_cohort(cohort_dir, group, writer)
                print(f"  -> {n} sample(s) added")
                total_rows += n
        elif direct_samples:
            group = input_dir.name
            print(f"Processing cohort: {group} (input_dir used as cohort)")
            n = process_cohort(input_dir, group, writer)
            print(f"  -> {n} sample(s) added")
            total_rows += n
        else:
            print(f"WARNING: no cohort or sample directories found under '{input_dir}'", file=sys.stderr)

    print(f"\nWrote {total_rows} row(s) to '{output_path}'")


if __name__ == "__main__":
    main()
