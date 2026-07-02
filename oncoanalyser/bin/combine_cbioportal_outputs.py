#!/usr/bin/env python3
"""Merge two cBioPortal output folders into one combined folder.

Usage:
    combine_cbioportal_outputs.py --input_dir_1 DIR1 --input_dir_2 DIR2 --output_dir OUT
    combine_cbioportal_outputs.py --input_dir_1 DIR1 --input_dir_2 DIR2 --output_dir OUT --study_id NEW_ID
"""

import argparse
import csv
import os
import shutil
import sys


# ---------------------------------------------------------------------------
# Dispatch tables
# ---------------------------------------------------------------------------

ROW_APPEND_FILES = [
    ("data_cna_hg38.seg", 1),
    ("data_cna_long.txt", 1),
    ("data_sv.txt", 1),
    ("util_linking_file.txt", 1),
]

MAF_FILE = "data_mutations_dna_rna_germline.txt"

WIDE_MATRIX_FILES = [
    ("data_expression.txt", ["Hugo_Symbol", "Entrez_Gene_Id"], "NA"),
    ("data_mutational_signatures_contribution_SBS.txt", ["ENTITY_STABLE_ID", "NAME", "DESCRIPTION"], "0"),
    ("data_mutational_signatures_counts_SBS.txt", ["ENTITY_STABLE_ID", "NAME"], "0"),
    ("data_mutational_signatures_contribution_DBS.txt", ["ENTITY_STABLE_ID", "NAME", "DESCRIPTION"], "0"),
    ("data_mutational_signatures_counts_DBS.txt", ["ENTITY_STABLE_ID", "NAME"], "0"),
    ("data_mutational_signatures_contribution_ID.txt", ["ENTITY_STABLE_ID", "NAME", "DESCRIPTION"], "0"),
    ("data_mutational_signatures_counts_ID.txt", ["ENTITY_STABLE_ID", "NAME"], "0"),
]

CLINICAL_FILES = [
    "data_clinical_sample.txt",
    "data_clinical_patient.txt",
]

META_FILES = [
    "meta_study.txt",
    "meta_cna_hg38.txt",
    "meta_cna_long.txt",
    "meta_sv.txt",
    "meta_expression.txt",
    "meta_sequenced.txt",
    "meta_clinical_sample.txt",
    "meta_clinical_patient.txt",
    "meta_mutational_signatures_contribution_SBS.txt",
    "meta_mutational_signatures_counts_SBS.txt",
    "meta_mutational_signatures_contribution_DBS.txt",
    "meta_mutational_signatures_counts_DBS.txt",
    "meta_mutational_signatures_contribution_ID.txt",
    "meta_mutational_signatures_counts_ID.txt",
]

SINGLE_VALUE_FILES = []

CASE_LISTS = ["cases_cnv.txt", "cases_sequenced.txt", "cases_sv.txt"]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def log(msg):
    """Print a log message to stderr."""
    print(f"[INFO] {msg}", file=sys.stderr)


def resolve(directory, filename):
    """Return the full path if the file exists, else None."""
    path = os.path.join(directory, filename)
    return path if os.path.isfile(path) else None


def parse_meta_file(path):
    """Parse a cBioPortal meta file (key: value format) into a dict."""
    meta = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            key, _, value = line.partition(":")
            meta[key.strip()] = value.strip()
    return meta


def collect_sample_ids(directory):
    """Extract sample IDs from a cBioPortal output folder."""
    linking = os.path.join(directory, "util_linking_file.txt")
    if os.path.isfile(linking):
        samples = set()
        with open(linking) as f:
            header = f.readline().strip().split("\t")
            idx = header.index("sample_id")
            for line in f:
                fields = line.strip().split("\t")
                if len(fields) > idx:
                    samples.add(fields[idx])
        return samples

    clinical = os.path.join(directory, "data_clinical_sample.txt")
    if os.path.isfile(clinical):
        samples = set()
        with open(clinical) as f:
            lines = f.readlines()
            # Skip comment lines (start with #) then read header
            data_start = 0
            for i, line in enumerate(lines):
                if not line.startswith("#"):
                    data_start = i
                    break
            header = lines[data_start].strip().split("\t")
            idx = header.index("SAMPLE_ID")
            for line in lines[data_start + 1:]:
                fields = line.strip().split("\t")
                if len(fields) > idx:
                    samples.add(fields[idx])
        return samples

    sys.exit(f"ERROR: cannot find util_linking_file.txt or data_clinical_sample.txt in {directory}")


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def validate_inputs(dir1, dir2, study_id_override=None):
    """Validate both directories and return the cancer_study_identifier to use.

    If study_id_override is provided, skip the mismatch check and return
    the override value. Otherwise require both dirs to share the same ID.
    """
    for d in (dir1, dir2):
        if not os.path.isdir(d):
            sys.exit(f"ERROR: directory does not exist: {d}")
        if not os.path.isfile(os.path.join(d, "meta_study.txt")):
            sys.exit(f"ERROR: meta_study.txt not found in {d}")

    if study_id_override:
        return study_id_override

    meta1 = parse_meta_file(os.path.join(dir1, "meta_study.txt"))
    meta2 = parse_meta_file(os.path.join(dir2, "meta_study.txt"))

    id1 = meta1.get("cancer_study_identifier", "")
    id2 = meta2.get("cancer_study_identifier", "")

    if id1 != id2:
        sys.exit(
            f"ERROR: cancer_study_identifier mismatch: '{id1}' (dir1) vs '{id2}' (dir2)"
        )

    return id1


def check_sample_overlap(samples1, samples2):
    """Exit with error if any sample IDs overlap between the two folders."""
    overlap = samples1 & samples2
    if overlap:
        ids = ", ".join(sorted(overlap))
        sys.exit(f"ERROR: overlapping sample IDs found in both folders: {ids}")


# ---------------------------------------------------------------------------
# Merge strategies
# ---------------------------------------------------------------------------

def merge_row_append(path1, path2, out_path, skip_header_lines=1):
    """Merge two TSV files by appending rows. Keep header from file 1 only."""
    if path1 is None and path2 is None:
        return
    if path1 is None or path2 is None:
        shutil.copy2(path1 or path2, out_path)
        return

    with open(out_path, "w") as out:
        with open(path1) as f1:
            for line in f1:
                out.write(line)
        with open(path2) as f2:
            for i, line in enumerate(f2):
                if i < skip_header_lines:
                    continue
                out.write(line)


def merge_maf(path1, path2, out_path):
    """Merge two MAF files. Line 1 is #version 2.4, line 2 is header."""
    merge_row_append(path1, path2, out_path, skip_header_lines=2)


def _read_wide_file(path, key_cols):
    """Read a wide-matrix TSV into an indexed structure.

    Returns:
        rows: dict mapping key-column tuple -> dict of {sample_col: value}
        sample_cols: list of sample column names
        key_order: list of key tuples in file order
    """
    rows = {}
    key_order = []
    sample_cols = []

    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        sample_cols = [c for c in reader.fieldnames if c not in key_cols]
        for row in reader:
            key = tuple(row[k] for k in key_cols)
            key_order.append(key)
            rows[key] = {sc: row[sc] for sc in sample_cols}

    return rows, sample_cols, key_order


def merge_wide_matrix(path1, path2, out_path, key_cols, fill_value):
    """Merge two wide-matrix TSVs by outer-joining on key columns."""
    if path1 is None and path2 is None:
        return
    if path1 is None or path2 is None:
        shutil.copy2(path1 or path2, out_path)
        return

    rows1, scols1, order1 = _read_wide_file(path1, key_cols)
    rows2, scols2, order2 = _read_wide_file(path2, key_cols)

    # Stable key ordering: file 1 order, then new keys from file 2
    seen = set(order1)
    all_keys = list(order1)
    for k in order2:
        if k not in seen:
            all_keys.append(k)
            seen.add(k)

    all_sample_cols = sorted(set(scols1) | set(scols2))

    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(key_cols + all_sample_cols)
        for key in all_keys:
            vals1 = rows1.get(key, {})
            vals2 = rows2.get(key, {})
            row = list(key)
            for sc in all_sample_cols:
                v = vals1.get(sc) or vals2.get(sc) or fill_value
                row.append(v)
            writer.writerow(row)


def merge_clinical(path1, path2, out_path):
    """Merge two clinical files with 4-line comment header."""
    if path1 is None and path2 is None:
        return
    if path1 is None or path2 is None:
        shutil.copy2(path1 or path2, out_path)
        return

    with open(path1) as f:
        lines1 = f.readlines()
    with open(path2) as f:
        lines2 = f.readlines()

    # Find where data starts (after # comment lines and column header)
    def split_header_data(lines):
        header = []
        for i, line in enumerate(lines):
            if line.startswith("#"):
                header.append(line)
            else:
                # This line is the column header
                header.append(line)
                return header, lines[i + 1:]
        return header, []

    header1, data1 = split_header_data(lines1)
    _, data2 = split_header_data(lines2)

    with open(out_path, "w") as out:
        for line in header1:
            out.write(line)
        for line in data1:
            out.write(line)
        for line in data2:
            out.write(line)


def copy_if_exists(path1, path2, out_path):
    """Copy a file from whichever folder has it. Prefer path1."""
    if path1 is None and path2 is None:
        return
    shutil.copy2(path1 or path2, out_path)


def _replace_study_id_in_file(path, new_study_id):
    """Rewrite cancer_study_identifier and stable_id prefix in a meta/case list file."""
    with open(path) as f:
        lines = f.readlines()

    old_study_id = None
    for line in lines:
        if line.startswith("cancer_study_identifier:"):
            _, _, value = line.partition(":")
            old_study_id = value.strip()
            break

    with open(path, "w") as f:
        for line in lines:
            if line.startswith("cancer_study_identifier:"):
                f.write(f"cancer_study_identifier: {new_study_id}\n")
            elif line.startswith("stable_id:") and old_study_id:
                _, _, value = line.partition(":")
                value = value.strip()
                if value.startswith(old_study_id):
                    suffix = value[len(old_study_id):]
                    f.write(f"stable_id: {new_study_id}{suffix}\n")
                else:
                    f.write(line)
            else:
                f.write(line)


def _rewrite_meta_study(path, new_study_id):
    """Rewrite meta_study.txt, replacing the old study ID in all fields."""
    with open(path) as f:
        lines = f.readlines()

    old_study_id = None
    for line in lines:
        if line.startswith("cancer_study_identifier:"):
            _, _, value = line.partition(":")
            old_study_id = value.strip()
            break

    if not old_study_id or old_study_id == new_study_id:
        return

    with open(path, "w") as f:
        for line in lines:
            f.write(line.replace(old_study_id, new_study_id))


def rewrite_study_id(out_dir, new_study_id):
    """Rewrite cancer_study_identifier in all meta and case list files."""
    meta_study_path = os.path.join(out_dir, "meta_study.txt")
    if os.path.isfile(meta_study_path):
        _rewrite_meta_study(meta_study_path, new_study_id)

    for filename in META_FILES:
        if filename == "meta_study.txt":
            continue
        path = os.path.join(out_dir, filename)
        if os.path.isfile(path):
            _replace_study_id_in_file(path, new_study_id)

    case_dir = os.path.join(out_dir, "case_lists")
    for filename in CASE_LISTS:
        path = os.path.join(case_dir, filename)
        if os.path.isfile(path):
            _replace_study_id_in_file(path, new_study_id)


def merge_case_list(path1, path2, out_path):
    """Merge two case list files by unioning sample IDs."""
    if path1 is None and path2 is None:
        return
    if path1 is None or path2 is None:
        shutil.copy2(path1 or path2, out_path)
        return

    def parse_case_list(path):
        metadata_lines = []
        sample_ids = set()
        with open(path) as f:
            for line in f:
                if line.startswith("case_list_ids:"):
                    _, _, ids_str = line.partition(":")
                    sample_ids = set(ids_str.strip().split("\t"))
                else:
                    metadata_lines.append(line)
        return metadata_lines, sample_ids

    meta_lines, ids1 = parse_case_list(path1)
    _, ids2 = parse_case_list(path2)

    merged_ids = sorted(ids1 | ids2)

    with open(out_path, "w") as out:
        for line in meta_lines:
            out.write(line)
        out.write("case_list_ids: " + "\t".join(merged_ids) + "\n")


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

def merge_folders(dir1, dir2, out_dir, study_id_override=None):
    """Merge two cBioPortal output folders into out_dir."""
    study_id = validate_inputs(dir1, dir2, study_id_override)
    log(f"Study ID: {study_id}")

    samples1 = collect_sample_ids(dir1)
    samples2 = collect_sample_ids(dir2)
    log(f"Dir 1 samples ({len(samples1)}): {', '.join(sorted(samples1))}")
    log(f"Dir 2 samples ({len(samples2)}): {', '.join(sorted(samples2))}")
    check_sample_overlap(samples1, samples2)
    log("No sample overlap detected")

    os.makedirs(out_dir, exist_ok=False)
    os.makedirs(os.path.join(out_dir, "case_lists"), exist_ok=False)

    # Row-append files
    for filename, skip in ROW_APPEND_FILES:
        p1 = resolve(dir1, filename)
        p2 = resolve(dir2, filename)
        if p1 or p2:
            action = "merging" if (p1 and p2) else "copying"
            log(f"{action} {filename}")
            merge_row_append(p1, p2, os.path.join(out_dir, filename), skip)

    # MAF file
    p1 = resolve(dir1, MAF_FILE)
    p2 = resolve(dir2, MAF_FILE)
    if p1 or p2:
        action = "merging" if (p1 and p2) else "copying"
        log(f"{action} {MAF_FILE}")
        merge_maf(p1, p2, os.path.join(out_dir, MAF_FILE))

    # Wide-matrix files
    for filename, key_cols, fill in WIDE_MATRIX_FILES:
        p1 = resolve(dir1, filename)
        p2 = resolve(dir2, filename)
        if p1 or p2:
            action = "merging" if (p1 and p2) else "copying"
            log(f"{action} {filename}")
            merge_wide_matrix(p1, p2, os.path.join(out_dir, filename), key_cols, fill)

    # Clinical files
    for filename in CLINICAL_FILES:
        p1 = resolve(dir1, filename)
        p2 = resolve(dir2, filename)
        if p1 or p2:
            action = "merging" if (p1 and p2) else "copying"
            log(f"{action} {filename}")
            merge_clinical(p1, p2, os.path.join(out_dir, filename))

    # Meta files
    for filename in META_FILES:
        p1 = resolve(dir1, filename)
        p2 = resolve(dir2, filename)
        if p1 or p2:
            log(f"copying {filename}")
            copy_if_exists(p1, p2, os.path.join(out_dir, filename))

    # Single-value files
    for filename in SINGLE_VALUE_FILES:
        p1 = resolve(dir1, filename)
        p2 = resolve(dir2, filename)
        if p1 or p2:
            log(f"copying {filename}")
            copy_if_exists(p1, p2, os.path.join(out_dir, filename))

    # Case lists
    for filename in CASE_LISTS:
        p1 = resolve(os.path.join(dir1, "case_lists"), filename)
        p2 = resolve(os.path.join(dir2, "case_lists"), filename)
        if p1 or p2:
            action = "merging" if (p1 and p2) else "copying"
            log(f"{action} case_lists/{filename}")
            merge_case_list(p1, p2, os.path.join(out_dir, "case_lists", filename))

    # Rewrite study ID in all meta and case list files if overridden
    if study_id_override:
        log(f"Overriding study ID to '{study_id_override}' in all meta/case list files")
        rewrite_study_id(out_dir, study_id_override)

    total = len(samples1) + len(samples2)
    log(f"Done. {total} samples merged into {out_dir}")


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Merge two cBioPortal output folders into one."
    )
    parser.add_argument(
        "--input_dir_1", required=True, help="First cBioPortal output directory"
    )
    parser.add_argument(
        "--input_dir_2", required=True, help="Second cBioPortal output directory"
    )
    parser.add_argument(
        "--output_dir", required=True, help="Output directory for merged result"
    )
    parser.add_argument(
        "--study_id",
        required=False,
        default=None,
        help="Override the cancer_study_identifier in all output files. "
        "When set, input folders are not required to share the same study ID.",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    dir1 = os.path.abspath(args.input_dir_1)
    dir2 = os.path.abspath(args.input_dir_2)
    out_dir = os.path.abspath(args.output_dir)

    if os.path.exists(out_dir):
        sys.exit(f"ERROR: output directory already exists: {out_dir}")

    merge_folders(dir1, dir2, out_dir, study_id_override=args.study_id)


if __name__ == "__main__":
    main()
