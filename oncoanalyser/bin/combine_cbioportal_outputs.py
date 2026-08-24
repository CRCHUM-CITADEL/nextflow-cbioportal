#!/usr/bin/env python3
"""Merge two cBioPortal output folders into one combined folder.

Inputs and output can be directories or .tar.gz archives containing the
same cBioPortal folder structure.  Archives are auto-extracted/created.

Handled content:
    row-append          data_cna_hg38.seg, data_cna_long.txt, data_sv.txt,
                        util_linking_file.txt
    MAF (2-line header) data_mutations_dna_rna_germline.txt
    wide matrices       data_expression.txt,
                        data_mutational_signatures_{contribution,counts}_{SBS,DBS,ID}.txt
    union-append        data_timeline.txt (column set varies with the event
                        types present in each batch)
    clinical            data_clinical_sample.txt, data_clinical_patient.txt
                        (4 '#' metadata rows + column header; column sets may
                        differ between batches and are unioned)
    headerless          cancer_type.txt (deduplicated on the first field)
    meta files          every meta_*.txt found in either input
    case lists          every case_lists/*.txt found in either input
    subject folders     every other sub-directory (per-subject genomic cache and
                        per-sample clinical splits) is union-copied

Not merged:
    machine_learning/   the processed tables are cohort-normalised (log2,
                        standardised expression), so appending rows would give
                        wrong values.  Regenerate them by re-running the ML
                        subworkflow over the combined cohort.

Anything else found in the inputs is reported as unhandled and left out of the
merged folder; --strict turns that report into a non-zero exit.

Usage:
    combine_cbioportal_outputs.py --input_dir_1 DIR1 --input_dir_2 DIR2 --output_dir OUT
    combine_cbioportal_outputs.py --input_dir_1 DIR1 --input_dir_2 DIR2 --output_dir OUT --study_id NEW_ID

    # Inputs and/or output as .tar.gz:
    combine_cbioportal_outputs.py --input_dir_1 batch1.tar.gz --input_dir_2 batch2.tar.gz --output_dir merged.tar.gz

    # In-place merge (output_dir == one of the input dirs):
    combine_cbioportal_outputs.py --input_dir_1 ACCUMULATED --input_dir_2 NEW_BATCH --output_dir ACCUMULATED

    # Fail if either input contains a file this script does not know how to merge:
    combine_cbioportal_outputs.py --input_dir_1 DIR1 --input_dir_2 DIR2 --output_dir OUT --strict
"""

import argparse
import csv
import glob
import os
import shutil
import sys
import tarfile
import tempfile


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

# Long-format files whose column set depends on which record types a batch
# contained.  gen_timeline.R writes the four common columns first, then the
# union of the remaining columns sorted alphabetically.
UNION_APPEND_FILES = [
    ("data_timeline.txt", ["PATIENT_ID", "START_DATE", "STOP_DATE", "EVENT_TYPE"], ""),
]

CLINICAL_FILES = [
    "data_clinical_sample.txt",
    "data_clinical_patient.txt",
]

# Headerless tab-separated files, deduplicated on the first field.
HEADERLESS_DEDUPE_FILES = [
    "cancer_type.txt",
]

# Sub-directories deliberately left out of the merged folder.
SKIPPED_DIRS = {"machine_learning"}

# Sub-directories with their own merge strategy (i.e. not subject folders).
STRUCTURAL_DIRS = {"case_lists"}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def log(msg):
    """Print a log message to stderr."""
    print(f"[INFO] {msg}", file=sys.stderr)


def warn(msg):
    """Print a warning to stderr."""
    print(f"[WARNING] {msg}", file=sys.stderr)


def resolve(directory, filename):
    """Return the full path if the file exists, else None."""
    path = os.path.join(directory, filename)
    return path if os.path.isfile(path) else None


def is_tarball(path):
    """Return True if path looks like a .tar.gz archive."""
    return path.endswith(".tar.gz") or path.endswith(".tgz")


def extract_tarball(tarball_path):
    """Extract a .tar.gz archive and return (content_dir, tmpdir_root).

    content_dir is the directory containing cBioPortal files (descends into
    a single top-level directory if there is one).  tmpdir_root is the
    extraction root that should be cleaned up.
    """
    tmpdir = tempfile.mkdtemp(prefix="cbio_tar_")
    log(f"Extracting {tarball_path} → {tmpdir}")
    with tarfile.open(tarball_path, "r:gz") as tf:
        try:
            tf.extractall(tmpdir, filter="data")
        except TypeError:
            # Python without the extraction filter API (< 3.8.17 / 3.11.4).
            tf.extractall(tmpdir)
    # If there is exactly one top-level directory, descend into it
    entries = os.listdir(tmpdir)
    if len(entries) == 1:
        candidate = os.path.join(tmpdir, entries[0])
        if os.path.isdir(candidate):
            return candidate, tmpdir
    return tmpdir, tmpdir


def create_tarball(source_dir, tarball_path):
    """Create a .tar.gz archive from a directory."""
    log(f"Creating archive {tarball_path}")
    with tarfile.open(tarball_path, "w:gz") as tf:
        tf.add(source_dir, arcname=os.path.basename(tarball_path).replace(".tar.gz", "").replace(".tgz", ""))


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


def list_meta_files(dir1, dir2):
    """Return the sorted union of meta_*.txt basenames across both folders."""
    names = set()
    for d in (dir1, dir2):
        for path in glob.glob(os.path.join(d, "meta_*.txt")):
            names.add(os.path.basename(path))
    return sorted(names)


def list_case_lists(dir1, dir2):
    """Return the sorted union of case_lists/*.txt basenames across both folders."""
    names = set()
    for d in (dir1, dir2):
        for path in glob.glob(os.path.join(d, "case_lists", "*.txt")):
            names.add(os.path.basename(path))
    return sorted(names)


def list_subject_dirs(directory):
    """Return sub-directories that hold per-subject output."""
    if not os.path.isdir(directory):
        return []
    return sorted(
        name
        for name in os.listdir(directory)
        if os.path.isdir(os.path.join(directory, name))
        and name not in SKIPPED_DIRS
        and name not in STRUCTURAL_DIRS
    )


def _read_tsv_rows(path):
    """Read a single-header TSV into (rows, columns). Missing cells become ''."""
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        columns = list(reader.fieldnames or [])
        rows = []
        for row in reader:
            rows.append({c: ("" if row.get(c) is None else row[c]) for c in columns})
    return rows, columns


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
    """Merge two TSV files by appending rows. Keep header from file 1 only.

    Not every generator terminates its last line: util_linking_file.txt is written
    by a Groovy `join("\n")` in workflows/genomic.nf, so a plain concatenation
    would glue file 1's last row onto file 2's first row.
    """
    if path1 is None and path2 is None:
        return
    if path1 is None or path2 is None:
        shutil.copy2(path1 or path2, out_path)
        return

    with open(out_path, "w") as out:
        last = ""
        with open(path1) as f1:
            for line in f1:
                out.write(line)
                last = line
        if last and not last.endswith("\n"):
            out.write("\n")
        with open(path2) as f2:
            for i, line in enumerate(f2):
                if i < skip_header_lines:
                    continue
                out.write(line)


def merge_maf(path1, path2, out_path):
    """Merge two MAF files. Line 1 is #version 2.4, line 2 is header."""
    merge_row_append(path1, path2, out_path, skip_header_lines=2)


def merge_union_append(path1, path2, out_path, leading_cols, fill_value):
    """Append rows of two long-format TSVs, unioning their column sets.

    The output header is leading_cols (those present in either file) followed by
    the remaining columns sorted alphabetically — the same ordering gen_timeline.R
    produces — so merging two batches gives the same layout as one run over both.
    """
    if path1 is None and path2 is None:
        return
    if path1 is None or path2 is None:
        shutil.copy2(path1 or path2, out_path)
        return

    rows1, cols1 = _read_tsv_rows(path1)
    rows2, cols2 = _read_tsv_rows(path2)

    union = set(cols1) | set(cols2)
    lead = [c for c in leading_cols if c in union]
    rest = sorted(union - set(lead))
    columns = lead + rest

    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        for row in rows1 + rows2:
            writer.writerow([row.get(c, fill_value) for c in columns])


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
            rows[key] = {sc: ("" if row.get(sc) is None else row[sc]) for sc in sample_cols}

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
                # Membership, not truthiness — a real value of "0" or "" must survive.
                if sc in vals1:
                    row.append(vals1[sc])
                elif sc in vals2:
                    row.append(vals2[sc])
                else:
                    row.append(fill_value)
            writer.writerow(row)


def _read_clinical(path):
    """Read a cBioPortal clinical file.

    Returns (col_meta, columns, rows, n_meta_rows) where col_meta maps a column
    name to the tuple of its '#' metadata values (display name, description,
    datatype, priority) and rows is a list of {column: value} dicts.
    """
    with open(path) as f:
        lines = f.read().splitlines()

    meta_rows = []
    idx = 0
    while idx < len(lines) and lines[idx].startswith("#"):
        meta_rows.append(lines[idx][1:].split("\t"))
        idx += 1

    columns = lines[idx].split("\t") if idx < len(lines) else []
    idx += 1

    rows = []
    for line in lines[idx:]:
        if not line.strip():
            continue
        fields = line.split("\t")
        fields += [""] * (len(columns) - len(fields))
        rows.append(dict(zip(columns, fields)))

    col_meta = {
        col: tuple(mr[i] if i < len(mr) else "" for mr in meta_rows)
        for i, col in enumerate(columns)
    }
    return col_meta, columns, rows, len(meta_rows)


def merge_clinical(path1, path2, out_path):
    """Merge two clinical files, unioning their column sets.

    clin_format.R drops columns whose source data is absent, so two batches can
    legitimately carry different attributes. Columns are file-1 order followed by
    the file-2-only columns; the '#' metadata rows are rebuilt per column and
    missing cells are filled with NA.
    """
    if path1 is None and path2 is None:
        return
    if path1 is None or path2 is None:
        shutil.copy2(path1 or path2, out_path)
        return

    meta1, cols1, rows1, n1 = _read_clinical(path1)
    meta2, cols2, rows2, n2 = _read_clinical(path2)

    n_meta = max(n1, n2)
    columns = list(cols1) + [c for c in cols2 if c not in cols1]

    added = [c for c in cols2 if c not in cols1]
    dropped = [c for c in cols1 if c not in cols2]
    if added or dropped:
        name = os.path.basename(out_path)
        log(f"{name}: column sets differ — only in dir1: {dropped or 'none'}; only in dir2: {added or 'none'}")

    col_meta = {}
    for col in columns:
        m1 = meta1.get(col)
        m2 = meta2.get(col)
        if m1 is not None and m2 is not None and m1 != m2:
            warn(f"{os.path.basename(out_path)}: metadata rows differ for column '{col}' — keeping dir1's")
        m = m1 if m1 is not None else (m2 or ())
        col_meta[col] = tuple(m) + ("",) * (n_meta - len(m))

    with open(out_path, "w") as out:
        for i in range(n_meta):
            out.write("#" + "\t".join(col_meta[c][i] for c in columns) + "\n")
        out.write("\t".join(columns) + "\n")
        for row in rows1 + rows2:
            out.write("\t".join(row.get(c, "NA") for c in columns) + "\n")


def merge_headerless_dedupe(path1, path2, out_path):
    """Concatenate two headerless TSVs, dropping rows whose first field repeats."""
    if path1 is None and path2 is None:
        return

    seen = set()
    out_lines = []
    for path in (path1, path2):
        if path is None:
            continue
        with open(path) as f:
            for line in f:
                line = line.rstrip("\n")
                if not line.strip():
                    continue
                key = line.split("\t")[0]
                if key in seen:
                    continue
                seen.add(key)
                out_lines.append(line)

    with open(out_path, "w") as out:
        for line in out_lines:
            out.write(line + "\n")


def copy_if_exists(path1, path2, out_path):
    """Copy a file from whichever folder has it. Prefer path1."""
    if path1 is None and path2 is None:
        return
    shutil.copy2(path1 or path2, out_path)


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


def merge_subject_dirs(dir1, dir2, out_dir):
    """Union-copy the per-subject folders of both inputs into out_dir."""
    subs1 = set(list_subject_dirs(dir1))
    subs2 = set(list_subject_dirs(dir2))

    for name in sorted(subs1 | subs2):
        dest = os.path.join(out_dir, name)
        if name in subs1 and name not in subs2:
            shutil.copytree(os.path.join(dir1, name), dest)
        elif name in subs2 and name not in subs1:
            shutil.copytree(os.path.join(dir2, name), dest)
        else:
            warn(f"subject '{name}' is present in both folders — keeping dir1's copy of shared files")
            shutil.copytree(os.path.join(dir1, name), dest)
            conflicts = []
            src_dir = os.path.join(dir2, name)
            for entry in sorted(os.listdir(src_dir)):
                src = os.path.join(src_dir, entry)
                tgt = os.path.join(dest, entry)
                if os.path.exists(tgt):
                    conflicts.append(entry)
                    continue
                if os.path.isdir(src):
                    shutil.copytree(src, tgt)
                else:
                    shutil.copy2(src, tgt)
            if conflicts:
                warn(f"subject '{name}': kept dir1's version of {', '.join(conflicts)}")

    return sorted(subs1 | subs2)


# ---------------------------------------------------------------------------
# Study ID rewriting
# ---------------------------------------------------------------------------

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

    if old_study_id is None:
        # e.g. meta_cancer_type.txt, which carries no study identifier at all.
        log(f"{os.path.basename(path)}: no cancer_study_identifier — left unchanged")
        return

    with open(path, "w") as f:
        for line in lines:
            if line.startswith("cancer_study_identifier:"):
                f.write(f"cancer_study_identifier: {new_study_id}\n")
            elif line.startswith("stable_id:"):
                _, _, value = line.partition(":")
                value = value.strip()
                if value.startswith(old_study_id):
                    suffix = value[len(old_study_id):]
                    f.write(f"stable_id: {new_study_id}{suffix}\n")
                else:
                    # Genomic/clinical meta stable_ids are fixed, not study-prefixed.
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


def rewrite_study_id(out_dir, new_study_id, meta_names, case_list_names):
    """Rewrite cancer_study_identifier in all meta and case list files."""
    meta_study_path = os.path.join(out_dir, "meta_study.txt")
    if os.path.isfile(meta_study_path):
        _rewrite_meta_study(meta_study_path, new_study_id)

    for filename in meta_names:
        if filename == "meta_study.txt":
            continue
        path = os.path.join(out_dir, filename)
        if os.path.isfile(path):
            _replace_study_id_in_file(path, new_study_id)

    case_dir = os.path.join(out_dir, "case_lists")
    for filename in case_list_names:
        path = os.path.join(case_dir, filename)
        if os.path.isfile(path):
            _replace_study_id_in_file(path, new_study_id)


# ---------------------------------------------------------------------------
# Unhandled-content report
# ---------------------------------------------------------------------------

def find_unhandled(dir1, dir2, handled_files, handled_case_lists, subject_dirs):
    """Return input entries that no merge strategy claimed."""
    unhandled = set()
    for d in (dir1, dir2):
        for entry in sorted(os.listdir(d)):
            path = os.path.join(d, entry)
            if os.path.isdir(path):
                if entry in SKIPPED_DIRS or entry in STRUCTURAL_DIRS or entry in subject_dirs:
                    continue
                unhandled.add(entry + "/")
            elif entry not in handled_files:
                unhandled.add(entry)

        case_dir = os.path.join(d, "case_lists")
        if os.path.isdir(case_dir):
            for entry in sorted(os.listdir(case_dir)):
                if entry not in handled_case_lists:
                    unhandled.add(f"case_lists/{entry}")

    return sorted(unhandled)


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

def merge_folders(dir1, dir2, out_dir, study_id_override=None, strict=False):
    """Merge two cBioPortal output folders into out_dir.

    When out_dir is the same path as dir1 or dir2, an in-place merge is
    performed: data is merged into a temporary directory, then the original
    is replaced with the result.
    """
    in_place = os.path.normpath(out_dir) in (
        os.path.normpath(dir1),
        os.path.normpath(dir2),
    )

    actual_out = None
    if in_place:
        actual_out = out_dir
        out_dir = tempfile.mkdtemp(
            prefix="cbio_merge_", dir=os.path.dirname(actual_out)
        )
        log(f"In-place merge detected — using temp directory: {out_dir}")

    study_id = validate_inputs(dir1, dir2, study_id_override)
    log(f"Study ID: {study_id}")

    samples1 = collect_sample_ids(dir1)
    samples2 = collect_sample_ids(dir2)
    log(f"Dir 1 samples ({len(samples1)}): {', '.join(sorted(samples1))}")
    log(f"Dir 2 samples ({len(samples2)}): {', '.join(sorted(samples2))}")
    check_sample_overlap(samples1, samples2)
    log("No sample overlap detected")

    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(os.path.join(out_dir, "case_lists"), exist_ok=True)

    handled_files = set()

    # Row-append files
    for filename, skip in ROW_APPEND_FILES:
        handled_files.add(filename)
        p1 = resolve(dir1, filename)
        p2 = resolve(dir2, filename)
        if p1 or p2:
            action = "merging" if (p1 and p2) else "copying"
            log(f"{action} {filename}")
            merge_row_append(p1, p2, os.path.join(out_dir, filename), skip)

    # MAF file
    handled_files.add(MAF_FILE)
    p1 = resolve(dir1, MAF_FILE)
    p2 = resolve(dir2, MAF_FILE)
    if p1 or p2:
        action = "merging" if (p1 and p2) else "copying"
        log(f"{action} {MAF_FILE}")
        merge_maf(p1, p2, os.path.join(out_dir, MAF_FILE))

    # Wide-matrix files
    for filename, key_cols, fill in WIDE_MATRIX_FILES:
        handled_files.add(filename)
        p1 = resolve(dir1, filename)
        p2 = resolve(dir2, filename)
        if p1 or p2:
            action = "merging" if (p1 and p2) else "copying"
            log(f"{action} {filename}")
            merge_wide_matrix(p1, p2, os.path.join(out_dir, filename), key_cols, fill)

    # Union-append files (timeline)
    for filename, leading_cols, fill in UNION_APPEND_FILES:
        handled_files.add(filename)
        p1 = resolve(dir1, filename)
        p2 = resolve(dir2, filename)
        if p1 or p2:
            action = "merging" if (p1 and p2) else "copying"
            log(f"{action} {filename}")
            merge_union_append(p1, p2, os.path.join(out_dir, filename), leading_cols, fill)

    # Clinical files
    for filename in CLINICAL_FILES:
        handled_files.add(filename)
        p1 = resolve(dir1, filename)
        p2 = resolve(dir2, filename)
        if p1 or p2:
            action = "merging" if (p1 and p2) else "copying"
            log(f"{action} {filename}")
            merge_clinical(p1, p2, os.path.join(out_dir, filename))

    # Headerless files
    for filename in HEADERLESS_DEDUPE_FILES:
        handled_files.add(filename)
        p1 = resolve(dir1, filename)
        p2 = resolve(dir2, filename)
        if p1 or p2:
            action = "merging" if (p1 and p2) else "copying"
            log(f"{action} {filename}")
            merge_headerless_dedupe(p1, p2, os.path.join(out_dir, filename))

    # Meta files — discovered, so new ones are never silently dropped
    meta_names = list_meta_files(dir1, dir2)
    for filename in meta_names:
        handled_files.add(filename)
        p1 = resolve(dir1, filename)
        p2 = resolve(dir2, filename)
        log(f"copying {filename}")
        copy_if_exists(p1, p2, os.path.join(out_dir, filename))

    # Case lists — also discovered
    case_list_names = list_case_lists(dir1, dir2)
    for filename in case_list_names:
        p1 = resolve(os.path.join(dir1, "case_lists"), filename)
        p2 = resolve(os.path.join(dir2, "case_lists"), filename)
        action = "merging" if (p1 and p2) else "copying"
        log(f"{action} case_lists/{filename}")
        merge_case_list(p1, p2, os.path.join(out_dir, "case_lists", filename))

    # Per-subject folders
    subject_dirs = merge_subject_dirs(dir1, dir2, out_dir)
    if subject_dirs:
        log(f"copied {len(subject_dirs)} subject folder(s): {', '.join(subject_dirs)}")

    # Rewrite study ID in all meta and case list files if overridden
    if study_id_override:
        log(f"Overriding study ID to '{study_id_override}' in all meta/case list files")
        rewrite_study_id(out_dir, study_id_override, meta_names, case_list_names)

    # Report anything the dispatch tables did not claim
    for skipped in sorted(SKIPPED_DIRS):
        if os.path.isdir(os.path.join(dir1, skipped)) or os.path.isdir(os.path.join(dir2, skipped)):
            warn(
                f"{skipped}/ was NOT merged — its tables are cohort-normalised. "
                "Regenerate them over the combined cohort."
            )

    unhandled = find_unhandled(dir1, dir2, handled_files, case_list_names, subject_dirs)
    if unhandled:
        warn("unhandled input(s), left out of the merged folder: " + ", ".join(unhandled))
        warn("add them to a dispatch table in combine_cbioportal_outputs.py")

    # Swap temp directory into place for in-place merges
    if in_place:
        backup = f"{actual_out}.bak_{os.getpid()}"
        shutil.move(actual_out, backup)
        try:
            shutil.move(out_dir, actual_out)
        except Exception:
            shutil.move(backup, actual_out)
            raise
        shutil.rmtree(backup, ignore_errors=True)
        out_dir = actual_out
        log(f"Replaced {actual_out} with merged result")

    total = len(samples1) + len(samples2)
    log(f"Done. {total} samples merged into {out_dir}")

    if unhandled and strict:
        sys.exit("ERROR: --strict and unhandled input(s) present: " + ", ".join(unhandled))


def check_archive_completeness(dir_path, tarball_path):
    """Warn when an extracted archive looks like a genomic-only package."""
    has_study = os.path.isfile(os.path.join(dir_path, "meta_study.txt"))
    has_clinical = os.path.isfile(os.path.join(dir_path, "data_clinical_sample.txt"))
    if has_study and not has_clinical:
        warn(
            f"{os.path.basename(tarball_path)} has no data_clinical_sample.txt — it looks like a "
            "genomic-only archive. Use the study output DIRECTORY instead to keep clinical data."
        )


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Merge two cBioPortal output folders into one."
    )
    parser.add_argument(
        "--input_dir_1", required=True, help="First cBioPortal output directory or .tar.gz archive"
    )
    parser.add_argument(
        "--input_dir_2", required=True, help="Second cBioPortal output directory or .tar.gz archive"
    )
    parser.add_argument(
        "--output_dir", required=True, help="Output directory or .tar.gz archive for merged result"
    )
    parser.add_argument(
        "--study_id",
        required=False,
        default=None,
        help="Override the cancer_study_identifier in all output files. "
        "When set, input folders are not required to share the same study ID.",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Exit non-zero if either input contains a file no merge strategy handles.",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    input1 = os.path.abspath(args.input_dir_1)
    input2 = os.path.abspath(args.input_dir_2)
    output = os.path.abspath(args.output_dir)

    # Extract tar.gz inputs to temp directories
    tmp_dirs = []
    if is_tarball(input1):
        if not os.path.isfile(input1):
            sys.exit(f"ERROR: archive not found: {input1}")
        dir1, tmproot1 = extract_tarball(input1)
        tmp_dirs.append(tmproot1)
        check_archive_completeness(dir1, input1)
    else:
        dir1 = input1

    if is_tarball(input2):
        if not os.path.isfile(input2):
            sys.exit(f"ERROR: archive not found: {input2}")
        dir2, tmproot2 = extract_tarball(input2)
        tmp_dirs.append(tmproot2)
        check_archive_completeness(dir2, input2)
    else:
        dir2 = input2

    tar_output = is_tarball(output)
    if tar_output:
        out_dir = tempfile.mkdtemp(prefix="cbio_out_")
        tmp_dirs.append(out_dir)
        if os.path.exists(output):
            sys.exit(f"ERROR: output archive already exists: {output}")
    else:
        out_dir = output

    try:
        in_place = os.path.normpath(out_dir) in (
            os.path.normpath(dir1),
            os.path.normpath(dir2),
        )

        if os.path.exists(out_dir) and not in_place and not tar_output:
            sys.exit(f"ERROR: output directory already exists: {out_dir}")

        merge_folders(
            dir1, dir2, out_dir,
            study_id_override=args.study_id,
            strict=args.strict,
        )

        if tar_output:
            create_tarball(out_dir, output)
    finally:
        for d in tmp_dirs:
            if os.path.isdir(d):
                shutil.rmtree(d, ignore_errors=True)


if __name__ == "__main__":
    main()
