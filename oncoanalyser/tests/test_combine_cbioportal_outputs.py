"""Tests for bin/combine_cbioportal_outputs.py."""

import subprocess
import sys
import tarfile
from pathlib import Path

SCRIPT = Path(__file__).parent.parent / "bin" / "combine_cbioportal_outputs.py"


def run(args, **kwargs):
    """Run the script and return CompletedProcess."""
    return subprocess.run(
        [sys.executable, str(SCRIPT)] + args,
        capture_output=True,
        text=True,
        **kwargs,
    )


# ── Fixture builders ──────────────────────────────────────────────────────────

def write_tsv(path, header, rows):
    """Write a single-header TSV."""
    lines = ["\t".join(header)] + ["\t".join(r) for r in rows]
    path.write_text("\n".join(lines) + "\n")


def write_clinical(path, columns, rows):
    """Write a cBioPortal clinical file: 4 '#' metadata rows + header + data."""
    meta_rows = [
        [f"{c} name" for c in columns],
        [f"{c} description" for c in columns],
        ["STRING" for _ in columns],
        ["1" for _ in columns],
    ]
    lines = ["#" + "\t".join(mr) for mr in meta_rows]
    lines.append("\t".join(columns))
    lines += ["\t".join(r) for r in rows]
    path.write_text("\n".join(lines) + "\n")


def read_clinical(path):
    """Return (meta_rows, columns, rows-as-dicts) from a clinical file."""
    lines = path.read_text().splitlines()
    meta_rows = [ln[1:].split("\t") for ln in lines if ln.startswith("#")]
    body = [ln for ln in lines if not ln.startswith("#")]
    columns = body[0].split("\t")
    rows = [dict(zip(columns, ln.split("\t"))) for ln in body[1:] if ln.strip()]
    return meta_rows, columns, rows


def read_tsv(path):
    """Return (columns, rows-as-dicts) from a single-header TSV."""
    lines = [ln for ln in path.read_text().splitlines() if ln.strip()]
    columns = lines[0].split("\t")
    rows = [dict(zip(columns, ln.split("\t"))) for ln in lines[1:]]
    return columns, rows


def write_case_list(path, study_id, label, sample_ids):
    """Write a cases_*.txt case list."""
    path.write_text(
        f"cancer_study_identifier: {study_id}\n"
        f"stable_id: {study_id}_{label}\n"
        "case_list_name: add_text\n"
        "case_list_description: ADD TEXT\n"
        "case_list_ids: " + "\t".join(sample_ids) + "\n"
    )


def make_study(
    root,
    name,
    pairs,
    study_id="test_study",
    timeline_columns=None,
    clinical_sample_columns=None,
    expression=None,
    cancer_type=True,
    timeline=True,
    subject_dirs=True,
    machine_learning=False,
    extra_case_lists=None,
    extra_files=None,
):
    """Build a minimal but realistic cBioPortal study folder.

    pairs is a list of (subject, sample) tuples.
    """
    d = root / name
    (d / "case_lists").mkdir(parents=True)
    subjects = [s for s, _ in pairs]
    samples = [s for _, s in pairs]

    (d / "meta_study.txt").write_text(
        f"type_of_cancer: mixed\ncancer_study_identifier: {study_id}\n"
        f"name: {name}\ndescription: test\nadd_global_case_list: true\nreference_genome: hg38\n"
    )
    for label, extra in [
        ("cna_long", "stable_id: cna\n"),
        ("sv", "stable_id: structural_variants\n"),
        ("clinical_sample", ""),
        ("clinical_patient", ""),
    ]:
        (d / f"meta_{label}.txt").write_text(
            f"cancer_study_identifier: {study_id}\ngenetic_alteration_type: X\n{extra}"
        )

    # Row-append data
    write_tsv(
        d / "data_cna_long.txt",
        ["Hugo_Symbol", "Entrez_Gene_Id", "Sample_Id", "Value"],
        [["TP53", "7157", s, "-1"] for s in samples],
    )
    write_tsv(
        d / "data_sv.txt",
        ["Sample_Id", "SV_Status", "Site1_Hugo_Symbol"],
        [[s, "SOMATIC", "BRAF"] for s in samples],
    )
    write_tsv(
        d / "util_linking_file.txt",
        ["subject_id", "sample_id"],
        [[sub, sam] for sub, sam in pairs],
    )

    # MAF (2-line header)
    (d / "data_mutations_dna_rna_germline.txt").write_text(
        "#version 2.4\nHugo_Symbol\tTumor_Sample_Barcode\n"
        + "".join(f"TP53\t{s}\n" for s in samples)
    )

    # Wide matrix
    if expression is None:
        expression = {"TP53": {s: "1.5" for s in samples}}
    genes = sorted(expression)
    write_tsv(
        d / "data_expression.txt",
        ["Hugo_Symbol", "Entrez_Gene_Id"] + samples,
        [[g, "7157"] + [expression[g].get(s, "NA") for s in samples] for g in genes],
    )

    # Clinical
    sample_cols = clinical_sample_columns or ["PATIENT_ID", "SAMPLE_ID", "SAMPLE_TYPE"]
    write_clinical(
        d / "data_clinical_sample.txt",
        sample_cols,
        [[sub, sam] + ["v"] * (len(sample_cols) - 2) for sub, sam in pairs],
    )
    write_clinical(
        d / "data_clinical_patient.txt",
        ["PATIENT_ID", "SEX"],
        [[sub, "Male"] for sub in subjects],
    )

    # Timeline
    if timeline:
        cols = timeline_columns or [
            "PATIENT_ID", "START_DATE", "STOP_DATE", "EVENT_TYPE", "SUBTYPE",
        ]
        rows = [[sub, "0", "10", "SURGERY"] + ["x"] * (len(cols) - 4) for sub in subjects]
        write_tsv(d / "data_timeline.txt", cols, rows)
        (d / "meta_timeline.txt").write_text(
            f"cancer_study_identifier: {study_id}\ngenetic_alteration_type: CLINICAL\n"
            "datatype: TIMELINE\ndata_filename: data_timeline.txt\n"
        )

    # Cancer type
    if cancer_type:
        (d / "cancer_type.txt").write_text("mixed\tMixed\tBlack\ttissue\n")
        (d / "meta_cancer_type.txt").write_text(
            "genetic_alteration_type: CANCER_TYPE\ndatatype: CANCER_TYPE\n"
            "data_filename: cancer_type.txt\n"
        )

    # Case lists
    write_case_list(d / "case_lists" / "cases_sequenced.txt", study_id, "sequenced", samples)
    write_case_list(d / "case_lists" / "cases_cnv.txt", study_id, "cnv", samples)
    for label in extra_case_lists or []:
        write_case_list(d / "case_lists" / f"cases_{label}.txt", study_id, label, samples)

    # Per-subject folders
    if subject_dirs:
        for sub, sam in pairs:
            sd = d / sub
            sd.mkdir(exist_ok=True)
            (sd / f"{sam}.tpm.tsv").write_text(f"gene\t{sam}\nTP53\t1.0\n")
            (sd / f"{sam}.data_clinical_sample.txt").write_text("stub\n")

    if machine_learning:
        ml = d / "machine_learning" / "processed"
        ml.mkdir(parents=True)
        (ml / "expression_processed_standardized.tsv").write_text("gene\tv\nTP53\t0.1\n")

    for fname, content in (extra_files or {}).items():
        (d / fname).write_text(content)

    return d


def merge(tmp_path, dir1, dir2, out_name="merged", extra_args=None):
    """Run a merge and assert it succeeded."""
    out = tmp_path / out_name
    result = run(
        ["--input_dir_1", str(dir1), "--input_dir_2", str(dir2), "--output_dir", str(out)]
        + (extra_args or [])
    )
    assert result.returncode == 0, result.stderr
    return out, result


# ── Happy paths ───────────────────────────────────────────────────────────────

def test_row_append_and_maf(tmp_path):
    d1 = make_study(tmp_path, "b1", [("P1", "S1")])
    d2 = make_study(tmp_path, "b2", [("P2", "S2")])
    out, _ = merge(tmp_path, d1, d2)

    _, sv = read_tsv(out / "data_sv.txt")
    assert [r["Sample_Id"] for r in sv] == ["S1", "S2"]

    maf = (out / "data_mutations_dna_rna_germline.txt").read_text().splitlines()
    assert maf[0] == "#version 2.4"
    assert maf.count("#version 2.4") == 1
    assert maf[1].startswith("Hugo_Symbol")
    assert len([ln for ln in maf if ln.startswith("TP53")]) == 2

    _, linking = read_tsv(out / "util_linking_file.txt")
    assert {r["sample_id"] for r in linking} == {"S1", "S2"}


def test_row_append_when_first_file_lacks_trailing_newline(tmp_path):
    # util_linking_file.txt is written with a Groovy join("\n") and has no final
    # newline, so a naive concatenation glues two rows together.
    d1 = make_study(tmp_path, "b1", [("P1", "S1")])
    d2 = make_study(tmp_path, "b2", [("P2", "S2")])
    for d, sub, sam in [(d1, "P1", "S1"), (d2, "P2", "S2")]:
        (d / "util_linking_file.txt").write_text(f"subject_id\tsample_id\n{sub}\t{sam}")

    out, _ = merge(tmp_path, d1, d2)
    lines = (out / "util_linking_file.txt").read_text().splitlines()
    assert lines == ["subject_id\tsample_id", "P1\tS1", "P2\tS2"]


def test_timeline_column_union(tmp_path):
    d1 = make_study(
        tmp_path, "b1", [("P1", "S1")],
        timeline_columns=["PATIENT_ID", "START_DATE", "STOP_DATE", "EVENT_TYPE", "SUBTYPE", "TEST"],
    )
    d2 = make_study(
        tmp_path, "b2", [("P2", "S2")],
        timeline_columns=["PATIENT_ID", "START_DATE", "STOP_DATE", "EVENT_TYPE", "SPECIMEN_TYPE", "STATUS"],
    )
    out, _ = merge(tmp_path, d1, d2)

    columns, rows = read_tsv(out / "data_timeline.txt")
    # Four common columns first, then the union sorted alphabetically — the
    # ordering gen_timeline.R produces for a single run over both batches.
    assert columns == [
        "PATIENT_ID", "START_DATE", "STOP_DATE", "EVENT_TYPE",
        "SPECIMEN_TYPE", "STATUS", "SUBTYPE", "TEST",
    ]
    assert len(rows) == 2
    p1 = next(r for r in rows if r["PATIENT_ID"] == "P1")
    p2 = next(r for r in rows if r["PATIENT_ID"] == "P2")
    # Values stay under their own columns; the other batch's columns are blank.
    assert p1["SUBTYPE"] == "x" and p1["TEST"] == "x"
    assert p1["STATUS"] == "" and p1["SPECIMEN_TYPE"] == ""
    assert p2["STATUS"] == "x" and p2["SPECIMEN_TYPE"] == "x"
    assert p2["SUBTYPE"] == "" and p2["TEST"] == ""


def test_timeline_present_in_only_one_dir(tmp_path):
    d1 = make_study(tmp_path, "b1", [("P1", "S1")])
    d2 = make_study(tmp_path, "b2", [("P2", "S2")], timeline=False)
    out, _ = merge(tmp_path, d1, d2)

    _, rows = read_tsv(out / "data_timeline.txt")
    assert [r["PATIENT_ID"] for r in rows] == ["P1"]


def test_timeline_absent_from_both(tmp_path):
    d1 = make_study(tmp_path, "b1", [("P1", "S1")], timeline=False)
    d2 = make_study(tmp_path, "b2", [("P2", "S2")], timeline=False)
    out, _ = merge(tmp_path, d1, d2)
    assert not (out / "data_timeline.txt").exists()


def test_clinical_column_union(tmp_path):
    d1 = make_study(
        tmp_path, "b1", [("P1", "S1")],
        clinical_sample_columns=["PATIENT_ID", "SAMPLE_ID", "SAMPLE_TYPE", "PSA_LEVEL"],
    )
    d2 = make_study(
        tmp_path, "b2", [("P2", "S2")],
        clinical_sample_columns=["PATIENT_ID", "SAMPLE_ID", "SAMPLE_TYPE", "CEA"],
    )
    out, _ = merge(tmp_path, d1, d2)

    meta_rows, columns, rows = read_clinical(out / "data_clinical_sample.txt")
    assert columns == ["PATIENT_ID", "SAMPLE_ID", "SAMPLE_TYPE", "PSA_LEVEL", "CEA"]
    assert len(meta_rows) == 4
    # Every '#' row stays aligned with the merged header.
    for mr in meta_rows:
        assert len(mr) == len(columns)
    assert meta_rows[0][-1] == "CEA name"
    assert meta_rows[0][-2] == "PSA_LEVEL name"

    s1 = next(r for r in rows if r["SAMPLE_ID"] == "S1")
    s2 = next(r for r in rows if r["SAMPLE_ID"] == "S2")
    assert s1["PSA_LEVEL"] == "v" and s1["CEA"] == "NA"
    assert s2["CEA"] == "v" and s2["PSA_LEVEL"] == "NA"


def test_clinical_patient_rows_are_not_deduplicated(tmp_path):
    # Same patient, different samples in each batch: both rows are kept.
    d1 = make_study(tmp_path, "b1", [("P1", "S1")])
    d2 = make_study(tmp_path, "b2", [("P1", "S2")])
    out, _ = merge(tmp_path, d1, d2)

    _, _, rows = read_clinical(out / "data_clinical_patient.txt")
    assert [r["PATIENT_ID"] for r in rows] == ["P1", "P1"]


def test_wide_matrix_preserves_zero(tmp_path):
    # A real expression value of 0 must not be replaced by the fill value.
    d1 = make_study(
        tmp_path, "b1", [("P1", "S1")],
        expression={"TP53": {"S1": "0"}, "BRAF": {"S1": "2.0"}},
    )
    d2 = make_study(
        tmp_path, "b2", [("P2", "S2")],
        expression={"TP53": {"S2": "3.0"}},
    )
    out, _ = merge(tmp_path, d1, d2)

    columns, rows = read_tsv(out / "data_expression.txt")
    assert columns == ["Hugo_Symbol", "Entrez_Gene_Id", "S1", "S2"]
    tp53 = next(r for r in rows if r["Hugo_Symbol"] == "TP53")
    braf = next(r for r in rows if r["Hugo_Symbol"] == "BRAF")
    assert tp53["S1"] == "0"
    assert tp53["S2"] == "3.0"
    # Gene absent from batch 2 gets the fill value for that batch's sample.
    assert braf["S1"] == "2.0"
    assert braf["S2"] == "NA"


def test_cancer_type_deduplicated(tmp_path):
    d1 = make_study(tmp_path, "b1", [("P1", "S1")])
    d2 = make_study(tmp_path, "b2", [("P2", "S2")])
    out, _ = merge(tmp_path, d1, d2)

    lines = (out / "cancer_type.txt").read_text().splitlines()
    assert lines == ["mixed\tMixed\tBlack\ttissue"]


def test_meta_files_are_discovered_not_hardcoded(tmp_path):
    d1 = make_study(tmp_path, "b1", [("P1", "S1")],
                    extra_files={"meta_future_assay.txt": "cancer_study_identifier: test_study\n"})
    d2 = make_study(tmp_path, "b2", [("P2", "S2")])
    out, _ = merge(tmp_path, d1, d2)

    for name in ("meta_timeline.txt", "meta_cancer_type.txt", "meta_future_assay.txt"):
        assert (out / name).exists(), f"{name} missing from merged output"


def test_case_lists_are_discovered_not_hardcoded(tmp_path):
    d1 = make_study(tmp_path, "b1", [("P1", "S1")], extra_case_lists=["rna_seq"])
    d2 = make_study(tmp_path, "b2", [("P2", "S2")], extra_case_lists=["rna_seq"])
    out, _ = merge(tmp_path, d1, d2)

    text = (out / "case_lists" / "cases_rna_seq.txt").read_text()
    assert "case_list_ids: S1\tS2" in text
    sequenced = (out / "case_lists" / "cases_sequenced.txt").read_text()
    assert "case_list_ids: S1\tS2" in sequenced


def test_subject_dirs_unioned(tmp_path):
    d1 = make_study(tmp_path, "b1", [("P1", "S1")])
    d2 = make_study(tmp_path, "b2", [("P2", "S2")])
    out, _ = merge(tmp_path, d1, d2)

    assert (out / "P1" / "S1.tpm.tsv").exists()
    assert (out / "P1" / "S1.data_clinical_sample.txt").exists()
    assert (out / "P2" / "S2.tpm.tsv").exists()


def test_subject_dir_present_in_both_warns_and_keeps_union(tmp_path):
    d1 = make_study(tmp_path, "b1", [("P1", "S1")])
    d2 = make_study(tmp_path, "b2", [("P1", "S2")])
    out, result = merge(tmp_path, d1, d2)

    assert "subject 'P1' is present in both folders" in result.stderr
    assert (out / "P1" / "S1.tpm.tsv").exists()
    assert (out / "P1" / "S2.tpm.tsv").exists()


def test_machine_learning_is_skipped_with_warning(tmp_path):
    d1 = make_study(tmp_path, "b1", [("P1", "S1")], machine_learning=True)
    d2 = make_study(tmp_path, "b2", [("P2", "S2")], machine_learning=True)
    out, result = merge(tmp_path, d1, d2)

    assert not (out / "machine_learning").exists()
    assert "machine_learning/ was NOT merged" in result.stderr
    # A deliberate skip is not an unhandled input.
    assert "unhandled input" not in result.stderr


def test_study_id_override(tmp_path):
    d1 = make_study(tmp_path, "b1", [("P1", "S1")], study_id="old_study")
    d2 = make_study(tmp_path, "b2", [("P2", "S2")], study_id="other_study")
    out, _ = merge(tmp_path, d1, d2, extra_args=["--study_id", "new_study"])

    assert "cancer_study_identifier: new_study" in (out / "meta_study.txt").read_text()
    assert "cancer_study_identifier: new_study" in (out / "meta_timeline.txt").read_text()

    case_list = (out / "case_lists" / "cases_sequenced.txt").read_text()
    assert "cancer_study_identifier: new_study" in case_list
    assert "stable_id: new_study_sequenced" in case_list

    # Genomic meta stable_ids are fixed strings, not study-prefixed.
    assert "stable_id: cna" in (out / "meta_cna_long.txt").read_text()

    # meta_cancer_type.txt carries no study identifier and must survive untouched.
    assert (out / "meta_cancer_type.txt").read_text() == (
        "genetic_alteration_type: CANCER_TYPE\ndatatype: CANCER_TYPE\n"
        "data_filename: cancer_type.txt\n"
    )


def test_in_place_merge(tmp_path):
    d1 = make_study(tmp_path, "b1", [("P1", "S1")])
    d2 = make_study(tmp_path, "b2", [("P2", "S2")])
    result = run([
        "--input_dir_1", str(d1), "--input_dir_2", str(d2), "--output_dir", str(d1),
    ])
    assert result.returncode == 0, result.stderr

    _, sv = read_tsv(d1 / "data_sv.txt")
    assert [r["Sample_Id"] for r in sv] == ["S1", "S2"]
    assert (d1 / "P2" / "S2.tpm.tsv").exists()
    # No backup directory left behind.
    assert not list(tmp_path.glob("b1.bak_*"))


def test_tarball_roundtrip(tmp_path):
    d1 = make_study(tmp_path, "b1", [("P1", "S1")])
    d2 = make_study(tmp_path, "b2", [("P2", "S2")])
    for d in (d1, d2):
        with tarfile.open(str(tmp_path / f"{d.name}.tar.gz"), "w:gz") as tf:
            tf.add(str(d), arcname=d.name)

    out_tar = tmp_path / "merged.tar.gz"
    result = run([
        "--input_dir_1", str(tmp_path / "b1.tar.gz"),
        "--input_dir_2", str(tmp_path / "b2.tar.gz"),
        "--output_dir", str(out_tar),
    ])
    assert result.returncode == 0, result.stderr

    with tarfile.open(str(out_tar), "r:gz") as tf:
        names = {n.split("/", 1)[1] for n in tf.getnames() if "/" in n}
        tf.extractall(str(tmp_path / "unpacked"))
    assert "data_timeline.txt" in names
    assert "cancer_type.txt" in names
    assert "case_lists/cases_sequenced.txt" in names
    assert "P1/S1.tpm.tsv" in names

    _, sv = read_tsv(tmp_path / "unpacked" / "merged" / "data_sv.txt")
    assert [r["Sample_Id"] for r in sv] == ["S1", "S2"]


def test_genomic_only_archive_warns(tmp_path):
    d1 = make_study(tmp_path, "b1", [("P1", "S1")])
    (d1 / "data_clinical_sample.txt").unlink()
    with tarfile.open(str(tmp_path / "b1.tar.gz"), "w:gz") as tf:
        tf.add(str(d1), arcname="b1")
    d2 = make_study(tmp_path, "b2", [("P2", "S2")])

    result = run([
        "--input_dir_1", str(tmp_path / "b1.tar.gz"),
        "--input_dir_2", str(d2),
        "--output_dir", str(tmp_path / "merged"),
    ])
    assert result.returncode == 0, result.stderr
    assert "looks like a genomic-only archive" in result.stderr


# ── Error paths ───────────────────────────────────────────────────────────────

def test_unhandled_file_warns(tmp_path):
    d1 = make_study(tmp_path, "b1", [("P1", "S1")],
                    extra_files={"data_future_assay.txt": "a\tb\n1\t2\n"})
    d2 = make_study(tmp_path, "b2", [("P2", "S2")])
    out, result = merge(tmp_path, d1, d2)

    assert "unhandled input(s)" in result.stderr
    assert "data_future_assay.txt" in result.stderr
    assert not (out / "data_future_assay.txt").exists()


def test_unhandled_file_strict_exits_non_zero(tmp_path):
    d1 = make_study(tmp_path, "b1", [("P1", "S1")],
                    extra_files={"data_future_assay.txt": "a\tb\n1\t2\n"})
    d2 = make_study(tmp_path, "b2", [("P2", "S2")])
    result = run([
        "--input_dir_1", str(d1), "--input_dir_2", str(d2),
        "--output_dir", str(tmp_path / "merged"), "--strict",
    ])
    assert result.returncode != 0
    assert "ERROR" in result.stderr
    assert "data_future_assay.txt" in result.stderr


def test_clean_study_passes_strict(tmp_path):
    d1 = make_study(tmp_path, "b1", [("P1", "S1")], machine_learning=True)
    d2 = make_study(tmp_path, "b2", [("P2", "S2")], machine_learning=True)
    result = run([
        "--input_dir_1", str(d1), "--input_dir_2", str(d2),
        "--output_dir", str(tmp_path / "merged"), "--strict",
    ])
    assert result.returncode == 0, result.stderr


def test_overlapping_sample_ids_error(tmp_path):
    d1 = make_study(tmp_path, "b1", [("P1", "S1")])
    d2 = make_study(tmp_path, "b2", [("P2", "S1")])
    result = run([
        "--input_dir_1", str(d1), "--input_dir_2", str(d2),
        "--output_dir", str(tmp_path / "merged"),
    ])
    assert result.returncode != 0
    assert "ERROR: overlapping sample IDs" in result.stderr
    assert "S1" in result.stderr


def test_missing_meta_study_error(tmp_path):
    d1 = make_study(tmp_path, "b1", [("P1", "S1")])
    d2 = make_study(tmp_path, "b2", [("P2", "S2")])
    (d2 / "meta_study.txt").unlink()
    result = run([
        "--input_dir_1", str(d1), "--input_dir_2", str(d2),
        "--output_dir", str(tmp_path / "merged"),
    ])
    assert result.returncode != 0
    assert "meta_study.txt not found" in result.stderr


def test_study_id_mismatch_error(tmp_path):
    d1 = make_study(tmp_path, "b1", [("P1", "S1")], study_id="study_a")
    d2 = make_study(tmp_path, "b2", [("P2", "S2")], study_id="study_b")
    result = run([
        "--input_dir_1", str(d1), "--input_dir_2", str(d2),
        "--output_dir", str(tmp_path / "merged"),
    ])
    assert result.returncode != 0
    assert "cancer_study_identifier mismatch" in result.stderr


def test_existing_output_dir_error(tmp_path):
    d1 = make_study(tmp_path, "b1", [("P1", "S1")])
    d2 = make_study(tmp_path, "b2", [("P2", "S2")])
    (tmp_path / "merged").mkdir()
    result = run([
        "--input_dir_1", str(d1), "--input_dir_2", str(d2),
        "--output_dir", str(tmp_path / "merged"),
    ])
    assert result.returncode != 0
    assert "output directory already exists" in result.stderr


def test_missing_input_dir_error(tmp_path):
    d1 = make_study(tmp_path, "b1", [("P1", "S1")])
    result = run([
        "--input_dir_1", str(d1), "--input_dir_2", str(tmp_path / "nope"),
        "--output_dir", str(tmp_path / "merged"),
    ])
    assert result.returncode != 0
    assert "does not exist" in result.stderr
