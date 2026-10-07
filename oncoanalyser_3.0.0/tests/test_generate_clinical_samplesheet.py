"""Tests for bin/generate_clinical_samplesheet.py."""

import csv
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).parent.parent / "bin" / "generate_clinical_samplesheet.py"

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


def run(args, **kwargs):
    """Run the script and return CompletedProcess."""
    return subprocess.run(
        [sys.executable, str(SCRIPT)] + args,
        capture_output=True,
        text=True,
        **kwargs,
    )


def make_clinical_dir(tmp_path, filetypes=None):
    """Create a directory with stub clinical CSV files."""
    filetypes = filetypes if filetypes is not None else KNOWN_FILETYPES
    d = tmp_path / "clinical"
    d.mkdir()
    for ft in filetypes:
        (d / f"{ft}.csv").write_text(f"submitter_donor_id\nPATIENT1\n")
    return d


# ── Happy paths ───────────────────────────────────────────────────────────────

def test_all_filetypes_generated(tmp_path):
    """All 10 known files present → 10 rows, canonical order."""
    clinical_dir = make_clinical_dir(tmp_path)
    out = tmp_path / "samplesheet.csv"

    result = run([str(clinical_dir), "-o", str(out), "-d", "2025-05-20"])

    assert result.returncode == 0, result.stderr
    assert out.exists()

    with open(out) as f:
        rows = list(csv.DictReader(f))

    assert len(rows) == 10
    # Canonical order is preserved
    assert [r["filetype"] for r in rows] == KNOWN_FILETYPES


def test_partial_filetypes(tmp_path):
    """Only a subset of files present → only those rows written."""
    subset = ["donors", "primary_diagnoses", "sample_registrations"]
    clinical_dir = make_clinical_dir(tmp_path, filetypes=subset)
    out = tmp_path / "samplesheet.csv"

    result = run([str(clinical_dir), "-o", str(out), "-d", "2025-05-20"])

    assert result.returncode == 0
    with open(out) as f:
        rows = list(csv.DictReader(f))

    assert [r["filetype"] for r in rows] == subset


def test_extraction_date_written(tmp_path):
    """Specified date appears in every row."""
    clinical_dir = make_clinical_dir(tmp_path, filetypes=["donors"])
    out = tmp_path / "samplesheet.csv"

    run([str(clinical_dir), "-o", str(out), "-d", "2024-12-31"])

    with open(out) as f:
        rows = list(csv.DictReader(f))
    assert rows[0]["extraction_date"] == "2024-12-31"


def test_default_date_is_today(tmp_path):
    """Omitting --date uses today's date."""
    import datetime
    today = datetime.date.today().isoformat()
    clinical_dir = make_clinical_dir(tmp_path, filetypes=["donors"])
    out = tmp_path / "samplesheet.csv"

    run([str(clinical_dir), "-o", str(out)])

    with open(out) as f:
        rows = list(csv.DictReader(f))
    assert rows[0]["extraction_date"] == today


def test_study_id_written_to_group_id(tmp_path):
    """--study_id value appears in the group_id column."""
    clinical_dir = make_clinical_dir(tmp_path, filetypes=["donors"])
    out = tmp_path / "samplesheet.csv"

    run([str(clinical_dir), "-o", str(out), "-d", "2025-01-01", "--study_id", "MY_COHORT"])

    with open(out) as f:
        rows = list(csv.DictReader(f))
    assert rows[0]["group_id"] == "MY_COHORT"


def test_default_group_id_is_empty(tmp_path):
    """group_id is empty when --study_id is not given."""
    clinical_dir = make_clinical_dir(tmp_path, filetypes=["donors"])
    out = tmp_path / "samplesheet.csv"

    run([str(clinical_dir), "-o", str(out), "-d", "2025-01-01"])

    with open(out) as f:
        rows = list(csv.DictReader(f))
    assert rows[0]["group_id"] == ""


def test_absolute_paths_by_default(tmp_path):
    """Without --relative-to, file paths are absolute."""
    clinical_dir = make_clinical_dir(tmp_path, filetypes=["donors"])
    out = tmp_path / "samplesheet.csv"

    run([str(clinical_dir), "-o", str(out), "-d", "2025-01-01"])

    with open(out) as f:
        rows = list(csv.DictReader(f))
    assert Path(rows[0]["filepath"]).is_absolute()


def test_relative_to_flag(tmp_path):
    """--relative-to produces paths relative to the given base directory."""
    clinical_dir = make_clinical_dir(tmp_path, filetypes=["donors"])
    out = tmp_path / "samplesheet.csv"

    run([str(clinical_dir), "-o", str(out), "-d", "2025-01-01",
         "--relative-to", str(tmp_path)])

    with open(out) as f:
        rows = list(csv.DictReader(f))
    p = Path(rows[0]["filepath"])
    assert not p.is_absolute()
    assert p == Path("clinical") / "donors.csv"


def test_unknown_csv_files_ignored(tmp_path):
    """Unrecognised CSV files in the directory are silently skipped (warning only)."""
    clinical_dir = make_clinical_dir(tmp_path, filetypes=["donors"])
    (clinical_dir / "extra_data.csv").write_text("col\nval\n")
    out = tmp_path / "samplesheet.csv"

    result = run([str(clinical_dir), "-o", str(out), "-d", "2025-01-01"])

    assert result.returncode == 0
    assert "WARNING" in result.stderr
    with open(out) as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 1
    assert rows[0]["filetype"] == "donors"


def test_non_csv_files_ignored(tmp_path):
    """Non-CSV files are not included even if they match a filetype name."""
    clinical_dir = make_clinical_dir(tmp_path, filetypes=["donors"])
    (clinical_dir / "treatments.txt").write_text("ignored\n")
    out = tmp_path / "samplesheet.csv"

    run([str(clinical_dir), "-o", str(out), "-d", "2025-01-01"])

    with open(out) as f:
        rows = list(csv.DictReader(f))
    filetypes = [r["filetype"] for r in rows]
    assert "treatments" not in filetypes


def test_output_has_correct_columns(tmp_path):
    """Output CSV has exactly the expected header columns."""
    clinical_dir = make_clinical_dir(tmp_path, filetypes=["donors"])
    out = tmp_path / "samplesheet.csv"

    run([str(clinical_dir), "-o", str(out), "-d", "2025-01-01"])

    with open(out) as f:
        reader = csv.DictReader(f)
        assert reader.fieldnames == ["group_id", "filetype", "filepath", "extraction_date", "info"]


def test_info_column_is_empty(tmp_path):
    """The info column is always empty (the script does not populate it)."""
    clinical_dir = make_clinical_dir(tmp_path, filetypes=["donors"])
    out = tmp_path / "samplesheet.csv"

    run([str(clinical_dir), "-o", str(out), "-d", "2025-01-01"])

    with open(out) as f:
        rows = list(csv.DictReader(f))
    assert all(r["info"] == "" for r in rows)


def test_test_data_directory(tmp_path):
    """Run against the real test data directory and verify all 10 filetypes are found."""
    test_data_dir = (
        Path(__file__).parent.parent / "assets" / "test_data" / "clinical"
    )
    if not test_data_dir.exists():
        import pytest
        pytest.skip("test_data/clinical not found")

    out = tmp_path / "samplesheet.csv"
    result = run([str(test_data_dir), "-o", str(out), "-d", "2025-05-20"])

    assert result.returncode == 0
    with open(out) as f:
        rows = list(csv.DictReader(f))

    found_types = [r["filetype"] for r in rows]
    assert set(found_types) == set(KNOWN_FILETYPES)
    # Canonical order
    assert found_types == KNOWN_FILETYPES


# ── Error paths ───────────────────────────────────────────────────────────────

def test_missing_directory_exits(tmp_path):
    result = run([str(tmp_path / "nonexistent"), "-o", str(tmp_path / "out.csv")])
    assert result.returncode != 0
    assert "ERROR" in result.stderr


def test_empty_directory_exits(tmp_path):
    empty = tmp_path / "empty"
    empty.mkdir()
    result = run([str(empty), "-o", str(tmp_path / "out.csv")])
    assert result.returncode != 0
    assert "ERROR" in result.stderr


def test_invalid_date_exits(tmp_path):
    clinical_dir = make_clinical_dir(tmp_path, filetypes=["donors"])
    result = run([str(clinical_dir), "-o", str(tmp_path / "out.csv"), "-d", "not-a-date"])
    assert result.returncode != 0
    assert "ERROR" in result.stderr
