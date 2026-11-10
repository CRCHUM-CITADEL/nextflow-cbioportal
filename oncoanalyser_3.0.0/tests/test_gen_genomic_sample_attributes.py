"""Tests for bin/gen_genomic_sample_attributes.py."""

import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).parent.parent / "bin" / "gen_genomic_sample_attributes.py"

IDS = ["PATIENT_ID", "SAMPLE_ID", "TMB_PER_MB", "TMB_STATUS", "TML", "TML_STATUS",
       "HRD_SCORE", "HR_STATUS", "HRD_TYPE"]

PURITY = ("purity\tploidy\tgender\tstatus\tmsIndelsPerMb\tmsStatus\ttml\ttmlStatus\ttmbPerMb\ttmbStatus\n"
          "0.75\t2.0\tMALE\tPASS\t0.12\tMSS\t312\tHIGH\t12.5\tHIGH\n")
CHORD = ("sample\tBRCA1\tBRCA2\thrd\thr_status\thrd_type\tremarks_hr_status\tremarks_hrd_type\n"
         "T1-T\t0.02\t0.81\t0.83\tHR_deficient\tBRCA2_type\t\t\n")


def run(tmp_path, purity=None, chord=None):
    out = tmp_path / "out.tsv"
    args = [sys.executable, str(SCRIPT), "--patient", "P1", "--sample", "S1", "--output", str(out)]
    if purity is not None:
        (tmp_path / "purity.tsv").write_text(purity)
        args += ["--purity", str(tmp_path / "purity.tsv")]
    if chord is not None:
        (tmp_path / "chord.tsv").write_text(chord)
        args += ["--chord", str(tmp_path / "chord.tsv")]
    res = subprocess.run(args, capture_output=True, text=True)
    return res, out


def read(out):
    lines = out.read_text().rstrip("\n").split("\n")
    assert len(lines) == 6
    assert all(ln.startswith("#") for ln in lines[:4])
    header = lines[4].split("\t")
    assert all(len(ln.lstrip("#").split("\t")) == len(header) for ln in lines[:4])
    return lines, dict(zip(header, lines[5].split("\t")))


def test_purity_and_chord(tmp_path):
    res, out = run(tmp_path, PURITY, CHORD)
    assert res.returncode == 0, res.stderr
    lines, row = read(out)
    assert lines[4].split("\t") == IDS
    assert row == {"PATIENT_ID": "P1", "SAMPLE_ID": "S1", "TMB_PER_MB": "12.5", "TMB_STATUS": "HIGH",
                   "TML": "312", "TML_STATUS": "HIGH", "HRD_SCORE": "0.83",
                   "HR_STATUS": "HR_deficient", "HRD_TYPE": "BRCA2_type"}
    assert lines[2].lstrip("#").split("\t") == ["STRING", "STRING", "NUMBER", "STRING", "NUMBER",
                                                "STRING", "NUMBER", "STRING", "STRING"]


def test_missing_files_give_na(tmp_path):
    res, out = run(tmp_path, PURITY)
    assert res.returncode == 0, res.stderr
    _, row = read(out)
    assert row["TMB_PER_MB"] == "12.5"
    assert (row["HRD_SCORE"], row["HR_STATUS"], row["HRD_TYPE"]) == ("NA", "NA", "NA")

    res, out = run(tmp_path)
    assert res.returncode == 0, res.stderr
    _, row = read(out)
    assert all(row[c] == "NA" for c in IDS[2:])


def test_chord_p_hrd_header(tmp_path):
    # The CHORD README spells the probability column p_hrd.
    chord = CHORD.replace("\tBRCA1\tBRCA2\thrd\t", "\tp_BRCA1\tp_BRCA2\tp_hrd\t")
    res, out = run(tmp_path, chord=chord)
    assert res.returncode == 0, res.stderr
    assert read(out)[1]["HRD_SCORE"] == "0.83"


def test_empty_value_is_na(tmp_path):
    res, out = run(tmp_path, chord=CHORD.replace("BRCA2_type", ""))
    assert res.returncode == 0, res.stderr
    assert read(out)[1]["HRD_TYPE"] == "NA"


def test_purity_without_tmb_columns_fails(tmp_path):
    res, _ = run(tmp_path, "purity\tploidy\tgender\tstatus\n0.75\t2.0\tMALE\tPASS\n")
    assert res.returncode != 0
    assert "tmbPerMb" in res.stderr


def test_chord_without_hrd_column_fails(tmp_path):
    res, _ = run(tmp_path, chord="sample\thr_status\thrd_type\nT1-T\tHR_proficient\tnone\n")
    assert res.returncode != 0
    assert "hrd, p_hrd" in res.stderr
