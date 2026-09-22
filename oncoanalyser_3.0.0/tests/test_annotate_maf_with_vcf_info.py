"""Tests for bin/annotate_maf_with_vcf_info.py."""

import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).parent.parent / "bin" / "annotate_maf_with_vcf_info.py"

MAF_COLS = ["Hugo_Symbol", "Chromosome", "Start_Position",
            "Reference_Allele", "Tumor_Seq_Allele2"]


def run(args, **kwargs):
    return subprocess.run(
        [sys.executable, str(SCRIPT)] + args,
        capture_output=True, text=True, **kwargs,
    )


def write_vcf(path, records):
    """records: list of (chrom, pos, ref, alt, info)."""
    lines = ["##fileformat=VCFv4.2",
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
    lines += [f"{c}\t{p}\t.\t{r}\t{a}\t.\tPASS\t{i}" for c, p, r, a, i in records]
    path.write_text("\n".join(lines) + "\n")


def write_maf(path, rows):
    """rows: list of (hugo, chrom, start, ref, alt)."""
    lines = ["#version 2.4", "\t".join(MAF_COLS)]
    lines += ["\t".join(str(v) for v in r) for r in rows]
    path.write_text("\n".join(lines) + "\n")


def read_maf(path):
    """Return (header, [dict per row]) preserving blanks."""
    lines = path.read_text().rstrip("\n").split("\n")
    header = lines[1].split("\t")
    return header, [dict(zip(header, ln.split("\t"))) for ln in lines[2:]]


def annotate(tmp_path, vcf_records, maf_rows, extra=None):
    vcf, maf, out = tmp_path / "in.vcf", tmp_path / "in.maf", tmp_path / "out.maf"
    write_vcf(vcf, vcf_records)
    write_maf(maf, maf_rows)
    res = run(["--vcf", str(vcf), "--maf", str(maf), "--output", str(out)] + (extra or []))
    return res, out


# ── Coordinate conventions ────────────────────────────────────────────────────
# A MAF does not use VCF coordinates for indels. Keying the join on raw POS/REF/ALT
# would silently annotate only SNVs, so each class gets its own case.

def test_snv_is_annotated(tmp_path):
    res, out = annotate(
        tmp_path,
        [("chr1", 609395, "G", "A", "GND_FREQ=6.670e-03;TIER=LOW_CONFIDENCE;TNC=CGG")],
        [("GENE1", "1", 609395, "G", "A")],
    )
    assert res.returncode == 0, res.stderr
    _, rows = read_maf(out)
    assert rows[0]["gnomAD_AF"] == "6.670e-03"
    assert rows[0]["HMF.TIER"] == "LOW_CONFIDENCE"
    assert rows[0]["HMF.TNC"] == "CGG"


def test_deletion_is_annotated(tmp_path):
    # VCF 700000 AT>A  ->  MAF Start 700001, T > -
    res, out = annotate(
        tmp_path,
        [("chr1", 700000, "AT", "A", "GND_FREQ=1.2e-02")],
        [("GENE2", "1", 700001, "T", "-")],
    )
    assert res.returncode == 0, res.stderr
    _, rows = read_maf(out)
    assert rows[0]["gnomAD_AF"] == "1.2e-02"


def test_insertion_is_annotated(tmp_path):
    # VCF 800000 A>ATT  ->  MAF Start 800000 (flanking base), - > TT
    res, out = annotate(
        tmp_path,
        [("chr1", 800000, "A", "ATT", "GND_FREQ=3.4e-02")],
        [("GENE3", "1", 800000, "-", "TT")],
    )
    assert res.returncode == 0, res.stderr
    _, rows = read_maf(out)
    assert rows[0]["gnomAD_AF"] == "3.4e-02"


def test_multiallelic_alts_are_both_resolvable(tmp_path):
    # mafsmith picks one "effective" alt from a multi-allelic record; the index must
    # carry every alt so the join works whichever one it chose.
    res, out = annotate(
        tmp_path,
        [("chr1", 900000, "G", "A,T", "TIER=HOTSPOT")],
        [("GENE4", "1", 900000, "G", "T")],
    )
    assert res.returncode == 0, res.stderr
    _, rows = read_maf(out)
    assert rows[0]["HMF.TIER"] == "HOTSPOT"


def test_chr_prefix_mismatch_still_joins(tmp_path):
    # The MAF has its 'chr' prefix stripped downstream; both sides are normalized.
    res, out = annotate(
        tmp_path,
        [("chr7", 100, "C", "T", "TIER=PANEL")],
        [("GENE5", "chr7", 100, "C", "T")],
    )
    assert res.returncode == 0, res.stderr
    _, rows = read_maf(out)
    assert rows[0]["HMF.TIER"] == "PANEL"


# ── Width invariants ──────────────────────────────────────────────────────────
# Per-subject MAFs are merged by a group-level collectFile. A column set that varies
# with what a given VCF happened to declare would produce a ragged merged file.

def test_unmatched_row_gets_blanks_and_keeps_width(tmp_path):
    res, out = annotate(
        tmp_path,
        [("chr1", 100, "A", "G", "TIER=HOTSPOT")],
        [("HIT", "1", 100, "A", "G"), ("MISS", "1", 999, "C", "T")],
    )
    assert res.returncode == 0, res.stderr
    header, rows = read_maf(out)
    assert rows[0]["HMF.TIER"] == "HOTSPOT"
    assert rows[1]["HMF.TIER"] == ""
    assert all(len(r) == len(header) for r in rows)


def test_columns_emitted_even_when_vcf_declares_none(tmp_path):
    res, out = annotate(
        tmp_path,
        [("chr1", 100, "A", "G", "SOMETHING_ELSE=1")],
        [("GENE", "1", 100, "A", "G")],
    )
    assert res.returncode == 0, res.stderr
    header, rows = read_maf(out)
    assert "gnomAD_AF" in header and "HMF.TIER" in header
    assert rows[0]["gnomAD_AF"] == ""
    assert len(rows[0]) == len(header)


def test_flag_info_field_reads_as_present(tmp_path):
    res, out = annotate(
        tmp_path,
        [("chr1", 100, "A", "G", "MSG")],
        [("GENE", "1", 100, "A", "G")],
    )
    assert res.returncode == 0, res.stderr
    _, rows = read_maf(out)
    assert rows[0]["HMF.MSG"] == "1"


# ── Failure modes ─────────────────────────────────────────────────────────────

def test_zero_matches_is_a_hard_error(tmp_path):
    # A silently 0% join is the failure this guards: it means the key is wrong.
    res, _ = annotate(
        tmp_path,
        [("chr1", 100, "A", "G", "TIER=HOTSPOT")],
        [("GENE", "1", 500, "C", "T")],
    )
    assert res.returncode != 0
    assert "no MAF row matched" in res.stderr


def test_min_match_rate_is_enforced(tmp_path):
    res, _ = annotate(
        tmp_path,
        [("chr1", 100, "A", "G", "TIER=HOTSPOT")],
        [("HIT", "1", 100, "A", "G"), ("MISS", "1", 999, "C", "T")],
        extra=["--min-match-rate", "0.9"],
    )
    assert res.returncode != 0
    assert "below --min-match-rate" in res.stderr


def test_missing_required_maf_column_errors(tmp_path):
    vcf, maf, out = tmp_path / "in.vcf", tmp_path / "in.maf", tmp_path / "out.maf"
    write_vcf(vcf, [("chr1", 100, "A", "G", "TIER=HOTSPOT")])
    maf.write_text("#version 2.4\nHugo_Symbol\tChromosome\nGENE\t1\n")
    res = run(["--vcf", str(vcf), "--maf", str(maf), "--output", str(out)])
    assert res.returncode != 0
    assert "missing the Start_Position column" in res.stderr
