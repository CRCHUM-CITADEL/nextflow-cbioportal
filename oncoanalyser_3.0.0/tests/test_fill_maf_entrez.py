"""Tests for bin/fill_maf_entrez.py."""

import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).parent.parent / "bin" / "fill_maf_entrez.py"

ANNOT_HEADER = ["ensembl_id", "gene_symbol", "entrez_ncbi_id", "biotype"]
MAF_COLS = ["Hugo_Symbol", "Entrez_Gene_Id", "Chromosome", "Gene"]


def run(tmp_path, annot_rows, maf_rows, maf_cols=MAF_COLS):
    annot, maf, out = tmp_path / "annot.tsv", tmp_path / "in.maf", tmp_path / "out.maf"
    annot.write_text("\n".join("\t".join(r) for r in [ANNOT_HEADER] + annot_rows) + "\n")
    lines = ["#version 2.4", "\t".join(maf_cols)] + ["\t".join(r) for r in maf_rows]
    maf.write_text("\n".join(lines) + "\n")
    res = subprocess.run(
        [sys.executable, str(SCRIPT), "--maf", str(maf), "--annotations", str(annot),
         "--output", str(out)],
        capture_output=True, text=True,
    )
    return res, out


def entrez(out):
    lines = out.read_text().rstrip("\n").split("\n")
    header = lines[1].split("\t")
    return [dict(zip(header, ln.split("\t")))["Entrez_Gene_Id"] for ln in lines[2:]]


ANNOT = [
    ["ENSG00000012048", "BRCA1", "672", "protein_coding"],
    ["ENSG00000141510", "TP53", "7157.0", "protein_coding"],  # float-formatted id
    ["ENSG00000000001", "DUPSYM", "111", "protein_coding"],
    ["ENSG00000000002", "DUPSYM", "222", "protein_coding"],   # symbol maps to two ids
    ["ENSG00000000003", "NOID", "", "lncRNA"],
    ["ENSG00000000004", "NAID", "NA", "lncRNA"],
]


def test_filled_by_ensembl_gene_id(tmp_path):
    res, out = run(tmp_path, ANNOT, [["BRCA1", "0", "17", "ENSG00000012048"]])
    assert res.returncode == 0, res.stderr
    assert entrez(out) == ["672"]


def test_float_formatted_id_is_written_as_integer(tmp_path):
    res, out = run(tmp_path, ANNOT, [["TP53", "0", "17", "ENSG00000141510"]])
    assert res.returncode == 0, res.stderr
    assert entrez(out) == ["7157"]


def test_falls_back_to_hugo_symbol_without_gene_id(tmp_path):
    res, out = run(tmp_path, ANNOT, [["BRCA1", "0", "17", ""]])
    assert res.returncode == 0, res.stderr
    assert entrez(out) == ["672"]


def test_gene_id_wins_over_symbol(tmp_path):
    res, out = run(tmp_path, ANNOT, [["TP53", "0", "17", "ENSG00000012048"]])
    assert res.returncode == 0, res.stderr
    assert entrez(out) == ["672"]


def test_ambiguous_symbol_is_left_alone(tmp_path):
    res, out = run(tmp_path, ANNOT, [["DUPSYM", "0", "1", ""]])
    assert res.returncode == 0, res.stderr
    assert entrez(out) == ["0"]


def test_ambiguous_symbol_still_resolves_by_gene_id(tmp_path):
    res, out = run(tmp_path, ANNOT, [["DUPSYM", "0", "1", "ENSG00000000002"]])
    assert res.returncode == 0, res.stderr
    assert entrez(out) == ["222"]


def test_unmapped_rows_keep_their_value(tmp_path):
    rows = [
        ["Unknown", "0", "1", ""],                 # intergenic
        ["NOID", "0", "1", "ENSG00000000003"],     # gene has no Entrez id
        ["NAID", "0", "1", "ENSG00000000004"],     # Entrez id written as NA
    ]
    res, out = run(tmp_path, ANNOT, rows)
    assert res.returncode == 0, res.stderr
    assert entrez(out) == ["0", "0", "0"]


def test_other_columns_and_width_untouched(tmp_path):
    cols = MAF_COLS + ["HMF.TIER", "RNA.vaf"]
    rows = [["BRCA1", "0", "17", "ENSG00000012048", "", "0.25"]]
    res, out = run(tmp_path, ANNOT, rows, maf_cols=cols)
    assert res.returncode == 0, res.stderr
    lines = out.read_text().rstrip("\n").split("\n")
    assert lines[0] == "#version 2.4"
    assert lines[1].split("\t") == cols
    assert lines[2].split("\t") == ["BRCA1", "672", "17", "ENSG00000012048", "", "0.25"]


def test_works_without_gene_column(tmp_path):
    cols = ["Hugo_Symbol", "Entrez_Gene_Id"]
    res, out = run(tmp_path, ANNOT, [["BRCA1", "0"]], maf_cols=cols)
    assert res.returncode == 0, res.stderr
    assert entrez(out) == ["672"]


def test_missing_entrez_column_fails(tmp_path):
    res, _ = run(tmp_path, ANNOT, [["BRCA1", "17"]], maf_cols=["Hugo_Symbol", "Chromosome"])
    assert res.returncode != 0
    assert "Entrez_Gene_Id" in res.stderr
