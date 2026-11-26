"""Tests for bin/build_oncotree_map.py."""

import csv
import json
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).parent.parent / "bin" / "build_oncotree_map.py"


def node(code, name, main_type, parent, level):
    return {"code": code, "name": name, "mainType": main_type, "parent": parent, "level": level}


TREE = [
    node("PANCREAS", "Pancreas", "Pancreatic Cancer", "TISSUE", 1),
    node("PAAD", "Pancreatic Adenocarcinoma", "Pancreatic Cancer", "PANCREAS", 2),
    node("PANET", "Pancreatic Neuroendocrine Tumor", "Pancreatic Cancer", "PANCREAS", 2),
    node("OVARY", "Ovary/Fallopian Tube", "Ovarian Cancer", "TISSUE", 1),
    node("OVT", "Ovarian Epithelial Tumor", "Ovarian Cancer", "OVARY", 2),
    node("SOC", "Serous Ovarian Cancer", "Ovarian Cancer", "OVT", 3),
    node("MXOV", "Mixed Ovarian Carcinoma", "Ovarian Cancer", "OVT", 3),
    node("BOWEL", "Bowel", "Bowel Cancer", "TISSUE", 1),
    node("COAD", "Colon Adenocarcinoma", "Colorectal Cancer", "BOWEL", 2),
    node("APAD", "Appendiceal Adenocarcinoma", "Appendiceal Cancer", "BOWEL", 2),
    node("HEAD_NECK", "Head and Neck", "Head and Neck Cancer", "TISSUE", 1),
    node("HNSC", "Head and Neck Squamous Cell Carcinoma", "Head and Neck Cancer", "HEAD_NECK", 2),
    node("HNSCUP", "Head and Neck Squamous Cell Carcinoma of Unknown Primary", "Head and Neck Cancer", "HNSC", 3),
    node("OTHER", "Other", "Other Cancer", "TISSUE", 1),
    node("CUP", "Cancer of Unknown Primary", "Cancer of Unknown Primary", "OTHER", 2),
    node("CUPNOS", "Cancer of Unknown Primary, NOS", "Cancer of Unknown Primary", "CUP", 3),
]


def concept_map(pairs):
    return {"group": [{"element": [
        {"code": src, "target": [{"code": tgt, "equivalence": eq}]} for src, tgt, eq in pairs
    ]}]}


ICDO = concept_map([
    ("C25.0 8140/3", "PANCREAS", "inexact"),     # OncoTree's C25.9 row wins for the subsite
    ("C25.0 8246/3", "PANET", "inexact"),        # no OncoTree row: Luebeck fills the gap
    ("C56.9 8441/3", "MXOV", "inexact"),         # OncoTree says SOC
    ("C56.0 8310/3", "MXOV", "inexact"),         # disagrees with C56.9 8310/3 below
    ("C56.9 8310/3", "SOC", "inexact"),
    ("C56.9 8480/3", "MXOV", "inexact"),
    ("C18.1 8480/3", "COAD", "inexact"),         # C18.1 has its own OncoTree rows: no C18.9 fallback
    ("C25.4 8151/3", "PANCREAS", "inexact"),     # overridden site-wide
    ("C09.0 8070/3", "HNSCUP", "inexact"),       # unknown primary at a known site: lifted to HNSC
    ("C34.0 8000/3", "CUPNOS", "inexact"),       # wholly unknown primary at a known site: dropped
    ("C80.9 8000/3", "CUPNOS", "inexact"),       # unknown primary site: kept
    ("C56.9 ", "OVARY", "inexact"),              # topography only: skipped
    ("C22.0 8170/3", "BOWEL", "unmatched"),      # unmatched: skipped
    ("C18.7 8140/3", "TISSUE", "inexact"),       # root: skipped
    ("C18.2 8140/3", "GONE", "inexact"),         # not in tree: skipped
    ("C16.0 8140/3", "PAAD", "inexact"),         # one site...
    ("C16.0 9999/3", "SOC", "inexact"),          # ...whose codes only share the root
])

ICD10 = concept_map([
    ("C25.0", "PANCREAS", "relatedto"),
    ("C25.1", "PAAD", "relatedto"),
    ("C18.7", "BOWEL", "relatedto"),
])

MAPPINGS = "\n".join([
    "ONCOTREE_CODE\tNCIT_CODE\tUMLS_CODE\tICDO_TOPOGRAPHY_CODE\tICDO_MORPHOLOGY_CODE\tHEMEONC_CODE",
    "PAAD\t\t\tC25.9\t8140/3\t",
    "SOC\t\t\tC56.9\t8441/3\t",
    "COAD\t\t\tC18.9\t8480/3\t",
    "APAD\t\t\tC18.1\t8140/3\t",
    "PANCREAS\t\t\t\t\t",
    "",
])

OVERRIDES = "\n".join([
    "# comment",
    "ICD_CODE\tICDO_MORPHOLOGY\tONCOTREE_CODE\tREASON",
    "C25\t8151/3\tPANET\tinsulinoma",
    "C56.9\t8480/3\tOVT\tcode only",
    "",
])


def build(tmp_path, overrides=OVERRIDES, check=True):
    files = {}
    for name, data in [("icdo", ICDO), ("icd10", ICD10), ("tree", TREE)]:
        files[name] = tmp_path / f"{name}.json"
        files[name].write_text(json.dumps(data))
    files["mappings"] = tmp_path / "ontology_mappings.txt"
    files["mappings"].write_text(MAPPINGS)
    files["overrides"] = tmp_path / "overrides.tsv"
    files["overrides"].write_text(overrides)
    out = tmp_path / "out.tsv"
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--icdo", str(files["icdo"]), "--icd10", str(files["icd10"]),
         "--mappings", str(files["mappings"]), "--tumor-types", str(files["tree"]),
         "--overrides", str(files["overrides"]), "-o", str(out)],
        capture_output=True, text=True, check=check,
    )
    if result.returncode:
        return None, result.stderr
    with open(out) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return {(r["ICD_CODE"], r["ICDO_MORPHOLOGY"]): r for r in rows}, result.stderr


def code(rows, key):
    return rows[key]["ONCOTREE_CODE"], rows[key]["SOURCE"]


def test_oncotree_rows_win_over_luebeck(tmp_path):
    rows, _ = build(tmp_path)
    assert code(rows, ("C56.9", "8441/3")) == ("SOC", "oncotree")
    assert code(rows, ("C56", "8441/3")) == ("SOC", "oncotree")
    # A subsite OncoTree doesn't map takes the site's NOS (.9) row
    assert code(rows, ("C25.0", "8140/3")) == ("PAAD", "oncotree")
    assert code(rows, ("C25.0", "8246/3")) == ("PANET", "icd-o")


def test_specific_subsite_does_not_take_nos_rows(tmp_path):
    rows, _ = build(tmp_path)
    assert code(rows, ("C18.1", "8480/3")) == ("COAD", "icd-o")
    assert code(rows, ("C18.9", "8480/3")) == ("COAD", "oncotree")


def test_overrides_win_and_a_site_override_covers_subsites(tmp_path):
    rows, _ = build(tmp_path)
    assert code(rows, ("C25", "8151/3")) == ("PANET", "override")
    assert code(rows, ("C25.4", "8151/3")) == ("PANET", "override")
    assert code(rows, ("C56.9", "8480/3")) == ("OVT", "override")
    assert rows[("C56.9", "8480/3")]["CANCER_TYPE_DETAILED"] == "Ovarian Epithelial Tumor"


def test_override_outside_the_tree_fails(tmp_path):
    rows, stderr = build(tmp_path, "ICD_CODE\tICDO_MORPHOLOGY\tONCOTREE_CODE\tREASON\nC25\t8151/3\tGONE\tx\n", check=False)
    assert rows is None
    assert "GONE" in stderr


def test_unknown_primary_targets_at_known_sites(tmp_path):
    rows, _ = build(tmp_path)
    assert code(rows, ("C09.0", "8070/3")) == ("HNSC", "icd-o")
    assert ("C34.0", "8000/3") not in rows
    assert code(rows, ("C80.9", "8000/3")) == ("CUPNOS", "icd-o")


def test_conflicting_codes_resolve_to_common_ancestor(tmp_path):
    rows, _ = build(tmp_path)
    # C56.0 -> MXOV and C56.9 -> SOC for 8310/3: the pooled site takes their parent
    assert code(rows, ("C56.9", "8310/3")) == ("SOC", "icd-o")
    assert code(rows, ("C56", "8310/3")) == ("OVT", "icd-o")
    # PAAD and SOC share only the root, so the site gets no row
    assert ("C16", "*") not in rows


def test_any_morphology_rows(tmp_path):
    rows, _ = build(tmp_path)
    assert code(rows, ("C25.1", "*")) == ("PAAD", "icd-10")
    assert code(rows, ("C25", "*")) == ("PANCREAS", "icd-10")
    assert code(rows, ("C18.7", "*")) == ("BOWEL", "icd-10")
    # No ICD-10 row for C56: the site's resolved morphology rows
    assert code(rows, ("C56", "*")) == ("OVT", "icd-o")


def test_tissue_node_reports_main_type_as_detailed(tmp_path):
    rows, _ = build(tmp_path)
    assert rows[("C25.0", "*")]["CANCER_TYPE"] == "Pancreatic Cancer"
    assert rows[("C25.0", "*")]["CANCER_TYPE_DETAILED"] == "Pancreatic Cancer"


def test_unusable_targets_are_dropped(tmp_path):
    rows, stderr = build(tmp_path)
    assert ("C56.9", "") not in rows
    assert ("C22.0", "8170/3") not in rows
    assert ("C18.7", "8140/3") not in rows
    assert ("C18.2", "8140/3") not in rows
    assert "GONE" in stderr
