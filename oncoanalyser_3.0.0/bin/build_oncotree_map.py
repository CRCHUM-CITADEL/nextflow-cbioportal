#!/usr/bin/env python3
"""Build the kit's ICD -> OncoTree table.

Rows are keyed on (ICD_CODE, ICDO_MORPHOLOGY). ICD_CODE is a full code (C25.0) or its 3-character site (C25).
A morphology row takes, in order: our overrides (assets/oncotree_icdo_overrides.tsv), OncoTree's own ICD-O
mapping (ontology_mappings.txt), then the skfit-uni-luebeck ICD-O ConceptMap (a draft that often disagrees
with OncoTree's; its unknown-primary targets at known sites are lifted to a known-site ancestor or dropped).
A "*" row (any morphology) comes from the skfit-uni-luebeck ICD-10 ConceptMap, else, for a site, from the
site's morphology rows. Several OncoTree codes on one key resolve to their lowest common
ancestor; a key that only resolves to the root (TISSUE) is left out. Targets absent from the OncoTree version
are dropped, except in the overrides, where they fail the build.
"""

import argparse
import csv
import json
import sys
from collections import defaultdict

ROOT = "TISSUE"
# Ill-defined sites, lymph nodes and unknown primary: the only topographies where an unknown-primary type fits
UNKNOWN_PRIMARY_SITES = ("C76", "C77", "C80")
COLUMNS = ["ICD_CODE", "ICDO_MORPHOLOGY", "ONCOTREE_CODE", "CANCER_TYPE", "CANCER_TYPE_DETAILED", "SOURCE"]


def lineage(code, tree):
    """Code, its parent, ... up to the root."""
    chain = []
    while code in tree:
        chain.append(code)
        code = tree[code]["parent"]
    return chain + [ROOT]


def common_ancestor(codes, tree):
    lineages = [lineage(c, tree) for c in sorted(codes)]
    shared = set(lineages[0]).intersection(*lineages[1:])
    return next(c for c in lineages[0] if c in shared)


def labels(code, tree):
    node = tree[code]
    # A tissue node's name ("Pancreas") is not a cancer type; its main type ("Pancreatic Cancer") is.
    detailed = node["mainType"] if node["level"] <= 1 else node["name"]
    return node["mainType"], detailed


def concept_map(path, tree):
    """(source code, OncoTree target) for every matched element whose target is in the tree."""
    with open(path) as fh:
        groups = json.load(fh)["group"]
    for group in groups:
        for element in group["element"]:
            for target in element.get("target", []):
                code = target["code"].strip()
                if target.get("equivalence") == "unmatched" or code == ROOT:
                    continue
                if code not in tree:
                    print(f"skipping {element['code']!r} -> {code}: not in this OncoTree version", file=sys.stderr)
                    continue
                yield element["code"].strip().upper(), code


def known_site_code(code, tree):
    """Nearest ancestor that is not a Cancer of Unknown Primary type (HNSCUP -> HNSC), or None."""
    for c in lineage(code, tree):
        if c in (ROOT, "OTHER"):
            return None
        if "CUP" not in lineage(c, tree) and "unknown primary" not in tree[c]["name"].lower():
            return c
    return None


def luebeck_pairs(path, tree):
    """(code or site, morphology) -> OncoTree codes from the ICD-O ConceptMap.
    It often maps tumours of a known site to an unknown-primary type: those take known_site_code()."""
    pairs = defaultdict(set)
    for pair, code in concept_map(path, tree):
        parts = pair.split()
        if len(parts) != 2:
            continue
        if parts[0][:3] not in UNKNOWN_PRIMARY_SITES:
            code = known_site_code(code, tree)
        if code:
            for key in {parts[0], parts[0][:3]}:
                pairs[(key, parts[1])].add(code)
    return pairs


def oncotree_pairs(path, tree):
    """(topography, morphology) -> OncoTree codes from OncoTree's ontology_mappings.txt."""
    pairs = defaultdict(set)
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            code = row["ONCOTREE_CODE"].strip()
            topography = row["ICDO_TOPOGRAPHY_CODE"].strip().upper()
            morphology = row["ICDO_MORPHOLOGY_CODE"].strip()
            if topography and morphology and code in tree:
                pairs[(topography, morphology)].add(code)
    return pairs


def oncotree_codes(key, morphology, pairs, specific):
    """OncoTree's codes for a key: its own row, else its site's NOS (.9) row, unless the key is a subsite
    OncoTree maps on its own (C18.1 appendix must not take C18.9 colon rows)."""
    if (key, morphology) in pairs:
        return pairs[(key, morphology)]
    if key in specific:
        return None
    return pairs.get((key[:3] + ".9", morphology))


def read_overrides(path, tree):
    """(code or site, morphology) -> OncoTree code; a site row covers every subsite."""
    overrides = {}
    with open(path) as fh:
        for row in csv.DictReader((line for line in fh if not line.startswith("#")), delimiter="\t"):
            key = (row["ICD_CODE"].strip().upper(), row["ICDO_MORPHOLOGY"].strip())
            code = row["ONCOTREE_CODE"].strip()
            if code not in tree:
                sys.exit(f"ERROR: override {key} -> {code}: not in this OncoTree version")
            overrides[key] = code
    return overrides


def build(icdo_path, icd10_path, mappings_path, tree_path, overrides_path):
    with open(tree_path) as fh:
        tree = {t["code"]: t for t in json.load(fh)}

    luebeck = luebeck_pairs(icdo_path, tree)
    official = oncotree_pairs(mappings_path, tree)
    overrides = read_overrides(overrides_path, tree)
    specific = {t for t, _ in official if t[:3] != t and not t.endswith(".9")}

    keys = set(luebeck) | set(overrides)
    for topography, morphology in official:
        keys |= {(topography, morphology), (topography[:3], morphology)}
    keys |= {(key, m) for key, _ in luebeck for site, m in overrides if len(site) == 3 and key[:3] == site}

    resolved = {}
    for key, morphology in keys:
        override = overrides.get((key, morphology)) or overrides.get((key[:3], morphology))
        codes = oncotree_codes(key, morphology, official, specific)
        if override:
            resolved[(key, morphology)] = ({override}, "override")
        elif codes:
            resolved[(key, morphology)] = (codes, "oncotree")
        elif (key, morphology) in luebeck:
            resolved[(key, morphology)] = (luebeck[(key, morphology)], "icd-o")

    site_codes = defaultdict(set)
    for (key, _), (codes, _) in resolved.items():
        site_codes[key[:3]] |= codes

    icd10 = defaultdict(set)
    for icd_code, code in concept_map(icd10_path, tree):
        for key in {icd_code, icd_code[:3]}:
            icd10[key].add(code)
    resolved.update({(key, "*"): (codes, "icd-10") for key, codes in icd10.items()})
    for site, codes in site_codes.items():
        resolved.setdefault((site, "*"), (codes, "icd-o"))

    rows = []
    for (icd_code, morphology), (codes, source) in resolved.items():
        code = common_ancestor(codes, tree)
        if code != ROOT:
            rows.append([icd_code, morphology, code, *labels(code, tree), source])
    return sorted(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--icdo", required=True, help="skfit-uni-luebeck concept-map-icd-o-oncotree.json")
    parser.add_argument("--icd10", required=True, help="skfit-uni-luebeck concept-map-icd-10-oncotree.json")
    parser.add_argument("--mappings", required=True,
                        help="OncoTree scripts/ontology_to_ontology_mapping_tool/ontology_mappings.txt")
    parser.add_argument("--tumor-types", required=True, help="OncoTree /api/tumorTypes?version=<version> JSON")
    parser.add_argument("--overrides", required=True, help="assets/oncotree_icdo_overrides.tsv")
    parser.add_argument("-o", "--output", required=True, help="Output TSV")
    args = parser.parse_args()

    rows = build(args.icdo, args.icd10, args.mappings, args.tumor_types, args.overrides)
    with open(args.output, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(COLUMNS)
        writer.writerows(rows)
    print(f"wrote {len(rows)} rows to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
