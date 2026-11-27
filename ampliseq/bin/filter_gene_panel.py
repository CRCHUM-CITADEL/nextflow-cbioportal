#!/usr/bin/env python3
"""
Flag rows on/off the gene panel and drop rows of excluded genes.

Usage: filter_gene_panel.py <data_file> --type mutation|cna|sv --out <out> [--panel <gene_panel.tsv>] [--exclude G1,G2]

Gene panel TSV: a `gene` and a `type` column (mutation, cna, sv; `fusion` = sv), one row per
gene and type, `#` lines ignored. Genes and types are compared case-insensitively.

--panel adds (or overwrites) an `info.on_panel` column, true/false; no row is dropped for it:
- mutation / cna: the row's Hugo_Symbol is on the panel.
- sv: Site1_Hugo_Symbol and Site2_Hugo_Symbol are both on the panel (an empty partner is not).
The panel must list at least one gene for the requested type, so a typo cannot silently flag
every row false.

--exclude drops a row whose gene (sv: either partner) is listed, whatever the panel says.

An empty --panel / --exclude means unset; with neither, the file is copied unchanged. Lines
before the header starting with `#` are kept; rows are rewritten only to set `info.on_panel`.
"""

import argparse
import collections
import csv
import sys

TYPES = {"mutation": "mutation", "cna": "cna", "sv": "sv", "fusion": "sv"}
GENE_COLUMNS = {
    "mutation": ["Hugo_Symbol"],
    "cna": ["Hugo_Symbol"],
    "sv": ["Site1_Hugo_Symbol", "Site2_Hugo_Symbol"],
}
ON_PANEL_COLUMN = "info.on_panel"


def read_panel(path):
    panel = collections.defaultdict(set)
    with open(path, newline="") as fh:
        lines = (line for line in fh if line.strip() and not line.startswith("#"))
        rows = csv.reader(lines, delimiter="\t")
        header = [h.strip().lower() for h in next(rows, [])]
        if "gene" not in header or "type" not in header:
            sys.exit(f"ERROR: {path} needs a 'gene' and a 'type' column, got: {header}")
        gi, ti = header.index("gene"), header.index("type")
        for n, row in enumerate(rows, start=2):
            gene = row[gi].strip().upper() if gi < len(row) else ""
            kind = row[ti].strip().lower() if ti < len(row) else ""
            if not gene:
                continue
            if kind not in TYPES:
                sys.exit(f"ERROR: {path} row {n}: type '{kind}' is not one of {', '.join(TYPES)}")
            panel[TYPES[kind]].add(gene)
    return panel


def parse_genes(text):
    return {g.strip().upper() for g in (text or "").split(",") if g.strip()}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("data_file")
    parser.add_argument("--panel", default="")
    parser.add_argument("--exclude", default="")
    parser.add_argument("--type", required=True, choices=sorted(TYPES))
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    kind = TYPES[args.type]
    exclude = parse_genes(args.exclude)
    genes = None
    if args.panel:
        genes = read_panel(args.panel)[kind]
        if not genes:
            sys.exit(f"ERROR: {args.panel} lists no '{kind}' genes; refusing to flag every {kind} row off-panel")

    counts = collections.Counter()
    off_genes, excluded_genes = collections.Counter(), collections.Counter()
    with open(args.data_file) as fh, open(args.out, "w") as out:
        idx = flag = width = None
        for line in fh:
            if idx is None:
                if line.startswith("#") or not line.strip():
                    out.write(line)
                    continue
                header = line.rstrip("\r\n").split("\t")
                missing = [c for c in GENE_COLUMNS[kind] if c not in header]
                if missing:
                    sys.exit(f"ERROR: {args.data_file} has no {', '.join(missing)} column")
                idx = [header.index(c) for c in GENE_COLUMNS[kind]]
                if genes is not None:
                    if ON_PANEL_COLUMN not in header:
                        header.append(ON_PANEL_COLUMN)
                    flag, width = header.index(ON_PANEL_COLUMN), len(header)
                    out.write("\t".join(header) + "\n")
                else:
                    out.write(line)
                continue
            cols = line.rstrip("\r\n").split("\t")
            row_genes = [cols[i].strip().upper() if i < len(cols) else "" for i in idx]
            label = "/".join(g or "-" for g in row_genes)
            if any(g in exclude for g in row_genes if g):
                counts["excluded"] += 1
                excluded_genes[label] += 1
                continue
            if genes is None:
                out.write(line)
                counts["kept"] += 1
                continue
            on_panel = all(g in genes for g in row_genes)
            counts["on" if on_panel else "off"] += 1
            if not on_panel:
                off_genes[label] += 1
            cols += [""] * (width - len(cols))
            cols[flag] = "true" if on_panel else "false"
            out.write("\t".join(cols) + "\n")

    if genes is None:
        print(f"gene exclusion ({kind}): kept {counts['kept']}, excluded {counts['excluded']} rows", file=sys.stderr)
    else:
        print(
            f"gene panel ({kind}, {len(genes)} genes): on panel {counts['on']}, off panel {counts['off']}, "
            f"excluded {counts['excluded']} rows",
            file=sys.stderr,
        )
    for gene, n in excluded_genes.most_common(20):
        print(f"  excluded {gene}: {n}", file=sys.stderr)
    for gene, n in off_genes.most_common(20):
        print(f"  off panel {gene}: {n}", file=sys.stderr)


if __name__ == "__main__":
    main()
