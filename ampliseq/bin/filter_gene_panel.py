#!/usr/bin/env python3
"""
Keep only rows whose gene is on the gene panel for that data type.

Usage: filter_gene_panel.py <data_file> --panel <gene_panel.tsv> --type mutation|cna|sv --out <out>

Gene panel TSV: a `gene` and a `type` column (mutation, cna, sv; `fusion` = sv), one row per
gene and type, `#` lines ignored. Genes and types are compared case-insensitively.

- mutation / cna: the row's Hugo_Symbol must be on the panel.
- sv: Site1_Hugo_Symbol or Site2_Hugo_Symbol must be on the panel (a fusion's partner is often
  off-panel, e.g. EML4 in EML4-ALK).

The panel must list at least one gene for the requested type, so a typo cannot silently empty a
data type. Rows are written byte-for-byte; lines before the header starting with `#` are kept.
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


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("data_file")
    parser.add_argument("--panel", required=True)
    parser.add_argument("--type", required=True, choices=sorted(TYPES))
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    kind = TYPES[args.type]
    genes = read_panel(args.panel)[kind]
    if not genes:
        sys.exit(f"ERROR: {args.panel} lists no '{kind}' genes; refusing to drop every {kind} row")

    kept = dropped = 0
    dropped_genes = collections.Counter()
    with open(args.data_file) as fh, open(args.out, "w") as out:
        idx = None
        for line in fh:
            if idx is None:
                out.write(line)
                if line.startswith("#") or not line.strip():
                    continue
                header = line.rstrip("\r\n").split("\t")
                missing = [c for c in GENE_COLUMNS[kind] if c not in header]
                if missing:
                    sys.exit(f"ERROR: {args.data_file} has no {', '.join(missing)} column")
                idx = [header.index(c) for c in GENE_COLUMNS[kind]]
                continue
            cols = line.rstrip("\r\n").split("\t")
            row_genes = [cols[i].strip().upper() if i < len(cols) else "" for i in idx]
            if any(g in genes for g in row_genes):
                out.write(line)
                kept += 1
            else:
                dropped += 1
                dropped_genes["/".join(g or "-" for g in row_genes)] += 1

    print(f"gene panel ({kind}, {len(genes)} genes): kept {kept}, dropped {dropped} rows", file=sys.stderr)
    for gene, n in dropped_genes.most_common(20):
        print(f"  dropped {gene}: {n}", file=sys.stderr)


if __name__ == "__main__":
    main()
