#!/usr/bin/env python3
"""Fill a mafsmith MAF's Entrez_Gene_Id from the Ensembl annotations TSV.

mafsmith writes Entrez_Gene_Id = 0 on every row: it carries no Ensembl-to-Entrez
table. Each row is looked up by its Ensembl `Gene` id first, then by `Hugo_Symbol`.
A row with no match, or whose key maps to more than one Entrez id, keeps its value,
and cBioPortal then resolves the gene by Hugo_Symbol.
"""

import argparse
import re
import sys

ENTREZ = re.compile(r"^(\d+)(?:\.0+)?$")  # BioMart exports can write ids as floats


def read_annotations(path):
    """Return ({ensembl_id: entrez}, {symbol: entrez}), keeping unambiguous keys only."""
    by_gene, by_symbol = {}, {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for col in ("ensembl_id", "gene_symbol", "entrez_ncbi_id"):
            if col not in idx:
                sys.exit(f"ERROR: annotations TSV is missing the {col} column")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(header):
                continue
            m = ENTREZ.match(f[idx["entrez_ncbi_id"]].strip())
            if not m:
                continue
            for table, k in ((by_gene, f[idx["ensembl_id"]].strip()),
                             (by_symbol, f[idx["gene_symbol"]].strip())):
                if k:
                    table.setdefault(k, set()).add(m.group(1))
    return unique(by_gene), unique(by_symbol)


def unique(table):
    """Drop keys that map to more than one Entrez id."""
    return {k: next(iter(v)) for k, v in table.items() if len(v) == 1}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--maf", required=True)
    p.add_argument("--annotations", required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    by_gene, by_symbol = read_annotations(args.annotations)

    total = filled = 0
    with open(args.maf) as fh, open(args.output, "w") as out:
        version = fh.readline()
        header_line = fh.readline()
        out.write(version)
        out.write(header_line)
        header = header_line.rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        if "Entrez_Gene_Id" not in idx:
            sys.exit("ERROR: MAF is missing the Entrez_Gene_Id column")
        gene_i, sym_i, ent_i = idx.get("Gene"), idx.get("Hugo_Symbol"), idx["Entrez_Gene_Id"]

        for line in fh:
            row = line.rstrip("\n").split("\t")
            total += 1
            entrez = None
            if gene_i is not None and gene_i < len(row):
                entrez = by_gene.get(row[gene_i])
            if entrez is None and sym_i is not None and sym_i < len(row):
                entrez = by_symbol.get(row[sym_i])
            if entrez is not None:
                row[ent_i] = entrez
                filled += 1
            out.write("\t".join(row) + "\n")

    print(f"Filled Entrez_Gene_Id on {filled}/{total} MAF rows", file=sys.stderr)


if __name__ == "__main__":
    main()
