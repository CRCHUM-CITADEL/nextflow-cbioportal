#!/usr/bin/env python3
"""Usage: build_ensembl_annotations.py <release>.gff3.gz <release>.entrez.tsv.gz <out.tsv> (kit tool: BioMart-layout ensembl_annotations)"""
import csv
import gzip
import sys
from collections import defaultdict
from urllib.parse import unquote

HEADER = ["ensembl_id", "gene_symbol", "entrez_ncbi_id", "biotype", "chr", "start", "end",
          "strand", "description"]


def read_entrez(path):
    entrez = defaultdict(set)
    with gzip.open(path, "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["db_name"] == "EntrezGene":
                entrez[row["gene_stable_id"]].add(row["xref"])
    return entrez


def gene_rows(gff3, entrez):
    with gzip.open(gff3, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            # Gene records only (gene / ncRNA_gene / pseudogene); *_gene_segment rows are transcripts.
            if not attrs.get("ID", "").startswith("gene:"):
                continue
            gene_id = attrs["gene_id"]
            for eid in sorted(entrez.get(gene_id, {""}), key=lambda x: (len(x), x)):
                yield [gene_id, unquote(attrs.get("Name", "")), eid, attrs.get("biotype", ""),
                       f[0], f[3], f[4], "1" if f[6] == "+" else "-1",
                       unquote(attrs.get("description", ""))]


def main(gff3, entrez_tsv, out):
    # Sorted by ensembl_id (stable: a gene's Entrez ids keep their numeric order).
    rows = sorted(gene_rows(gff3, read_entrez(entrez_tsv)), key=lambda r: r[0])
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(HEADER)
        w.writerows(rows)
    print(f"{len(rows)} rows, {len({r[0] for r in rows})} genes, "
          f"{sum(1 for r in rows if r[2])} with Entrez", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:4])
