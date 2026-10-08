#!/usr/bin/env python3
"""
Format analysis_<number>_export.tsv, extracting fusion rows into data_sv.txt (cBioPortal SV format,
full column set from cbio_sv.SV_COLUMNS).
Usage: python3 format_tsv_to_sv.py <tsv_file> <sample_id> [--gene-loci TSV] [--min-supporting-reads N]

Fusion rows (CNV/SV + FUSION): Site1 = `Genes` at `Chr`:`Start`, Site2 = `Breakend Genes` at `End`.
`End` is on the partner's chromosome, which the export omits, so it is looked up in --gene-loci.

Filter mode: keep all rows of a (Genes, Breakend Genes) pair whose summed Supporting Reads >= N.
"""

import argparse
import gzip
import os
import re
import sys

import pandas as pd

from cbio_sv import NCBI_BUILD, sv_length, strip_chr, write_rows

DEFAULT_LOCI = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "assets", "grch37_gene_loci.tsv.gz")
LOCUS_PAD = 10_000  # breakpoint may sit just outside the annotated gene span


def load_loci(path):
    """Hugo_Symbol -> [(chrom, start, end)] from a Hugo_Symbol/Chromosome/Start/End TSV (gzip ok)."""
    loci = {}
    with (gzip.open if path.endswith(".gz") else open)(path, "rt") as fh:
        next(fh)
        for line in fh:
            sym, chrom, start, end = line.rstrip("\n").split("\t")[:4]
            loci.setdefault(sym, []).append((strip_chr(chrom), int(start), int(end)))
    return loci


def partner_chrom(loci, genes, pos, chrom1):
    """Chromosome of the `genes` locus containing pos, else the gene's only chromosome, else None."""
    cands = [lc for g in re.split(r"[,;]", genes) for lc in loci.get(g.strip(), [])]
    hits = [c for c, s, e in cands if s - LOCUS_PAD <= pos <= e + LOCUS_PAD]
    if hits:
        return chrom1 if chrom1 in hits else hits[0]
    chroms = {c for c, _, _ in cands}
    return chroms.pop() if len(chroms) == 1 else None


def region_number(value):
    return value if isinstance(value, str) and value.strip().isdigit() else ""


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("tsv_file")
    parser.add_argument("sample_id")
    parser.add_argument("--gene-loci", default=DEFAULT_LOCI)
    parser.add_argument("--min-supporting-reads", type=int, default=None)
    args = parser.parse_args()
    sample_id = args.sample_id

    df = pd.read_csv(args.tsv_file, sep="\t", dtype=str)
    df.columns = df.columns.str.strip()
    needed = ["Chr", "Start", "End", "Variant Type", "Variant Subtype", "Genes", "Breakend Genes",
              "Supporting Reads"]
    missing = [c for c in needed if c not in df.columns]
    if missing:
        sys.exit(f"ERROR: {args.tsv_file} lacks column(s) {', '.join(missing)}")
    sv_df = df[df["Variant Type"].isin(["CNV", "SV"]) & df["Variant Subtype"].isin(["FUSION"])]

    if args.min_supporting_reads is not None:
        read_counts = pd.to_numeric(sv_df["Supporting Reads"], errors="coerce").fillna(0)
        pair_total = read_counts.groupby([sv_df["Genes"].fillna(""), sv_df["Breakend Genes"].fillna("")]) \
            .transform("sum")
        sv_df = sv_df[pair_total >= args.min_supporting_reads]

    loci = load_loci(args.gene_loci) if len(sv_df) else {}
    rows = []
    for _, r in sv_df.iterrows():
        start, end = int(r["Start"]), int(r["End"])
        gene1 = r["Genes"] if pd.notna(r["Genes"]) else ""
        gene2 = r["Breakend Genes"] if pd.notna(r["Breakend Genes"]) else ""
        chrom = strip_chr(r["Chr"])
        chrom2 = partner_chrom(loci, gene2, end, chrom)
        if chrom2 is None:
            print(f"WARNING: cannot place {gene2 or 'partner'} at {end}; Site2_Chromosome left NA", file=sys.stderr)
        reads = r["Supporting Reads"] if pd.notna(r["Supporting Reads"]) else ""
        rows.append({
            "Sample_Id":           sample_id,
            "SV_Status":           "SOMATIC",
            "Site1_Hugo_Symbol":   gene1,
            "Site1_Region_Number": region_number(r.get("Exons")),
            "Site1_Chromosome":    chrom,
            "Site1_Position":      start,
            "Site2_Hugo_Symbol":   gene2,
            "Site2_Region_Number": region_number(r.get("Breakend Exon")),
            "Site2_Chromosome":    chrom2,
            "Site2_Position":      end,
            "NCBI_Build":          NCBI_BUILD,
            "Class":               "FUSION",
            "Length":              sv_length(chrom, start, chrom2, end),
            "DNA_Support":         "No",
            "RNA_Support":         "Yes",
            "Tumor_Variant_Count": reads,
            "Connection_Type":     "5to3",
            "Event_Info":          f"RNA-seq Fusion: {gene1}--{gene2}",
            "Annotation":          f"{gene1} - {gene2} fusion" if gene1 and gene2 else "",
        })

    write_rows(rows)


if __name__ == "__main__":
    main()
