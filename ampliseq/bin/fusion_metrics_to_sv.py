#!/usr/bin/env python3
"""
Format one sample's fusions from the cohort-wide fusion-metrics TSV into data_sv.txt (cBioPortal SV
format, full column set from cbio_sv.SV_COLUMNS).
Usage: python3 fusion_metrics_to_sv.py <metrics_tsv> <sample_id> [--min-total-reads N] [--min-zr-score Z]

Layout: `###` preamble lines, then a `#Sample` header line (#Sample, Source, File extension(s), Fusion,
GeneA, GeneB, ExonA, ExonB, Transcript_idA, Transcript_idB, BreakpointA, BreakpointB, Junction Read Count,
DepthA, DepthB, Total Reads, Z-Score, R-score, ZR-score), one row per fusion. Rows are matched to the
sample on `Sample`. Breakpoints are `[chr]<chrom>:<pos>[:strand]`.

Filter mode: keep a fusion when `Total Reads >= N` and `ZR-score >= Z`; a non-numeric value fails the
filter. Every other metric is ignored.
"""

import argparse
import io
import re
import sys

import pandas as pd

from cbio_sv import NCBI_BUILD, sv_length, strip_chr, write_rows

BREAKPOINT = re.compile(r"^(?:chr)?([0-9XYMT]+):(\d+)", re.IGNORECASE)
NEEDED = ["Sample", "GeneA", "GeneB", "BreakpointA", "BreakpointB", "Total Reads", "ZR-score"]


def read_metrics(path):
    """The metrics table, `###` preamble skipped and the `#` stripped from the header."""
    with open(path) as fh:
        lines = [line for line in fh if not line.startswith("###")]
    if not lines or not lines[0].startswith("#"):
        sys.exit(f"ERROR: {path} has no '#Sample' header line after its '###' preamble")
    lines[0] = lines[0][1:]
    df = pd.read_csv(io.StringIO("".join(lines)), sep="\t", dtype=str, keep_default_na=False)
    df.columns = df.columns.str.strip()
    # The header has been seen with BreakpointA twice; pandas renames the second one
    if "BreakpointB" not in df.columns and "BreakpointA.1" in df.columns:
        print("WARNING: no BreakpointB column; using the second BreakpointA as BreakpointB", file=sys.stderr)
        df = df.rename(columns={"BreakpointA.1": "BreakpointB"})
    missing = [c for c in NEEDED if c not in df.columns]
    if missing:
        sys.exit(f"ERROR: {path} lacks column(s) {', '.join(missing)}")
    return df


def breakpoint(value):
    """(chrom, pos) of a `chr7:55249071[:+]` breakpoint, (None, None) when unparseable."""
    match = BREAKPOINT.match(str(value).strip())
    if not match:
        return None, None
    return strip_chr(match.group(1)), int(match.group(2))


def clean(value):
    """Stripped value, '' for a missing-value marker."""
    value = value.strip() if isinstance(value, str) else ""
    return "" if value in (".", "N/A", "NA") else value


def region_number(value):
    return clean(value) if clean(value).isdigit() else ""


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("metrics_tsv")
    parser.add_argument("sample_id")
    parser.add_argument("--min-total-reads", type=float, default=None)
    parser.add_argument("--min-zr-score", type=float, default=None)
    args = parser.parse_args()
    sample_id = args.sample_id

    df = read_metrics(args.metrics_tsv)
    df = df[df["Sample"].str.strip() == sample_id]
    if df.empty:
        print(f"WARNING: no fusions for {sample_id} in {args.metrics_tsv}", file=sys.stderr)

    total_reads = pd.to_numeric(df["Total Reads"], errors="coerce")
    if args.min_total_reads is not None:
        df = df[total_reads >= args.min_total_reads]
    if args.min_zr_score is not None:
        df = df[pd.to_numeric(df["ZR-score"], errors="coerce") >= args.min_zr_score]

    rows = []
    for _, r in df.iterrows():
        gene1, gene2 = clean(r["GeneA"]), clean(r["GeneB"])
        chrom1, pos1 = breakpoint(r["BreakpointA"])
        chrom2, pos2 = breakpoint(r["BreakpointB"])
        for label, value, chrom in (("BreakpointA", r["BreakpointA"], chrom1), ("BreakpointB", r["BreakpointB"], chrom2)):
            if chrom is None:
                print(f"WARNING: cannot parse {label} '{value}' of {gene1}--{gene2}; left NA", file=sys.stderr)
        fusion = clean(r.get("Fusion")) or f"{gene1}--{gene2}"
        rows.append({
            "Sample_Id":                   sample_id,
            "SV_Status":                   "SOMATIC",
            "Site1_Hugo_Symbol":           gene1,
            "Site1_Ensembl_Transcript_Id": clean(r.get("Transcript_idA")),
            "Site1_Region_Number":         region_number(r.get("ExonA")),
            "Site1_Chromosome":            chrom1,
            "Site1_Position":              pos1,
            "Site2_Hugo_Symbol":           gene2,
            "Site2_Ensembl_Transcript_Id": clean(r.get("Transcript_idB")),
            "Site2_Region_Number":         region_number(r.get("ExonB")),
            "Site2_Chromosome":            chrom2,
            "Site2_Position":              pos2,
            "NCBI_Build":                  NCBI_BUILD,
            "Class":                       "FUSION",
            "Length":                      sv_length(chrom1, pos1, chrom2, pos2),
            "DNA_Support":                 "No",
            "RNA_Support":                 "Yes",
            "Tumor_Variant_Count":         clean(r["Total Reads"]),
            "Connection_Type":             "5to3",
            "Event_Info":                  f"RNA-seq Fusion: {fusion}",
            "Annotation":                  f"{gene1} - {gene2} fusion" if gene1 and gene2 else "",
        })

    write_rows(rows)


if __name__ == "__main__":
    main()
