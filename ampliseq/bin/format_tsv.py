#!/usr/bin/env python3
"""
Format analysis_<number>_export.tsv, extracting fusion rows into data_sv.txt (cBioPortal SV format,
full column set from cbio_sv.SV_COLUMNS).
Usage: python3 format_tsv.py <tsv_file> <sample_id>

Fusion rows are `Variant Type` CNV/SV with `Variant Subtype` FUSION. The export has one `Chr`
column, so both breakpoints are on that chromosome: Site1 = Start (`Genes`), Site2 = End
(`Breakend Genes`). End < Start marks an inversion.
"""

import sys

import pandas as pd

from cbio_sv import NCBI_BUILD, sv_length, strip_chr, write_rows


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <tsv_file> <sample_id>", file=sys.stderr)
        sys.exit(1)

    tsv_file = sys.argv[1]
    sample_id = sys.argv[2]

    df = pd.read_csv(tsv_file, sep="\t", dtype=str)
    sv_df = df[df["Variant Type"].isin(["CNV", "SV"]) & df["Variant Subtype"].isin(["FUSION"])]

    rows = []
    for _, r in sv_df.iterrows():
        start, end = int(r["Start"]), int(r["End"])
        gene1 = r["Genes"] if pd.notna(r["Genes"]) else ""
        gene2 = r["Breakend Genes"] if pd.notna(r["Breakend Genes"]) else ""
        chrom = strip_chr(r["Chr"])
        sv_class = "INVERSION" if end < start else "FUSION"
        reads = r["Supporting Reads"] if pd.notna(r["Supporting Reads"]) else ""
        rows.append({
            "Sample_Id":           sample_id,
            "SV_Status":           "SOMATIC",
            "Site1_Hugo_Symbol":   gene1,
            "Site1_Chromosome":    chrom,
            "Site1_Position":      start,
            "Site2_Hugo_Symbol":   gene2,
            "Site2_Chromosome":    chrom,
            "Site2_Position":      end,
            "NCBI_Build":          NCBI_BUILD,
            "Class":               sv_class,
            "Length":              sv_length(chrom, start, chrom, end),
            "DNA_Support":         "No",
            "RNA_Support":         "Yes",
            "Tumor_Variant_Count": reads,
            "Connection_Type":     "5to3",
            "Event_Info":          f"RNA-seq {sv_class.capitalize()}: {gene1}--{gene2}",
            "Annotation":          f"{gene1} - {gene2} fusion" if gene1 and gene2 else "",
        })

    write_rows(rows)


if __name__ == "__main__":
    main()
