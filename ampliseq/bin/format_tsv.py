#!/usr/bin/env python3
"""
Format analysis_<number>_export.tsv, extracting CNV/SV rows into data_sv.txt.
Usage: python3 format_tsv.py <tsv_file> <sample_id>
"""

import sys
import os
import pandas as pd

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <tsv_file> <sample_id>", file=sys.stderr)
        sys.exit(1)

    tsv_file = sys.argv[1]
    sample_id = sys.argv[2]

    df = pd.read_csv(tsv_file, sep="\t", dtype=str)

    # Keep only CNV and SV rows
    sv_df = df[(df["Variant Type"].isin(["CNV", "SV"])) & (df["Variant Subtype"].isin(["FUSION"]))].copy()
    sv_df["_start"] = sv_df["Start"].astype(int)
    sv_df["_end"]   = sv_df["End"].astype(int)
    is_inversion = sv_df["_end"] < sv_df["_start"]

    out = pd.DataFrame({
        "Sample_Id":             sample_id,
        "SV_Status":             "Somatic",
        "Site1_Hugo_Symbol":     sv_df["Genes"],
        "Site1_Chromosome":      sv_df["Chr"],
        "Site1_Region":          sv_df["Start"],
        "Site2_Hugo_Symbol":     sv_df["Breakend Genes"],
        "Class":                 sv_df["Variant Subtype"].where(~is_inversion, "INVERSION"),
        "Tumor_Variant_Count":   sv_df["Supporting Reads"],
        "SV_Length":             (sv_df["_start"] - sv_df["_end"]).where(
                                     is_inversion,
                                     sv_df["_end"] - sv_df["_start"]
                                 ),
    })

    out_path = os.path.join(os.getcwd(), "data_sv.txt")
    write_header = not os.path.exists(out_path)
    out.to_csv(out_path, sep="\t", index=False, mode="a", header=write_header)
    print(f"{'Written' if write_header else 'Appended'}: {out_path} ({len(out)} rows)")

if __name__ == "__main__":
    main()
