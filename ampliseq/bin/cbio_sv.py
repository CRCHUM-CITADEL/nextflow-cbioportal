"""
Shared cBioPortal structural-variant (data_sv.txt) column layout for the ampliseq SV writers
(format_tsv.py, fusion_vcf_to_sv.py).

Every per-sample _sv.txt is written with exactly SV_COLUMNS, in this order, so per-sample files
from either writer can be merged into one data_sv.txt. The set follows oncoanalyser's
gen_*_to_cbioportal.R writers plus the cBioPortal fields ampliseq inputs can populate
(Length, Site1/2_Description, External_Annotation).
"""

import csv
import os

SV_COLUMNS = [
    "Sample_Id",
    "SV_Status",
    "Site1_Hugo_Symbol",
    "Site1_Ensembl_Transcript_Id",
    "Site1_Region_Number",
    "Site1_Region",
    "Site1_Chromosome",
    "Site1_Position",
    "Site1_Description",
    "Site2_Hugo_Symbol",
    "Site2_Ensembl_Transcript_Id",
    "Site2_Region_Number",
    "Site2_Region",
    "Site2_Chromosome",
    "Site2_Position",
    "Site2_Description",
    "Site2_Effect_On_Frame",
    "NCBI_Build",
    "Class",
    "Length",
    "DNA_Support",
    "RNA_Support",
    "Tumor_Variant_Count",
    "Tumor_Split_Read_Count",
    "Tumor_Paired_End_Read_Count",
    "Tumor_Read_Count",
    "Connection_Type",
    "Breakpoint_Type",
    "Event_Info",
    "Annotation",
    "External_Annotation",
    "Comments",
]

NCBI_BUILD = "GRCh37"
MISSING = "NA"


def strip_chr(chrom):
    chrom = str(chrom)
    return chrom[3:] if chrom.startswith("chr") else chrom


def sv_length(chrom1, pos1, chrom2, pos2):
    """Distance between breakpoints on the same chromosome; NA for translocations."""
    try:
        if strip_chr(chrom1) != strip_chr(chrom2):
            return MISSING
        return abs(int(pos2) - int(pos1))
    except (TypeError, ValueError):
        return MISSING


def write_rows(rows, out_path="data_sv.txt"):
    """Append rows (dicts keyed by SV_COLUMNS) to out_path, writing the header only when the file
    is new or empty, so bin/run_pipeline.sh can call a writer once per sample into one file.
    Unknown keys are an error; missing or empty values are written as NA."""
    new_file = not os.path.exists(out_path) or os.path.getsize(out_path) == 0
    with open(out_path, "a", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=SV_COLUMNS, delimiter="\t",
                                restval=MISSING, extrasaction="raise", lineterminator="\n")
        if new_file:
            writer.writeheader()
        for row in rows:
            writer.writerow({k: (MISSING if v is None or v == "" else v) for k, v in row.items()})
    print(f"{'Written' if new_file else 'Appended'}: {os.path.abspath(out_path)} ({len(rows)} rows)")
