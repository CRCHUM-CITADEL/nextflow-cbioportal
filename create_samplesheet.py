#!/usr/bin/env python3

"""
This script can be used to automatically create a samplesheet given a directory of standard
DRAGEN pipeline tumor pair output. (i.e., a dna/ folder with somatic and germinal sequencing files,
and a rna/ folder with tumor sequencing output)
"""

import os, re
import argparse
import pandas as pd
from pathlib import Path, PurePosixPath
from typing import List

DATAFRAME = {"group_id" : [],"subject_id": [],"sample_id" : [],"sample_type" : [],"sequence_data" : [],"filetype" : [],"info" : [],"filepath" : []}

def yield_hook(generator_output, function):
    for val in generator_output:
        yield function(val)

def capture_folders_by_regex(root : PurePosixPath, pattern: str) -> List[str]:

    root = Path(root)

    files = root.glob("*")

    list_of_files = []
    for f in files:
        list_of_files.append(re.search(pattern, str(f)).group())

    return list_of_files

def parse_rna(input_dir : str) -> pd.DataFrame:


    rna_dir = Path(input_dir) / 'rna'

    pattern = r"[^/]+-[^-]+-\d+-\d+-\d+-[^/]+$"

    list_of_folders = capture_folders_by_regex(rna_dir, pattern)

    for folder in list_of_folders:

        sv_file_pattern = r".*.fusion_candidates.final$"
        expression_file_pattern = r".*.quant.genes.sf$"
        hard_filtered_pattern = r".*_somatic.hard-filtered.vcf.gz$"
        all_rna_files = rna_dir.glob(f"{folder}/{folder}.RNASeq_somatic/*")


        for f in all_rna_files:
            if re.search(sv_file_pattern, str(f)):
                DATAFRAME['filetype'].append('sv')
            elif re.search(expression_file_pattern,str(f)):
                DATAFRAME['filetype'].append('expression')
            elif re.search(hard_filtered_pattern,str(f)):
                DATAFRAME['filetype'].append('hard_filtered')
            else:
                break

            DATAFRAME['filepath'].append(str(f))

            group_pattern = r"^((?:[^-]*-){2}[^-]*)"
            
            group = re.search(group_pattern, str(folder)).group(1)

            DATAFRAME['group_id'].append(group)

            subject_id_pattern = r"(?:[^-]*-){3}([^-]*)"

            subject_id = re.search(subject_id_pattern, str(folder)).group(1)

            DATAFRAME['subject_id'].append(f"{group}-{subject_id}")

            sequencing_type_pattern = r"-(?:\d+)?([A-Za-z]+)$"

            sequencing_type = re.search(sequencing_type_pattern, str(folder)).group(1)

            if sequencing_type == "RT":
                sequencing_type = "DT"

            run_number_pattern = r"-(\d+)([A-Za-z]+)$"

            run_number = re.search(run_number_pattern, str(folder)).group(1)

            DATAFRAME['sample_id'].append(f"{group}-{subject_id}.{run_number}{sequencing_type}")

            DATAFRAME['sample_type'].append("somatic")

            DATAFRAME['sequence_data'].append("rna")

            DATAFRAME['info'].append(pd.NA)


    return True

def parse_dna(input_dir : str) -> pd.DataFrame:

    dna_dir = Path(input_dir) / 'dna'


    # Capture group ids

    pattern = r"[^/]+-[^-]+-\d+-\d+-\d+-[^/]+$"

    list_of_folders = capture_folders_by_regex(dna_dir, pattern)

    for folder in list_of_folders:

        cna_file_pattern = r".*.cnv.vcf.gz$"
        tumor_hard_filtered_pattern = r".*tumor_normal.hard-filtered.vcf.gz$"
        germinal_hard_filtered_pattern = r".*WGS_germinal.hard-filtered.vcf.gz$"
        all_dna_files = dna_dir.glob(f"{folder}/{folder}.WGS_*/*")


        for f in all_dna_files:
            if re.search(cna_file_pattern, str(f)):
                DATAFRAME['filetype'].append('cnv')
            if re.search(tumor_hard_filtered_pattern,str(f)):
                DATAFRAME['filetype'].append('hard_filtered')
            if re.search(germinal_hard_filtered_pattern,str(f)):
                DATAFRAME['filetype'].append('hard_filtered')

            DATAFRAME['filepath'].append(str(f))

            group_pattern = r"((?:[^-]*-){2}[^-]*)"

            group = re.search(group_pattern, str(folder)).group(1)

            DATAFRAME['group_id'].append(group)
            subject_id_pattern = r"(?:[^-]*-){3}([^-]*)"

            subject_id = re.search(subject_id_pattern, str(folder)).group(1)

            DATAFRAME['subject_id'].append(f"{group}-{subject_id}")
       
            sequencing_type_pattern = r"-(?:\d+)?([A-Za-z]+)$"

            sequencing_type = re.search(sequencing_type_pattern, str(folder)).group(1)
            run_number_pattern = r"-(\d+)([A-Za-z]+)$"

            run_number = re.search(run_number_pattern, str(folder)).group(1)

            DATAFRAME['sample_id'].append(f"{group}-{subject_id}.{run_number}{sequencing_type}")

            if sequencing_type == "DN":
                DATAFRAME['sample_type'].append("germinal")
            elif sequencing_type == "DT":
                DATAFRAME['sample_type'].append("somatic")

            DATAFRAME['sequence_data'].append("dna")

            DATAFRAME['info'].append(pd.NA)


    return True
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--input_dir", required = True, help="Directory containing rna/ and dna/ folders from DRAGEN results.")
    parser.add_argument("--output_dir", required = True, help="Directory that will contain resulting samplesheet file.")

    args = parser.parse_args()

    parse_rna(args.input_dir)
    
    parse_dna(args.input_dir)

    print(DATAFRAME)

    print(pd.DataFrame(DATAFRAME).head())
    pd.DataFrame(DATAFRAME).to_csv("samplesheet.csv", index = False)
