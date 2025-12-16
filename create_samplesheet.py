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
from typing import Any, List
from collections.abc import Mapping
from collections import Counter

DATAFRAME = {"group_id" : [],"subject_id": [],"sample_id" : [],"sample_type" : [],"sequence_data" : [],"info" : [],"filepath" : []}

BASE_FOLDER_PATTERN = r"[^/]+-[^-]+-\d+-\d+-\d+-[^/]+$"
GROUP_PATTERN = r"((?:[^-]*-){2}[^-]*)"
SUBJECT_ID_PATTERN = r"(?:[^-]*-){3}([^-]*)"
SEQUENCING_TYPE_PATTERN = r"-(?:\d+)([A-Za-z]+)"
RUN_NUMBER_PATTERN = r"-(\d+)([A-Za-z]+)"
 
def capture_folders_by_regex(root : PurePosixPath, pattern: str) -> List[Any]:

    root = Path(root)

    # should list the last folders before containers the files (so for each sample, 1 DT folder, 1RT folder and another DN folder)
    folders = root.glob("*")

    list_of_folders = []
    for f in folders:
        list_of_folders.append(re.search(pattern, str(f)).group())

    return list_of_folders

def parse_folders(input_dir : str) -> bool:

    # 3 rows per sample, DT (somatic dna), DN (germinal DNA) and RT (somatic RNA)

    print(input_dir)

    rna_dir = Path(input_dir) / 'rna'
    dna_dir = Path(input_dir) / 'dna'

    list_of_rna_folders = capture_folders_by_regex(rna_dir, BASE_FOLDER_PATTERN)
    list_of_dna_folders = capture_folders_by_regex(dna_dir, BASE_FOLDER_PATTERN)

    for folder in list_of_rna_folders:
        parse_rna(folder, input_dir)

    for folder in list_of_dna_folders:
        parse_dna(folder, input_dir)

    return True

def parse_rna(folder: str, root_path : str) -> bool:


    group = re.search(GROUP_PATTERN, folder).group(1)
    DATAFRAME['group_id'].append(group)

    subject_id = re.search(SUBJECT_ID_PATTERN, folder).group(1)
    DATAFRAME['subject_id'].append(f"{group}-{subject_id}")

    sequencing_type = re.search(SEQUENCING_TYPE_PATTERN, folder).group(1)

    # special function that will get the somatic tumor run number for the equivalent subject 
    run_number = get_somatic_tumor_dna_run_number(Path(f"{root_path}/dna/"), f"{group}-{subject_id}")
    sequencing_type = "DT"

    # run_number = re.search(RUN_NUMBER_PATTERN, folder).group(1)

    sample_id = f"{group}-{subject_id}.{run_number}{sequencing_type}"

    DATAFRAME['sample_id'].append(sample_id)
    
    DATAFRAME['filepath'].append(root_path + "/rna/" + folder + "/" + folder + ".RNASeq_somatic")

    DATAFRAME['sample_type'].append("somatic")

    DATAFRAME['sequence_data'].append("rna")

    DATAFRAME['info'].append(pd.NA)

    return True

def parse_dna(folder : str, root_path : str) -> bool:

    group = re.search(GROUP_PATTERN, folder).group(1)
    DATAFRAME['group_id'].append(group)

    subject_id = re.search(SUBJECT_ID_PATTERN, folder).group(1)
    DATAFRAME['subject_id'].append(f"{group}-{subject_id}")

    sequencing_type = re.search(SEQUENCING_TYPE_PATTERN, folder).group(1)
    run_number = re.search(RUN_NUMBER_PATTERN, folder).group(1)
    
    sample_id = f"{group}-{subject_id}.{run_number}{sequencing_type}"
    DATAFRAME['sample_id'].append(sample_id)

    if sequencing_type == "DN":
        DATAFRAME['sample_type'].append("germinal")
        DATAFRAME['filepath'].append(root_path + "/dna/" + folder + "/" + folder + ".WGS_germinal")

    elif sequencing_type == "DT":
        DATAFRAME['sample_type'].append("somatic")
        DATAFRAME['filepath'].append(root_path + "/dna/" + folder + "/" + folder + ".WGS_somatic-tumor_normal")

    DATAFRAME['sequence_data'].append("dna")

    DATAFRAME['info'].append(pd.NA)

    return True

def get_somatic_tumor_dna_run_number(base_path : Path,  subject_id : str) -> str:

    dna_files = base_path.glob(f"{subject_id}*DT")

    # take the first one
    for dna_file in dna_files:
    
        run_number = re.search(RUN_NUMBER_PATTERN, str(dna_file)).group(1)
        
        return run_number

## tests
def are_dicts_equal(d1, d2):
    if d1.keys() != d2.keys():
        return False

    for key in d1:
        v1, v2 = d1[key], d2[key]
        
        if isinstance(v1, Mapping) and isinstance(v2, Mapping):  # Recursively check dictionaries
            if not are_dicts_equal(v1, v2):
                return False
        elif isinstance(v1, list) and isinstance(v2, list):  # Check lists ignoring order
            if Counter(v1) != Counter(v2):
                return False
        else:
            if v1 != v2:
                return False
    return True

def run_tests():

    current_directory = Path(__file__).resolve().parent

    return_value = [
            Path(f"{current_directory}/path/to/rna/MoHQ-CM-3-264-592928-1RT/MoHQ-CM-3-264-592928-1RT.RNASeq_somatic"),
            Path(f"{current_directory}/path/to/rna/MoHQ-CM-3-261-581463-1RT/MoHQ-CM-3-261-581463-1RT.RNASeq_somatic"),
            Path(f"{current_directory}/path/to/dna/MoHQ-CM-3-264-592928-2DT/MoHQ-CM-3-264-592928-2DT.WGS_somatic-tumor_normal"),
            Path(f"{current_directory}/path/to/dna/MoHQ-CM-3-261-581463-1DT/MoHQ-CM-3-261-581463-1DT.WGS_somatic-tumor_normal"),
            Path(f"{current_directory}/path/to/dna/MoHQ-CM-3-264-620280-1DN/MoHQ-CM-3-264-620280-1DN.WGS_germinal"),
            Path(f"{current_directory}/path/to/dna/MoHQ-CM-3-261-620167-1DN/MoHQ-CM-3-261-620167-1DN.WGS_germinal")
        ]
    
    # create tmp files
    for path in return_value:
        path.parent.mkdir(parents=True, exist_ok=True)

    root_path = f"{current_directory}/tmp/path/to"

    print(current_directory)

    parse_folders(root_path)

    correct_df = {
        "group_id" : [
            "MoHQ-CM-3",
            "MoHQ-CM-3",
            "MoHQ-CM-3",
            "MoHQ-CM-3",
            "MoHQ-CM-3",
            "MoHQ-CM-3",
        ],
        "subject_id": [
            "MoHQ-CM-3-264",
            "MoHQ-CM-3-261",
            "MoHQ-CM-3-264",
            "MoHQ-CM-3-261",
            "MoHQ-CM-3-264",
            "MoHQ-CM-3-261"
        ],
        "sample_id" : [
            "MoHQ-CM-3-264.2DT",
            "MoHQ-CM-3-261.1DT",
            "MoHQ-CM-3-264.2DT",
            "MoHQ-CM-3-261.1DT",
            "MoHQ-CM-3-264.1DN",
            "MoHQ-CM-3-261.1DN",
        ],
        "sample_type" : [
            "somatic",
            "somatic",
            "somatic",
            "somatic",
            "germinal",
            "germinal"
        ],
        "sequence_data" : [
            "rna",
            "rna",
            "dna",
            "dna",
            "dna",
            "dna",
        ],
        "info" : [
            pd.NA,
            pd.NA,
            pd.NA,
            pd.NA,
            pd.NA,
            pd.NA,
        ],
        "filepath" : [
            f"{current_directory}/tmp/path/to/rna/MoHQ-CM-3-264-592928-1RT/MoHQ-CM-3-264-592928-1RT.RNASeq_somatic",
            f"{current_directory}/tmp/path/to/rna/MoHQ-CM-3-261-581463-1RT/MoHQ-CM-3-261-581463-1RT.RNASeq_somatic",
            f"{current_directory}/tmp/path/to/dna/MoHQ-CM-3-264-592928-2DT/MoHQ-CM-3-264-592928-2DT.WGS_somatic-tumor_normal",
            f"{current_directory}/tmp/path/to/dna/MoHQ-CM-3-261-581463-1DT/MoHQ-CM-3-261-581463-1DT.WGS_somatic-tumor_normal",
            f"{current_directory}/tmp/path/to/dna/MoHQ-CM-3-264-620280-1DN/MoHQ-CM-3-264-620280-1DN.WGS_germinal",
            f"{current_directory}/tmp/path/to/dna/MoHQ-CM-3-261-620167-1DN/MoHQ-CM-3-261-620167-1DN.WGS_germinal"
        ]
    }

    assert are_dicts_equal(DATAFRAME, correct_df), f"not the correct final dataframe output: {DATAFRAME} vs {correct_df}"
    print("TEST PASSED!")

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--input_dir", required = False, help="Directory containing rna/ and dna/ folders from DRAGEN results.")
    parser.add_argument("--output_dir", required = False, help="Directory that will contain resulting samplesheet file.")
    parser.add_argument("--test", action="store_true", required = False, help="Run tests for script")
    args = parser.parse_args()

    if args.test:
        run_tests()
        quit()

    success = parse_folders(args.input_dir)

    if not success:
        raise Exception("unable to save dataframe")

    print(DATAFRAME)

    print(pd.DataFrame(DATAFRAME).head())
    pd.DataFrame(DATAFRAME).to_csv("samplesheet.csv", index = False)



