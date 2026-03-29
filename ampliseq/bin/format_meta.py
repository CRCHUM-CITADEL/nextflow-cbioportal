#!/usr/bin/env python3
"""
Generate cBioPortal meta files for all data types.
Usage: python3 format_meta.py <study_id> [out_dir]
"""

import sys
import os

TEMPLATES = {
    "meta_study.txt": """\
type_of_cancer: mixed
cancer_study_identifier: {study_id}
name: {study_id}
short_name: {study_id}
description: cBioPortal study for {study_id}
add_global_case_list: true
reference_genome: hg19
""",
    "meta_mutations.txt": """\
cancer_study_identifier: {study_id}
genetic_alteration_type: MUTATION_EXTENDED
datatype: MAF
stable_id: mutations
show_profile_in_analysis_tab: true
profile_name: Mutations
profile_description: Mutation data for {study_id}
data_filename: data_mutations.txt
gene_panel: gene_panel
""",
    "meta_sv.txt": """\
cancer_study_identifier: {study_id}
genetic_alteration_type: STRUCTURAL_VARIANT
datatype: SV
stable_id: structural_variants
show_profile_in_analysis_tab: true
profile_name: Structural Variants
profile_description: Structural variant data for {study_id}
data_filename: data_sv.txt
gene_panel: ampliseq_focus_fusion
""",
    "meta_cna.txt": """\
cancer_study_identifier: {study_id}
genetic_alteration_type: COPY_NUMBER_ALTERATION
datatype: DISCRETE_LONG
stable_id: cna
show_profile_in_analysis_tab: true
profile_name: Copy-number alterations
profile_description: Discrete copy number data for {study_id}. Values: -2 = homozygous deletion; -1 = hemizygous deletion; 0 = neutral / no change; 1 = gain; 2 = high level amplification.
data_filename: data_cna.txt
gene_panel: ampliseq_focus_cnv
""",
    "meta_clinical_patient.txt": """\
cancer_study_identifier: {study_id}
genetic_alteration_type: CLINICAL
datatype: PATIENT_ATTRIBUTES
data_filename: data_clinical_patient.txt
""",
    "meta_clinical_sample.txt": """\
cancer_study_identifier: {study_id}
genetic_alteration_type: CLINICAL
datatype: SAMPLE_ATTRIBUTES
data_filename: data_clinical_sample.txt
""",
}


def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <study_id> [out_dir]", file=sys.stderr)
        sys.exit(1)

    study_id = sys.argv[1]
    out_dir = sys.argv[2] if len(sys.argv) > 2 else "."

    os.makedirs(out_dir, exist_ok=True)

    for filename, template in TEMPLATES.items():
        path = os.path.join(out_dir, filename)
        with open(path, "w") as f:
            f.write(template.format(study_id=study_id))
        print(f"Written: {path}")


if __name__ == "__main__":
    main()
