process GENERATE_CLINICAL_TEMPLATE {
    publishDir "${params.outdir}/${group}", mode: 'copy'

    input:
    val group
    val sample_lines   // newline-separated "subject\tsample" lines
    val oncotree_code

    output:
    tuple val(group), path("data_clinical_sample.txt"), path("data_clinical_patient.txt")

    script:
    def cancer_type = oncotree_code ? oncotree_code.toLowerCase() : "NA"
    """
    # data_clinical_sample.txt
    printf '#Patient Identifier\\tSample Identifier\\tSample Type\\tCancer Type\\n' > data_clinical_sample.txt
    printf '#Identifier to uniquely specify a patient.\\tA unique sample identifier.\\tThe type of sample (i.e., normal, primary, met, recurrence).\\tCancer Type.\\n' >> data_clinical_sample.txt
    printf '#STRING\\tSTRING\\tSTRING\\tSTRING\\n' >> data_clinical_sample.txt
    printf '#1\\t1\\t1\\t1\\n' >> data_clinical_sample.txt
    printf 'PATIENT_ID\\tSAMPLE_ID\\tSAMPLE_TYPE\\tCANCER_TYPE\\n' >> data_clinical_sample.txt

    echo '${sample_lines}' | while IFS=\$'\\t' read -r patient sample; do
        [ -z "\$patient" ] && continue
        printf '%s\\t%s\\tPrimary\\t${cancer_type}\\n' "\$patient" "\$sample" >> data_clinical_sample.txt
    done

    # data_clinical_patient.txt
    printf '#Patient Identifier\\n' > data_clinical_patient.txt
    printf '#Identifier to uniquely specify a patient.\\n' >> data_clinical_patient.txt
    printf '#STRING\\n' >> data_clinical_patient.txt
    printf '#1\\n' >> data_clinical_patient.txt
    printf 'PATIENT_ID\\n' >> data_clinical_patient.txt

    echo '${sample_lines}' | while IFS=\$'\\t' read -r patient sample; do
        [ -z "\$patient" ] && continue
        printf '%s\\n' "\$patient" >> data_clinical_patient.txt
    done
    """
}
