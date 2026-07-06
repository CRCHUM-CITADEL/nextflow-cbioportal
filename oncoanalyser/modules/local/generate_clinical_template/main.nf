process GENERATE_CLINICAL_TEMPLATE {
    publishDir "${params.outdir}/${group}", mode: 'copy'

    input:
    val group
    val sample_lines   // newline-separated "subject\tsample" lines
    val icd_code

    output:
    tuple val(group), path("data_clinical_sample.txt")

    script:
    def cancer_type_code = icd_code ?: "NA"
    """
    printf '#Patient Identifier\\tSample Identifier\\tSample Type\\tCancer Type\\tCancer Type Code\\n' > data_clinical_sample.txt
    printf '#Identifier to uniquely specify a patient.\\tA unique sample identifier.\\tThe type of sample (i.e., normal, primary, met, recurrence).\\tCancer Type.\\tCancer Type Code.\\n' >> data_clinical_sample.txt
    printf '#STRING\\tSTRING\\tSTRING\\tSTRING\\tSTRING\\n' >> data_clinical_sample.txt
    printf '#1\\t1\\t1\\t1\\t1\\n' >> data_clinical_sample.txt
    printf 'PATIENT_ID\\tSAMPLE_ID\\tSAMPLE_TYPE\\tCANCER_TYPE\\tCANCER_TYPE_CODE\\n' >> data_clinical_sample.txt

    echo '${sample_lines}' | while IFS=\$'\\t' read -r patient sample; do
        [ -z "\$patient" ] && continue
        printf '%s\\t%s\\tPrimary\\tNA\\t${cancer_type_code}\\n' "\$patient" "\$sample" >> data_clinical_sample.txt
    done
    """
}
