process GENERATE_CLINICAL_TEMPLATE {
    publishDir { "${params.outdir}/${group}" }, mode: 'copy'

    input:
    val group
    val sample_lines   // newline-separated "subject\tsample\tsample_type" lines

    output:
    tuple val(group), path("data_clinical_sample.txt")

    script:
    """
    printf '#Patient Identifier\\tSample Identifier\\tSample Type\\n' > data_clinical_sample.txt
    printf '#Identifier to uniquely specify a patient.\\tA unique sample identifier.\\tThe type of sample (i.e., normal, primary, met, recurrence).\\n' >> data_clinical_sample.txt
    printf '#STRING\\tSTRING\\tSTRING\\n' >> data_clinical_sample.txt
    printf '#1\\t1\\t1\\n' >> data_clinical_sample.txt
    printf 'PATIENT_ID\\tSAMPLE_ID\\tSAMPLE_TYPE\\n' >> data_clinical_sample.txt

    echo '${sample_lines}' | while IFS=\$'\\t' read -r patient sample sample_type; do
        [ -z "\$patient" ] && continue
        printf '%s\\t%s\\t%s\\n' "\$patient" "\$sample" "\${sample_type:-Primary}" >> data_clinical_sample.txt
    done
    """
}
