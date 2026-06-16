process ANONYMIZE_CLINICAL {
    label 'python'
    stageInMode 'copy'
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path(clinical_file)
    path(linking_file)
    path(patient_map)

    output:
    path("${clinical_file}")

    script:
    """
    anonymize_clinical.py ${clinical_file} ${linking_file} ${patient_map}
    """
}
