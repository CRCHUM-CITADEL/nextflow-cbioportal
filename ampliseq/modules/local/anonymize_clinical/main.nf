process ANONYMIZE_CLINICAL {
    label 'python'
    stageInMode 'copy'
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path(clinical_file)
    path(orig_linking)
    path(anon_linking)

    output:
    path("${clinical_file}")

    script:
    """
    anonymize_clinical.py ${clinical_file} ${orig_linking} ${anon_linking}
    """
}
