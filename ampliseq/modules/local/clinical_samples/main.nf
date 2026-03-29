process CLINICAL_SAMPLES {
    label 'python'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(sample_file)

    output:
    path("data_clinical_sample.txt")

    script:
    """
    clinical_sample_format.py ${sample_file}
    """
}
