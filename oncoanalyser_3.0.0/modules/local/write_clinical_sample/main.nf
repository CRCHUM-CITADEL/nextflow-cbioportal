process WRITE_CLINICAL_SAMPLE {
    publishDir { "${params.outdir}/${group}" }, mode: 'copy'

    container params.container_r

    tag { group }

    input:
        tuple val(group), path(clinical_merged)

    output:
        tuple val(group), path("data_clinical_sample.txt"), emit: ch_clinical_sample

    script:
    """
    write_clinical_table.R --input ${clinical_merged} --mode sample -o data_clinical_sample.txt
    """

    stub:
    """
    touch data_clinical_sample.txt
    """
}
