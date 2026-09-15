process WRITE_CLINICAL_PATIENT {
    publishDir { "${params.outdir}/${group}" }, mode: 'copy'

    container params.container_r

    tag { group }

    input:
        tuple val(group), path(clinical_merged)

    output:
        tuple val(group), path("data_clinical_patient.txt"), emit: ch_clinical_patient

    script:
    """
    write_clinical_table.R --input ${clinical_merged} --mode patient -o data_clinical_patient.txt
    """

    stub:
    """
    touch data_clinical_patient.txt
    """
}
