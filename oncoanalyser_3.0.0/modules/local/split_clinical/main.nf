process SPLIT_CLINICAL {
    tag "$meta.sample"
    label 'process_low'

    container params.container_r

    publishDir { "${params.outdir}/${meta.group}/${meta.subject}" }, mode: 'copy'

    input:
        tuple val(meta), path(clinical_sample), path(clinical_patient)

    output:
        tuple val(meta), path("${meta.sample}.data_clinical_sample.txt"), path("${meta.sample}.data_clinical_patient.txt")

    script:
    """
    head -5 ${clinical_sample} > ${meta.sample}.data_clinical_sample.txt
    awk -F'\t' -v id="${meta.sample}" '\$2 == id' ${clinical_sample} >> ${meta.sample}.data_clinical_sample.txt

    head -5 ${clinical_patient} > ${meta.sample}.data_clinical_patient.txt
    awk -F'\t' -v id="${meta.subject}" '\$1 == id' ${clinical_patient} >> ${meta.sample}.data_clinical_patient.txt
    """

    stub:
    """
    touch ${meta.sample}.data_clinical_sample.txt
    touch ${meta.sample}.data_clinical_patient.txt
    """
}
