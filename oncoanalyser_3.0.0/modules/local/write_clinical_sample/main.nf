process WRITE_CLINICAL_SAMPLE {
    publishDir { "${params.outdir}/${group}" }, mode: 'copy'

    container params.container_r

    tag { group }

    input:
        tuple val(group), path(clinical_merged)
        path genomic_attributes // merged GENOMIC_SAMPLE_ATTRIBUTES (both mode), or [] when absent

    output:
        tuple val(group), path("data_clinical_sample.txt"), emit: ch_clinical_sample

    script:
    // cBioPortal allows one SAMPLE_ATTRIBUTES file, so TMB/HRD are joined in here by SAMPLE_ID.
    def join_cmd = genomic_attributes
        ? "mv data_clinical_sample.txt clinical_only.txt && join_genomic_sample_attributes.R --clinical clinical_only.txt --attributes ${genomic_attributes} --output data_clinical_sample.txt"
        : ''
    """
    write_clinical_table.R --input ${clinical_merged} --mode sample -o data_clinical_sample.txt

    ${join_cmd}
    """

    stub:
    """
    touch data_clinical_sample.txt
    """
}
