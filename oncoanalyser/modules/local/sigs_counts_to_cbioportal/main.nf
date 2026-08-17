process SIGS_COUNTS_TO_CBIOPORTAL {
    tag "$meta.sample"
    label 'process_single'

    container params.container_r

    publishDir { "${params.outdir}/${meta.group}/${meta.subject}" }, mode: 'copy'

    input:
        tuple val(meta), path(snv_counts)

    output:
        tuple val(meta), path("${meta.sample}.data_mutational_signatures_counts_SBS.txt"), emit: sigs_counts

    when:
        task.ext.when == null || task.ext.when

    script:
    """
    Rscript ${projectDir}/bin/gen_sigs_counts_to_cbioportal.R \\
        --input  ${snv_counts} \\
        --sample ${meta.sample} \\
        --output ${meta.sample}.data_mutational_signatures_counts_SBS.txt
    """

    stub:
    """
    touch ${meta.sample}.data_mutational_signatures_counts_SBS.txt
    """
}
