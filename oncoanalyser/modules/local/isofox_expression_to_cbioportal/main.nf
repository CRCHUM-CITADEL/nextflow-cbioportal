process ISOFOX_EXPRESSION_TO_CBIOPORTAL {
    tag "$meta.sample"
    label 'process_low'

    container params.container_r

    publishDir { "${params.outdir}/${meta.group}/${meta.subject}" }, mode: 'copy'

    input:
        tuple val(meta), path(exp_tsv)
        path ensembl_annotations

    output:
        tuple val(meta), path("${meta.sample}.tpm.tsv"), emit: tpm

    when:
        task.ext.when == null || task.ext.when

    script:
    """
    Rscript ${projectDir}/bin/gen_isofox_expression_to_cbioportal.R \\
        --input    ${exp_tsv} \\
        --gene_map ${ensembl_annotations} \\
        --sample_id ${meta.sample} \\
        --output   ${meta.sample}.tpm.tsv
    """
}
