process PURPLE_CNV_TO_CBIOPORTAL {
    tag "$meta.sample"
    label 'process_low'

    container params.container_r

    publishDir { "${params.outdir}/${meta.group}/${meta.subject}" }, mode: 'copy'

    input:
        tuple val(meta), path(purple_cnv_somatic), path(purple_cnv_gene)
        path ensembl_annotations

    output:
        tuple val(meta), path("${meta.sample}_data_cna_hg38.seg"), emit: seg
        tuple val(meta), path("${meta.sample}_data_cna_long.txt"), emit: long

    when:
        task.ext.when == null || task.ext.when

    script:
    """
    #take the first 19 columns as there are 2 buggy columns
    cut -f1-19 ${purple_cnv_gene} > ${meta.sample}.cnv.gene.cleaned.tsv

    Rscript ${projectDir}/bin/gen_purple_cnv_to_cbioportal.R \\
        --purple_cnv_somatic ${purple_cnv_somatic} \\
        --purple_cnv_gene    ${meta.sample}.cnv.gene.cleaned.tsv \\
        --sample_id          ${meta.sample} \\
        --ensembl_annotations ${ensembl_annotations} \\
        --output_seg         ${meta.sample}_data_cna_hg38.seg \\
        --output_long        ${meta.sample}_data_cna_long.txt
    """
}
