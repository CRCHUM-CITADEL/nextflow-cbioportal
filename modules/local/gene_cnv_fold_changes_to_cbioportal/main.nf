process GENE_CNV_FOLD_CHANGES_TO_CBIOPORTAL {
    publishDir "${params.outdir}/${meta.group}/${meta.sample}", mode: 'copy'

    // use meta.sample_id for logging
    tag { meta.sample }

    container params.container_r

    input:
    tuple val(meta), path(somatic_cnv_vcf), path(fold_changes_per_gene_cnv)

    output:
    tuple val(meta), path("${meta.sample}_data_cna_hg38.seg"), emit : seg
    tuple val(meta), path("${meta.sample}_data_cna_long.txt"), emit : long

    script:
    """
    gen_cbioportal_converter.R \
      --vcf $somatic_cnv_vcf \
      --tsv $fold_changes_per_gene_cnv \
      --sample_id ${meta.sample} \
      --output_dir .
    """
}
