process GENE_CNV_FOLD_CHANGES_TO_CBIOPORTAL {
    publishDir "${params.outdir}/${meta.group}/${meta.subject}", mode: 'copy'

    // use meta.sample_id for logging
    tag { meta.subject }

    container params.container_r

    input:
    tuple val(meta), path(somatic_cnv_vcf), path(fold_changes_per_gene_cnv)

    output:
    tuple val(meta), path("${meta.subject}_data_cna_hg38.seg"), emit : seg
    tuple val(meta), path("${meta.subject}_data_cna_long.txt"), emit : long

    script:
    """
    gen_cbioportal_converter.R \
      --vcf $somatic_cnv_vcf \
      --tsv $fold_changes_per_gene_cnv \
      --sample_id ${meta.subject} \
      --output_dir .
    """
}
