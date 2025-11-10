// modules/local/gene_cnv_fold_changes/main.nf

process EXTRACT_GENE_CNV_FOLD_CHANGES {
    tag { meta.sample }   // helps logging/tracing per sample

    container params.container_r

    input:
      tuple val(meta), path(somatic_cnv_vcf)      // one sample id + corresponding vcf.gz file
      path gene_annotations                             // one gene annotations file

    output:
      path "${meta.sample}.genes.cnv.tsv"


    script:
    """
    zcat $somatic_cnv_vcf  | grep "#" > ${meta.sample}.somatic.cnv.vcf
    zcat $somatic_cnv_vcf  | grep PASS >> ${meta.sample}.somatic.cnv.vcf

    gen_gene_cnv_fold_changes.R \
      --vcf ${meta.sample}.somatic.cnv.vcf \
      --annotation $gene_annotations \
      --output ${meta.sample}.genes.cnv.tsv
    """
}
