process GET_TPM {
    tag { meta.sample }   // helps logging/tracing per sample

    container params.container_r

    input:
        tuple val(meta), path(somatic_expression_file)      // one sample id + corresponding .quant.genes.sf  file
        path ensembl_annotations                            // one gene annotations file (biomart ensembl)

    output:
        tuple val(meta), path("${meta.sample}.tpm.tsv")

    script:
    """
    awk -F'\\t' 'NR>1{print \$1"\\t"\$9"\\t"\$2}' $ensembl_annotations \
    	> gene_id_to_name.tsv

    gen_get_tpm.R \
        --input $somatic_expression_file \
        --gene_map gene_id_to_name.tsv \
        --sample ${meta.sample} \
        --output ${meta.sample}.tpm.tsv
    """
}
