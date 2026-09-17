process MAFSMITH {
    tag "$meta.sample"
    label 'process_medium_memory'
    container "${params.container_mafsmith}"

    input:
        tuple val(meta), path(vcf) // compressed (.vcf.gz) or uncompressed (.vcf)
        path mafsmith_data         // mafsmith home: pre-fetched VEP data + fastvep

    output:
        tuple val(meta), path("${meta.sample}.maf"), emit: maf

    when:
        task.ext.when == null || task.ext.when

    script:
    """
    # mafsmith looks under \$HOME/.mafsmith for its pre-fetched data and the
    # fastvep executable, so point HOME at the task directory and link the bundle in.
    export HOME=./
    ln -s ${mafsmith_data} .mafsmith

    # Handle compressed VCF files
    if [[ $vcf == *.gz ]]; then
        tmp=\$(mktemp --suffix=.vcf)
        rm -f "\$tmp"
        gunzip -c "$vcf" > "\$tmp"
        INPUT_VCF="\$tmp"
    else
        INPUT_VCF="$vcf"
    fi

    cat \$INPUT_VCF | grep "#" > tmp.${meta.sample}.somatic.vcf
    cat \$INPUT_VCF | grep PASS >> tmp.${meta.sample}.somatic.vcf

    TMP_NORMAL_ID=\$(grep '^#CHROM' \$INPUT_VCF | awk '{print \$10}')
    TMP_TUMOR_ID=\$(grep '^#CHROM' \$INPUT_VCF | awk '{print \$11}')

    ID_ARGS=""
    [ -n "\$TMP_TUMOR_ID" ] && ID_ARGS="\$ID_ARGS --tumor-id \$TMP_TUMOR_ID"
    [ -n "\$TMP_NORMAL_ID" ] && ID_ARGS="\$ID_ARGS --normal-id \$TMP_NORMAL_ID"

    mafsmith vcf2maf \\
        \$ID_ARGS \\
        --input-vcf tmp.${meta.sample}.somatic.vcf \\
        --output-maf tmp.${meta.sample}.maf

    # Keep the #version line + column header, then rewrite column 16
    # (Tumor_Sample_Barcode) to the cBioPortal sample id for every data row.
    head -2 tmp.${meta.sample}.maf > ${meta.sample}.maf
    tail -n +3 tmp.${meta.sample}.maf | awk -v col16="${meta.sample}" 'BEGIN {FS=OFS="\\t"} {\$16=col16; print}' >> ${meta.sample}.maf
    """

    stub:
    """
    touch ${meta.sample}.maf
    """
}
