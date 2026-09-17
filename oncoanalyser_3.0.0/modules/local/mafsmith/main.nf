process MAFSMITH {
    tag "$meta.sample"
    label 'process_medium_memory'
    container "${params.container_mafsmith}"

    input:
        tuple val(meta), path(vcf) // Now accepts both compressed (.vcf.gz) and uncompressed (.vcf) files
        path mafsmith_data         // Required for VEP running. 

    output:
        tuple val(meta), path("${meta.sample}.maf"), emit: maf
        // path "versions.yml"           , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
    """
    # mafsmith will look here for its pre-fetched data and fastvep executable
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

    head -2 tmp.${meta.sample}.maf > ${meta.sample}.maf
    tail -n +3 tmp.${meta.sample}.maf | awk -v col16="${meta.sample}" 'BEGIN {FS=OFS="\\t"} {\$16=col16; print}' >> ${meta.sample}.maf
    """
}
