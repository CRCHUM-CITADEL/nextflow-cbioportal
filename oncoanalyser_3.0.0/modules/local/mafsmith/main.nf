process MAFSMITH {
    tag "$meta.sample"
    label 'process_medium_memory'
    container "${params.container_mafsmith}"

    input:
        tuple val(meta), path(vcf) // compressed (.vcf.gz) or uncompressed (.vcf)
        path mafsmith_data         // mafsmith home: reference data (+ optionally fastvep)

    output:
        tuple val(meta), path("${meta.sample}.maf"), emit: maf

    when:
        task.ext.when == null || task.ext.when

    script:
    """
    # mafsmith looks under \$HOME/.mafsmith for its reference data, so point HOME at the
    # task directory and link the bundle in.
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

    # Sample columns are the two that follow FORMAT on the #CHROM line (normal, then
    # tumor). Locating them by that header keeps this correct regardless of how many
    # fixed columns the VCF carries.
    SAMPLE_IDS=\$(awk 'BEGIN {FS=OFS="\\t"} /^#CHROM/ {
        for (i = 1; i <= NF; i++) if (\$i == "FORMAT") { print \$(i+1), \$(i+2); exit }
    }' \$INPUT_VCF)
    TMP_NORMAL_ID=\$(printf '%s' "\$SAMPLE_IDS" | cut -f1)
    TMP_TUMOR_ID=\$(printf '%s' "\$SAMPLE_IDS" | cut -f2)

    ID_ARGS=""
    [ -n "\$TMP_TUMOR_ID" ] && ID_ARGS="\$ID_ARGS --tumor-id \$TMP_TUMOR_ID"
    [ -n "\$TMP_NORMAL_ID" ] && ID_ARGS="\$ID_ARGS --normal-id \$TMP_NORMAL_ID"

    mafsmith vcf2maf \\
        \$ID_ARGS \\
        --input-vcf tmp.${meta.sample}.somatic.vcf \\
        --output-maf tmp.${meta.sample}.maf

    # Rewrite Tumor_Sample_Barcode to the cBioPortal sample id. Line 1 of a MAF is
    # #version and line 2 the column header, so the header names are read from line 2 and
    # the column is located by name rather than assumed to sit at a fixed position.
    SAMPLE_COL=\$(sed -n '2p' tmp.${meta.sample}.maf | tr '\\t' '\\n' | grep -n -x 'Tumor_Sample_Barcode' | cut -d: -f1)
    if [ -z "\$SAMPLE_COL" ]; then
        echo "ERROR: no Tumor_Sample_Barcode column in the header of tmp.${meta.sample}.maf" >&2
        exit 1
    fi

    head -2 tmp.${meta.sample}.maf > ${meta.sample}.maf
    tail -n +3 tmp.${meta.sample}.maf \\
        | awk -v col="\$SAMPLE_COL" -v sample="${meta.sample}" 'BEGIN {FS=OFS="\\t"} {\$col=sample; print}' \\
        >> ${meta.sample}.maf
    """

    stub:
    """
    touch ${meta.sample}.maf
    """
}
