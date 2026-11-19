process MAFSMITH {
    tag "$meta.sample"
    label 'process_medium_memory'
    container "${params.container_mafsmith}"

    input:
        tuple val(meta), path(vcf) // compressed (.vcf.gz) or uncompressed (.vcf)
        path mafsmith_data         // mafsmith home: reference data (+ optionally fastvep)
        path ensembl_annotations   // BioMart TSV: fills Entrez_Gene_Id, which mafsmith leaves 0

    output:
        tuple val(meta), path("${meta.sample}.maf"), emit: maf

    when:
        task.ext.when == null || task.ext.when

    script:
    // Extra CSQ subfields kept as MAF columns; names must match fastVEP's DEFAULT_CSQ_FIELDS (unknown → empty column)
    def retain_ann = params.mafsmith_retain_ann ? "--retain-ann ${params.mafsmith_retain_ann}" : ''
    def args       = task.ext.args ?: ''
    """
    # mafsmith reads only \$HOME/.mafsmith, so HOME is the task dir with the bundle linked in
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

    # Sample columns follow FORMAT on the #CHROM line (normal, then tumor), found by that header
    SAMPLE_IDS=\$(awk 'BEGIN {FS=OFS="\\t"} /^#CHROM/ {
        for (i = 1; i <= NF; i++) if (\$i == "FORMAT") { print \$(i+1), \$(i+2); exit }
    }' \$INPUT_VCF)
    TMP_NORMAL_ID=\$(printf '%s' "\$SAMPLE_IDS" | cut -f1)
    TMP_TUMOR_ID=\$(printf '%s' "\$SAMPLE_IDS" | cut -f2)

    # Without --vcf-tumor-id/--vcf-normal-id mafsmith takes the FIRST sample (here the normal) as tumour
    ID_ARGS=""
    [ -n "\$TMP_TUMOR_ID" ] && ID_ARGS="\$ID_ARGS --tumor-id \$TMP_TUMOR_ID --vcf-tumor-id \$TMP_TUMOR_ID"
    [ -n "\$TMP_NORMAL_ID" ] && ID_ARGS="\$ID_ARGS --normal-id \$TMP_NORMAL_ID --vcf-normal-id \$TMP_NORMAL_ID"

    mafsmith vcf2maf \\
        \$ID_ARGS \\
        ${retain_ann} \\
        ${args} \\
        --input-vcf tmp.${meta.sample}.somatic.vcf \\
        --output-maf tmp.${meta.sample}.maf

    # Join back the INFO fields --retain-ann drops (gnomAD, TIER, PURPLE_*, ...); always written so the merge stays rectangular
    annotate_maf_with_vcf_info.py \\
        --vcf tmp.${meta.sample}.somatic.vcf \\
        --maf tmp.${meta.sample}.maf \\
        --output tmp.${meta.sample}.annotated.maf
    mv tmp.${meta.sample}.annotated.maf tmp.${meta.sample}.maf

    # mafsmith has no Ensembl-to-Entrez table and writes Entrez_Gene_Id = 0 on every row.
    fill_maf_entrez.py \\
        --maf tmp.${meta.sample}.maf \\
        --annotations ${ensembl_annotations} \\
        --output tmp.${meta.sample}.entrez.maf
    mv tmp.${meta.sample}.entrez.maf tmp.${meta.sample}.maf

    # Tumor_Sample_Barcode = cBioPortal sample id, Mutation_Status = Somatic (mafsmith leaves it empty); columns by name on line 2
    col_of() { sed -n '2p' tmp.${meta.sample}.maf | awk -F"\\t" -v name="\$1" '{for (i = 1; i <= NF; i++) if (\$i == name) { print i; exit }}'; }
    SAMPLE_COL=\$(col_of Tumor_Sample_Barcode)
    STATUS_COL=\$(col_of Mutation_Status)
    if [ -z "\$SAMPLE_COL" ] || [ -z "\$STATUS_COL" ]; then
        echo "ERROR: no Tumor_Sample_Barcode or Mutation_Status column in the header of tmp.${meta.sample}.maf" >&2
        exit 1
    fi

    head -2 tmp.${meta.sample}.maf > ${meta.sample}.maf
    tail -n +3 tmp.${meta.sample}.maf \\
        | awk -v col="\$SAMPLE_COL" -v scol="\$STATUS_COL" -v sample="${meta.sample}" 'BEGIN {FS=OFS="\\t"} {\$col=sample; \$scol="Somatic"; print}' \\
        >> ${meta.sample}.maf
    """

    stub:
    """
    touch ${meta.sample}.maf
    """
}
