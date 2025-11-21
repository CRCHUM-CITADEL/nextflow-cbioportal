// modified to allow vcf.gz and vcf (10/10/2025)
process VCF2MAF {
    tag "$meta.sample"
    label 'process_medium_memory'
    container "${params.container_vcf2maf}"

    input:
        tuple val(meta), path(vcf) // Now accepts both compressed (.vcf.gz) and uncompressed (.vcf) files
        path fasta                 // Required
        path vep_data // Required for VEP running. A default of /.vep is supplied.

    output:
        tuple val(meta), path("${meta.sample}.maf"), emit: maf

    script:
    def VEP_CMD       = "--vep-path ${vep_data} --vep-data ${vep_data}/cache ${params.vep_params}"
    def VERSION       = '1.6.22' 
    """

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

    ## TODO: is DN always first?
    TMP_NORMAL_ID=\$(grep "^#CHROM" \$INPUT_VCF | awk '{print \$10}')
    TMP_TUMOR_ID=\$(grep "^#CHROM" \$INPUT_VCF | awk '{print \$11}')
    
    vcf2maf.pl \\
        --tumor-id \$TMP_TUMOR_ID \\
        --normal-id \$TMP_NORMAL_ID \\
        $VEP_CMD \\
        --ref-fasta $fasta \\
        --input-vcf tmp.${meta.sample}.somatic.vcf \\
        --output-maf tmp.${meta.sample}.maf

    head -2 tmp.${meta.sample}.maf > ${meta.sample}.maf
    tail -n +3 tmp.${meta.sample}.maf | awk -v col16="${meta.sample}" -v col17="${meta.germinal_sample}" 'BEGIN {FS=OFS="\\t"} NR==0 {print; next} {\$16=col16; \$17=col17; print}' >> ${meta.sample}.maf
    """
}
