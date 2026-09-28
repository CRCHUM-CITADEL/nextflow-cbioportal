process CONVERT_CPSR_TO_MAF {
    publishDir { "${params.outdir}/${maf_meta.group}/${maf_meta.subject}" }, mode: 'copy'

    tag { maf_meta.sample }
    label "process_medium_memory"
    container params.container_r

    input:
        tuple val(join_key), val(maf_meta), path(som_dna_rna_maf), val(ger_meta), path(ger_dna_tsv_gz)

    output:
        tuple val(maf_meta), path("${maf_meta.sample}.somatic_rna_germline.maf")


    script:
    """
    zcat $ger_dna_tsv_gz > tmp.tsv
    head -1 tmp.tsv > tmp.germline.cpsr.tsv

    # Keep ClinVar P/LP/VUS; column found by name since its position and spelling vary across cpsr versions.
    CLINVAR_COL=\$(head -1 tmp.tsv | tr '\\t' '\\n' | grep -n -x 'CLINVAR_CLASSIFICATION' | cut -d: -f1)
    if [ -z "\$CLINVAR_COL" ]; then
        echo "ERROR: no CLINVAR_CLASSIFICATION column in the header of $ger_dna_tsv_gz" >&2
        exit 1
    fi
    awk -F"\\t" -v col="\$CLINVAR_COL" 'NR>1 && (\$col=="Pathogenic" || \$col=="Likely_Pathogenic" || \$col=="Likely Pathogenic" || \$col=="VUS")' tmp.tsv >> tmp.germline.cpsr.tsv

    rm tmp.tsv # to reduce size of work dir

    gen_convert_cpsr_to_maf.R \
       tmp.germline.cpsr.tsv \
       $som_dna_rna_maf \
       tmp.${maf_meta.sample}.somatic_rna_germline.maf

    head -n2 tmp.${maf_meta.sample}.somatic_rna_germline.maf > ${maf_meta.sample}.somatic_rna_germline.maf
    awk -F'\t' 'NR>2{if(\$9!="Intron" && \$9!="IGR"){print \$0}}' tmp.${maf_meta.sample}.somatic_rna_germline.maf >> ${maf_meta.sample}.somatic_rna_germline.maf
    """

    stub:
    """
    touch "${maf_meta.sample}.somatic_rna_germline.maf"
    """
}
