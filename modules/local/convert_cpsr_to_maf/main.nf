process CONVERT_CPSR_TO_MAF {
    publishDir "${params.outdir}/${maf_meta.group}/${maf_meta.sample}", mode: 'copy'

    tag { subject_id }

    container params.container_r

    input:
        tuple val(subject_id), val(maf_meta), path(som_dna_rna_maf), val(ger_meta), path(ger_dna_tsv_gz)

    output:
        tuple val(maf_meta), path("${subject_id}.somatic_rna_germline.maf")


    script:
    """
    zcat $ger_dna_tsv_gz > tmp.tsv
    head -1 tmp.tsv > tmp.germline.cpsr.tsv
    awk -F"\\t" '\$52=="Pathogenic" || \$52=="Likely_Pathogenic" || \$52=="VUS"' tmp.tsv >> tmp.germline.cpsr.tsv

    rm tmp.tsv # to reduce size of work dir

    gen_convert_cpsr_to_maf.R \
        tmp.germline.cpsr.tsv \
        $som_dna_rna_maf \
        tmp.${subject_id}.somatic_rna_germline.maf

    head -n2 tmp.${subject_id}.somatic_rna_germline.maf > ${subject_id}.somatic_rna_germline.maf
    awk -F'\t' 'NR>2{if(\$9!="Intron" && \$9!="IGR"){print \$0}}' tmp.${subject_id}.somatic_rna_germline.maf >> ${subject_id}.somatic_rna_germline.maf
    """

}
