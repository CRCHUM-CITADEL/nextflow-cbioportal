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

    # Keep P/LP/VUS on CPSR's final call, not CLINVAR_CLASSIFICATION (empty for variants absent from ClinVar).
    # cpsr >= 2.3 names it CLASSIFICATION, earlier FINAL_CLASSIFICATION.
    # A header-only file is PCGR's no-variants placeholder: nothing to filter.
    if [ \$(wc -l < tmp.tsv) -gt 1 ]; then
        CLASS_COL=""
        for name in CLASSIFICATION FINAL_CLASSIFICATION; do
            CLASS_COL=\$(head -1 tmp.tsv | awk -F"\\t" -v name="\$name" '{for (i = 1; i <= NF; i++) if (\$i == name) { print i; exit }}')
            if [ -n "\$CLASS_COL" ]; then break; fi
        done
        if [ -z "\$CLASS_COL" ]; then
            echo "ERROR: no CLASSIFICATION or FINAL_CLASSIFICATION column in the header of $ger_dna_tsv_gz" >&2
            exit 1
        fi
        awk -F"\\t" -v col="\$CLASS_COL" 'NR>1 && (\$col=="Pathogenic" || \$col=="Likely_Pathogenic" || \$col=="Likely Pathogenic" || \$col=="VUS")' tmp.tsv >> tmp.germline.cpsr.tsv
    fi

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
