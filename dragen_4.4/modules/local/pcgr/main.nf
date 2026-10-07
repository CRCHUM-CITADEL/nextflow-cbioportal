process PCGR {
    tag { meta.sample }
    label 'process_medium_memory'

    container params.container_pcgr

    input:
        tuple val(meta), path(ger_dna_vcf), path(ger_dna_vcf_tbi)
        path vep_data
        path ref_data

    output:
        tuple val(meta), path("${meta.sample}.cpsr.grch38.classification.tsv.gz")

    script:
    """
    # rename needed to match with index file (skip if already matching)
    if [ "${ger_dna_vcf}" != "${meta.sample}.vcf.gz" ]; then
        mv ${ger_dna_vcf} ${meta.sample}.vcf.gz
    fi


    cpsr \
    --input_vcf ${meta.sample}.vcf.gz \
    --vep_dir ${vep_data} \
    --refdata_dir $ref_data \
    --output_dir . \
    --genome_assembly grch38 \
    --panel_id 0 \
    --no_html \
    --vep_buffer_size 1000 \
    --sample_id ${meta.sample}

    # cpsr exits silently when the VCF has no variants in the selected cancer
    # predisposition geneset (common on chr21-only test data). Emit an empty
    # header-only file so the process still satisfies the declared output.
    [ ! -f ${meta.sample}.cpsr.grch38.classification.tsv.gz ] && echo -e "SAMPLE_ID\\tGENOMIC_CHANGE\\tGENOME_VERSION\\tVCF_SAMPLE_ID\\tVARIANT_CLASS\\tSYMBOL\\tGENE_BIOTYPE\\tCODING_STATUS\\tEXONIC_STATUS\\tCONSEQUENCE\\tPROTEIN_CHANGE\\tHGVSp\\tHGVSc\\tCDNA_CHANGE\\tTRANSCRIPT_START\\tPFAM_DOMAIN\\tPFAM_DOMAIN_NAME\\tCDS_CHANGE\\tEFFECT_PREDICTIONS\\tMUTATION_HOTSPOT\\tRMSK_HIT\\tCALL_CONFIDENCE\\tDP_TUMOR\\tAF_TUMOR\\tDP_CONTROL\\tAF_CONTROL\\tCONTROL_SAMPLE\\tAF_GNOMAD_COMBINED\\tAF_GNOMAD_AFR\\tAF_GNOMAD_AMR\\tAF_GNOMAD_EAS\\tAF_GNOMAD_SAS\\tAF_GNOMAD_NFE\\tAF_GNOMAD_FIN\\tAF_GNOMAD_OTH\\tDBSNP_RSID\\tCLINVAR_CLASSIFICATION\\tCLINVAR_MSID\\tCLINVAR_VARIANT_ORIGIN\\tCLINVAR_CONFLICTED\\tCLINVAR_PHENOTYPE\\tCPSR_CLASSIFICATION\\tCPSR_PATHOGENICITY_SCORE\\tCPSR_CLASSIFICATION_CODE\\tCPSR_CLASSIFICATION_DOC\\tPANEL_OF_NORMALS\\tVIRTUAL_PANEL_ID\\tPREDICTED_EFFECT" | gzip > ${meta.sample}.cpsr.grch38.classification.tsv.gz || true
    """
}
