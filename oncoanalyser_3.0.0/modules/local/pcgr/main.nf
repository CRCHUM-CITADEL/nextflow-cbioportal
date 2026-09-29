process PCGR {
    tag { meta.sample }
    label 'process_medium_memory'

    container params.container_pcgr

    publishDir { "${params.outdir}/${meta.group}/${meta.subject}" }, mode: 'copy', pattern: "*.cpsr_debug"

    input:
        tuple val(meta), path(ger_dna_vcf), path(ger_dna_vcf_tbi)
        path vep_data
        path ref_data

    output:
        tuple val(meta), path("${meta.sample}.cpsr.grch38.classification.tsv.gz"), emit: tsv
        tuple val(meta), path("${meta.sample}.cpsr_debug"),                         emit: debug

    script:
    def cpsr_id = "${meta.subject}-N"
    def debug   = "${meta.sample}.cpsr_debug"
    def args    = task.ext.args ?: ''   // e.g. '--debug' keeps CPSR's intermediate validation VCFs
    """
    mkdir -p ${debug}

    # CPSR annotates every variant itself and aborts when an input INFO tag clashes with one of its
    # own (PAVE writes IMPACT), so it gets the calls with all INFO removed; FILTER and genotypes stay.
    bcftools annotate -x INFO -Oz -o cpsr_input.vcf.gz ${ger_dna_vcf}
    bcftools index -t cpsr_input.vcf.gz

    # What CPSR is given: sample columns, genotypes (CPSR drops 0/0), FILTER values, contig names.
    summarise_vcf() {
        echo "records: \$(bcftools view -H "\$1" | wc -l)"
        echo "samples: \$(bcftools query -l "\$1" | tr '\\n' ' ')"
        echo "-- GT per sample"; bcftools query -f '[%SAMPLE\\t%GT\\n]' "\$1" | sort | uniq -c
        echo "-- FILTER"; bcftools query -f '%FILTER\\n' "\$1" | sort | uniq -c
        echo "-- contigs"; bcftools query -f '%CHROM\\n' "\$1" | uniq | awk 'NR <= 30'
    }
    { echo "## CPSR input (${ger_dna_vcf}, INFO removed)"; summarise_vcf cpsr_input.vcf.gz; } > ${debug}/input_summary.txt 2>&1 || true

    cpsr \\
    --input_vcf cpsr_input.vcf.gz \\
    --vep_dir ${vep_data} \\
    --refdata_dir $ref_data \\
    --output_dir . \\
    --genome_assembly grch38 \\
    --panel_id 0 \\
    --vep_buffer_size 1000 \\
    --no_html \\
    --sample_id ${cpsr_id} \\
    ${args} 2>&1 | tee ${debug}/cpsr.log

    # CPSR exits 0 even when it fails (e.g. input validation): an error must stop the task, never
    # pass for "no variants" and become the header-only placeholder.
    if grep -q -- '- ERROR -' ${debug}/cpsr.log; then
        echo "CPSR failed for ${meta.sample}:" >&2
        grep -- '- ERROR -' ${debug}/cpsr.log >&2
        exit 1
    fi

    # CPSR's other outputs (annotated PASS VCF/TSV, config, --debug intermediates) for inspection.
    for f in ${cpsr_id}.*; do
        [ -e "\$f" ] && [ "\$f" != "${cpsr_id}.cpsr.grch38.classification.tsv.gz" ] && cp -r "\$f" ${debug}/ || true
    done
    for f in ${cpsr_id}*.pass.vcf.gz; do
        [ -f "\$f" ] || continue
        { echo; echo "## after CPSR validation, panel filter and annotation: \$f"; summarise_vcf "\$f"; } >> ${debug}/input_summary.txt 2>&1 || true
    done

    # CPSR names its output after --sample_id (the normal); downstream expects the tumour sample name.
    if [ -f ${cpsr_id}.cpsr.grch38.classification.tsv.gz ] && [ "${cpsr_id}" != "${meta.sample}" ]; then
        mv ${cpsr_id}.cpsr.grch38.classification.tsv.gz ${meta.sample}.cpsr.grch38.classification.tsv.gz
    fi

    # CPSR writes no file when it finds nothing: emit a header-only placeholder.
    if [ ! -f ${meta.sample}.cpsr.grch38.classification.tsv.gz ]; then
        echo "CPSR wrote no classification file: header-only placeholder used" >> ${debug}/cpsr.log
        echo -e "SAMPLE_ID\tGENOMIC_CHANGE\tGENOME_VERSION\tVCF_SAMPLE_ID\tVARIANT_CLASS\tSYMBOL\tGENE_BIOTYPE\tCODING_STATUS\tEXONIC_STATUS\tCONSEQUENCE\tPROTEIN_CHANGE\tHGVSp\tHGVSc\tCDNA_CHANGE\tTRANSCRIPT_START\tPFAM_DOMAIN\tPFAM_DOMAIN_NAME\tCDS_CHANGE\tEFFECT_PREDICTIONS\tMUTATION_HOTSPOT\tRMSK_HIT\tCALL_CONFIDENCE\tDP_TUMOR\tAF_TUMOR\tDP_CONTROL\tAF_CONTROL\tCONTROL_SAMPLE\tAF_GNOMAD_COMBINED\tAF_GNOMAD_AFR\tAF_GNOMAD_AMR\tAF_GNOMAD_EAS\tAF_GNOMAD_SAS\tAF_GNOMAD_NFE\tAF_GNOMAD_FIN\tAF_GNOMAD_OTH\tDBSNP_RSID\tCLINVAR_CLASSIFICATION\tCLINVAR_MSID\tCLINVAR_VARIANT_ORIGIN\tCLINVAR_CONFLICTED\tCLINVAR_PHENOTYPE\tCPSR_CLASSIFICATION\tCPSR_PATHOGENICITY_SCORE\tCPSR_CLASSIFICATION_CODE\tCPSR_CLASSIFICATION_DOC\tPANEL_OF_NORMALS\tVIRTUAL_PANEL_ID\tPREDICTED_EFFECT" | gzip > ${meta.sample}.cpsr.grch38.classification.tsv.gz
    fi
    """

    stub:
    """
    mkdir -p "${meta.sample}.cpsr_debug"
    echo | gzip > "${meta.sample}.cpsr.grch38.classification.tsv.gz"
    """
}
