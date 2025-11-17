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
    cpsr \
    --input_vcf $ger_dna_vcf \
    --vep_dir ${vep_data}/cache \
    --refdata_dir $ref_data \
    --output_dir . \
    --genome_assembly grch38 \
    --panel_id 0 \
    --sample ${meta.sample}
    """
}
