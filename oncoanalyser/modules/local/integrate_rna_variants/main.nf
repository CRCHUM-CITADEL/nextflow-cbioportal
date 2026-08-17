process INTEGRATE_RNA_VARIANTS {
    publishDir { "${params.outdir}/${dna_meta.group}/${dna_meta.subject}" }, mode: 'copy'
    tag { dna_meta.sample }


    label "process_medium_memory"
    container params.container_r

    input:
        tuple val(join_key), val(rna_meta), path(som_rna_vcf), val(dna_meta), path(som_dna_maf)

    output:
        tuple val(dna_meta), path("${dna_meta.sample}.somatic_rna.maf")

    script:
    """
    zcat $som_rna_vcf | grep "#" > tmp.${dna_meta.sample}.vcf
    zcat $som_rna_vcf | grep PASS >> tmp.${dna_meta.sample}.vcf
    gen_integrate_rna_variants.R \
        -d $som_dna_maf \
        -r tmp.${dna_meta.sample}.vcf \
        -o ${dna_meta.sample}.somatic_rna.maf \
        --min_depth=3 --min_vaf=0.05
    """
}
