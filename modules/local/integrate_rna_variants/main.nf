process INTEGRATE_RNA_VARIANTS {
    tag { subject_id }

    container params.container_r

    input:
        tuple val(subject_id), val(rna_meta), path(som_rna_vcf), val(dna_meta), path(som_dna_maf)

    output:
        tuple val(rna_meta), path("${subject_id}.somatic_rna.maf")

    script:
    """
    zcat $som_rna_vcf | grep "#" > tmp.${subject_id}.vcf
    zcat $som_rna_vcf | grep PASS >> tmp.${subject_id}.vcf
    gen_integrate_rna_variants.R \
        -d $som_dna_maf \
        -r tmp.${subject_id}.vcf \
        -o ${subject_id}.somatic_rna.maf \
        --min_depth=3 --min_vaf=0.05
    """
}
