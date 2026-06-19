process STUB_MAF {
    tag "${meta.sample_id}"

    input:
    tuple val(meta), path(sample_folder)

    output:
    tuple val(meta), path("${meta.sample_id}.maf")

    script:
    """
    printf 'Hugo_Symbol\\tEntrez_Gene_Id\\tCenter\\tNCBI_Build\\tChromosome\\tStart_Position\\tEnd_Position\\tStrand\\tVariant_Classification\\tVariant_Type\\tReference_Allele\\tTumor_Seq_Allele1\\tTumor_Seq_Allele2\\tdbSNP_RS\\tdbSNP_Val_Status\\tTumor_Sample_Barcode\\n' > "${meta.sample_id}.maf"
    """
}
