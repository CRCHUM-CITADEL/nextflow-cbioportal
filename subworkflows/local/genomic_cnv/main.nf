include { EXTRACT_GENE_CNV_FOLD_CHANGES } from '../../../modules/local/extract_gene_cnv_fold_changes'
include { GENE_CNV_FOLD_CHANGES_TO_CBIOPORTAL } from '../../../modules/local/gene_cnv_fold_changes_to_cbioportal'
include { GENERATE_CASE_LIST } from '../../../modules/local/generate_case_list'
include { GENERATE_META_FILE } from '../../../modules/local/generate_meta_file'

workflow GENOMIC_CNV {
    take:
        cnv_vcf // tuple (sample_id, filepath)
        ensembl_annotations
    main:

        cna_case_list = GENERATE_CASE_LIST(
            "cnv",
            cnv_vcf.map { it[0]}.collect().map{ it.sort(false).join('\t') } // item at index 0 is samplename, join all by tabs in order to send a list
        )

        fold_change_per_gene_cnv = EXTRACT_GENE_CNV_FOLD_CHANGES(
            cnv_vcf,
            ensembl_annotations
            )

        cbioportal_genomic_cnv_files = GENE_CNV_FOLD_CHANGES_TO_CBIOPORTAL(
            cnv_vcf,
            fold_change_per_gene_cnv
            )

        cbioportal_genomic_cnv_merged = cbioportal_genomic_cnv_files
            .collectFile( name : 'data_cna_hg38.seg', storeDir: "${params.outdir}", keepHeader : true, skip: 1, sort: 'deep')

        meta_text = """cancer_study_identifier: add_text
genetic_alteration_type: COPY_NUMBER_ALTERATION
datatype: SEG
reference_genome_id: hg38
description: Somatic CNA data (copy number segment file)
data_filename: data_cna_hg38.seg
        """

        GENERATE_META_FILE(
            "cna_hg38",
            meta_text
        )

    emit:
        cna_case_list
        cbioportal_genomic_cnv_merged

}
