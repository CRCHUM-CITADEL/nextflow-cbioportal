include { EXTRACT_GENE_CNV_FOLD_CHANGES } from '../../../modules/local/extract_gene_cnv_fold_changes'
include { GENE_CNV_FOLD_CHANGES_TO_CBIOPORTAL } from '../../../modules/local/gene_cnv_fold_changes_to_cbioportal'

workflow GENOMIC_CNV {
    take:
        cnv_vcf // tuple (meta, filepath) 
        ensembl_annotations
    main:

        fold_change_per_gene_cnv = EXTRACT_GENE_CNV_FOLD_CHANGES(
            cnv_vcf,
            ensembl_annotations
            )

        cnv_vcf_with_fold_changes = cnv_vcf.join(fold_change_per_gene_cnv)

        cbioportal_genomic_cnv_files = GENE_CNV_FOLD_CHANGES_TO_CBIOPORTAL(
            cnv_vcf_with_fold_changes
            )

    emit:
        longfile   = cbioportal_genomic_cnv_files.long // channel [ long: [meta, seg_file], seg : [meta, long_file]]
        segfile    = cbioportal_genomic_cnv_files.seg
}
