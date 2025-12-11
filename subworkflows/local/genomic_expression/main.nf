include { GET_TPM                              } from '../../../modules/local/get_tpm'
include { MERGE_EXPRESSION_FILES_TO_CBIOPORTAL } from '../../../modules/local/merge_expression_files_to_cbioportal'
include { GENERATE_META_FILE                   } from '../../../modules/local/generate_meta_file'


workflow GENOMIC_EXPRESSION {
    take:
        somatic_expression // tuple (meta, filepath)
        ensembl_annotations_expr // gene annotation file

    main:

        all_groups = somatic_expression.map {meta, sample -> meta.group}.unique()

        tpm_file_ch = GET_TPM(
            somatic_expression,
            ensembl_annotations_expr
            )

    emit:
        out = tpm_file_ch
}
