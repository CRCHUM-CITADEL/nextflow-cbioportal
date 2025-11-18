include { GET_TPM } from '../../../modules/local/get_tpm'
include { MERGE_EXPRESSION_FILES_TO_CBIOPORTAL } from '../../../modules/local/merge_expression_files_to_cbioportal'
include { GENERATE_META_FILE } from '../../../modules/local/generate_meta_file'


workflow GENOMIC_EXPRESSION {
    take:
        somatic_expression // tuple (sample_id, filepath)
        gencode_annotations // gene annotation file

    main:

        all_groups = somatic_expression.map {meta, sample -> meta.group}.unique()

        tpm_file_ch = GET_TPM(
            somatic_expression,
            gencode_annotations
            )

        tpm_file_list = tpm_file_ch
            .map { meta, file -> tuple(meta.group, meta, file) }
            .groupTuple()
            .map { group, metas, files ->
                def meta = metas[0]  // Take first meta since they share the same group
                def sortedFiles = files.sort { a, b ->
                    def na = a.toString().split(/[\/\\]/).last()
                    def nb = b.toString().split(/[\/\\]/).last()
                    na <=> nb
                }
                if (sortedFiles.size() < 2) {
                    log.warn "GENOMIC_EXPRESSION: Found ${sortedFiles.size()} TPM file(s) for group ${group}. Need at least 2 files to merge. Skipping merge step."
                    return null
                }
                return tuple(meta, sortedFiles)
            }
            .filter { it != null }

        cbioportal_genomic_expression = MERGE_EXPRESSION_FILES_TO_CBIOPORTAL(
            tpm_file_list
            )

        meta_text = """cancer_study_identifier: add_text
genetic_alteration_type: MRNA_EXPRESSION
datatype: CONTINUOUS
stable_id: rna_seq_v2_mrna
show_profile_in_analysis_tab: true
profile_name: mRNA expression (RNA-Seq TPM)
profile_description: Expression levels (RNA-Seq TPM values)
data_filename: data_expression.txt
        """

        meta_file = GENERATE_META_FILE(
            all_groups,
            "expression",
            meta_text
        )

    emit:
        meta_file
        cbioportal_genomic_expression

}
