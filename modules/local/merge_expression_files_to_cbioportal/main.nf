process MERGE_EXPRESSION_FILES_TO_CBIOPORTAL {
    publishDir "${params.outdir}/${meta.group}", mode: 'copy'

    container params.container_r

    input:
        tuple val(meta), path(tpm_file_list)

    output:
        tuple val(meta), path("data_expression.txt")

    script:
    """
    gen_merge_expression_files_to_cbioportal.R \
    --input_files ${tpm_file_list.join(',')} \
    --output_file data_expression.txt \
    --fill_missing 0 \
    --strict
    """
}
