process FORMAT_CNA {
    tag "${meta.sample_id}"
    label 'python'
    publishDir { "${params.outdir}/samples/${meta.sample_id}" }, mode: 'copy'

    input:
    tuple val(meta), path(tsv)

    output:
    path("${meta.sample_id}_cna.txt")

    script:
    // Filter mode: HIGH-confidence CN >= 6 only (toBoolean: NF 26 CLI params are strings)
    def filter_args = params.filter_tsv_variants.toString().toBoolean()
        ? "--min-copy-number ${params.cna_min_copy_number} --confidence '${params.cna_confidence}'"
        : ''
    """
    format_cna.py "${tsv}" "${meta.sample_id}" ${filter_args}
    mv data_cna.txt "${meta.sample_id}_cna.txt"
    """
}
