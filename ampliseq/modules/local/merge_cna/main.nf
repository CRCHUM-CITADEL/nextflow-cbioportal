process MERGE_CNA {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(cna_files)

    output:
    path("data_cna.txt")

    script:
    """
    # Header-aware merge: per-sample files may differ in columns (older pipeline versions,
    # different writers), so rows are mapped onto the union of columns by name.
    merge_tsv_by_header.sh data_cna.txt *_cna.txt
    """
}
