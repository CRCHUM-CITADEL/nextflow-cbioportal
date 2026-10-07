process MERGE_SV {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(sv_files)

    output:
    path("data_sv.txt")

    script:
    """
    # Header-aware merge: per-sample files may differ in columns (older pipeline versions,
    # different writers), so rows are mapped onto the union of columns by name.
    merge_tsv_by_header.sh data_sv.txt *_sv.txt
    """
}
