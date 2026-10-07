process MERGE_SEG {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(seg_files)

    output:
    path("data_seg.txt")

    script:
    """
    # Header-aware merge: per-sample files may differ in columns (older pipeline versions,
    # different writers), so rows are mapped onto the union of columns by name.
    merge_tsv_by_header.sh data_seg.txt *_seg.txt
    """
}
