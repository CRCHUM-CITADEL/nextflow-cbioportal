process FILTER_LINKING {
    input:
    path(sample_ids_file) // one anonymized sample_id per line (no header)
    path(linking_file)

    output:
    path("linking_filtered.txt")

    script:
    """
    awk '
        NR==FNR { ids[\$1]=1; next }
        FNR==1  { print; next }
        \$1 in ids
    ' ${sample_ids_file} ${linking_file} > linking_filtered.txt
    """
}
