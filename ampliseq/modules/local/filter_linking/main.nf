process FILTER_LINKING {
    input:
    path(pairs_file)   // TSV: SUBJECT_ID(uppercase)\tsample_id (no header)
    path(linking_file)

    output:
    path("linking_filtered.txt")

    script:
    // Join samplesheet pairs with linking file on subject_id.
    // Linking file sample_id column is the subject-level anon ID (uppercase).
    // Output replaces it with the actual sample_id from the samplesheet.
    """
    awk '
        NR==FNR { subj_to_sample[\$1]=\$2; next }
        FNR==1  { print; next }
        \$1 in subj_to_sample { print subj_to_sample[\$1] "\\t" \$2 "\\t" \$3 }
    ' ${pairs_file} ${linking_file} > linking_filtered.txt
    """
}
