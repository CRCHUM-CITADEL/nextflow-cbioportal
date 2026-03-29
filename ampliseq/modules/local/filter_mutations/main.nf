process FILTER_MUTATIONS {
    tag "${meta.sample_id}"

    input:
    tuple val(meta), path(maf), path(tsv)

    output:
    path("${meta.sample_id}_mutations.txt")

    script:
    """
    # Build region set from TSV (Chr col has no chr prefix; MAF Chromosome col has chr prefix)
    # TSV: \$1=Chr \$2=Start \$3=End  |  MAF: \$5=Chromosome \$6=Start_Position \$7=End_Position
    awk '
        NR==FNR { regions["chr" \$1 ":" \$2 "-" \$3] = 1; next }
        FNR==1  { print; next }
        \$5 ":" \$6 "-" \$7 in regions
    ' "${tsv}" ${maf} > "${meta.sample_id}_mutations.txt"
    """
}
