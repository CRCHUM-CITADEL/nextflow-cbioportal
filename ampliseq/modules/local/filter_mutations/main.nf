process FILTER_MUTATIONS {
    tag "${meta.sample_id}"
    publishDir "${params.outdir}/samples/${meta.sample_id}", mode: 'copy'

    input:
    tuple val(meta), path(maf), path(tsv)

    output:
    path("${meta.sample_id}_mutations.txt")

    script:
    """
    # Build region set from TSV (\$1=Chr \$2=Start \$3=End). MAF columns are found by header
    # name, and the chr prefix is ignored on both sides so either contig naming matches.
    awk '
<<<<<<< HEAD
        NR==FNR { regions["chr" \$1 ":" \$2 "-" \$3] = 1; next }
        FNR==1  { print; next }
        \$5 ":" \$6 "-" \$7 in regions
    ' "${tsv}" ${maf} > ${meta.sample_id}_mutations_non_filtered.txt

    filter_mutations.R ${meta.sample_id}_mutations_non_filtered.txt ${meta.sample_id}_mutations.txt 
=======
        function nochr(c) { sub(/^chr/, "", c); return c }
        BEGIN { FS = "\\t" }
        NR==FNR { regions[nochr(\$1) ":" \$2 "-" \$3] = 1; next }
        FNR==1  {
            for (i = 1; i <= NF; i++) {
                if (\$i == "Chromosome") c = i
                else if (\$i == "Start_Position") s = i
                else if (\$i == "End_Position") e = i
            }
            if (!c || !s || !e) { print "ERROR: MAF lacks Chromosome/Start_Position/End_Position" > "/dev/stderr"; exit 1 }
            print; next
        }
        (nochr(\$c) ":" \$s "-" \$e) in regions
    ' "${tsv}" ${maf} > "${meta.sample_id}_mutations.txt"
>>>>>>> 3015a187e1e91a401619907b44ab9d15742d1b70
    """
}
