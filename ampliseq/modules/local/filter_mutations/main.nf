process FILTER_MUTATIONS {
    tag "${meta.sample_id}"
    publishDir { "${params.outdir}/samples/${meta.sample_id}" }, mode: 'copy'

    input:
    tuple val(meta), path(maf), path(tsv)

    output:
    path("${meta.sample_id}_mutations.txt")

    script:
    """
    # Keep MAF rows whose [Start_Position, End_Position] overlaps any TSV region (\$1=Chr
    # \$2=Start \$3=End, inclusive) on the same chromosome. TSV rows are regions (gene spans,
    # CNV/SV intervals), so an exact coordinate match would drop nearly every SNV/indel.
    # MAF columns are found by header name; the chr prefix is ignored on both sides.
    awk '
        function nochr(c) { sub(/^chr/, "", c); return c }
        BEGIN { FS = "\\t" }
        FILENAME == ARGV[1] {   # TSV regions (not NR==FNR: an empty TSV would swallow the MAF)
            if (\$2 !~ /^[0-9]+\$/ || \$3 !~ /^[0-9]+\$/) next   # header / malformed rows
            chr = nochr(\$1); k = ++n[chr]; rs[chr, k] = \$2 + 0; re[chr, k] = \$3 + 0
            next
        }
        FNR==1  {
            for (i = 1; i <= NF; i++) {
                if (\$i == "Chromosome") c = i
                else if (\$i == "Start_Position") s = i
                else if (\$i == "End_Position") e = i
            }
            if (!c || !s || !e) { print "ERROR: MAF lacks Chromosome/Start_Position/End_Position" > "/dev/stderr"; exit 1 }
            print; next
        }
        {
            chr = nochr(\$c); ms = \$s + 0; me = \$e + 0
            for (k = 1; k <= n[chr]; k++) if (ms <= re[chr, k] && me >= rs[chr, k]) { print; next }
        }
    ' "${tsv}" ${maf} > "${meta.sample_id}_mutations.txt"
    """
}
