process FILTER_MUTATIONS {
    tag "${meta.sample_id}"
    publishDir { "${params.outdir}/samples/${meta.sample_id}" }, mode: 'copy'

    input:
    tuple val(meta), path(maf), path(tsv)

    output:
    path("${meta.sample_id}_mutations.txt")

    script:
    """
    # Keep MAF rows overlapping a TSV row with Depth >= min_depth and VAF > min_vaf (N/A CNV/fusion rows never match)
    awk -v min_depth="${params.mutation_min_depth}" -v min_vaf="${params.mutation_min_vaf}" '
        function nochr(c) { sub(/^chr/, "", c); return c }
        function isnum(x) { return x ~ /^[0-9]+([.][0-9]*)?([eE][-+]?[0-9]+)?\$/ }
        BEGIN { FS = "\\t" }
        FILENAME == ARGV[1] {   # TSV regions (not NR==FNR: an empty TSV would swallow the MAF)
            sub(/\\r\$/, "")
            if (FNR == 1) {
                for (i = 1; i <= NF; i++) {
                    h = \$i; gsub(/^ +| +\$/, "", h)
                    if (h == "Chr") tc = i
                    else if (h == "Start") ts = i
                    else if (h == "End") te = i
                    else if (h == "Depth") td = i
                    else if (h == "VAF") tv = i
                }
                if (!tc || !ts || !te || !td || !tv) {
                    print "ERROR: TSV ${tsv} lacks one of Chr/Start/End/Depth/VAF (needed by --filter_tsv_variants)" > "/dev/stderr"
                    exit 1
                }
                next
            }
            if (\$ts !~ /^[0-9]+\$/ || \$te !~ /^[0-9]+\$/) next  # malformed rows
            if (!isnum(\$td) || !isnum(\$tv)) next  # N/A: CNV / fusion rows
            if (\$td + 0 < min_depth + 0 || \$tv + 0 <= min_vaf + 0) next  # fails QC
            chr = nochr(\$tc); k = ++n[chr]; rs[chr, k] = \$ts + 0; re[chr, k] = \$te + 0
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
