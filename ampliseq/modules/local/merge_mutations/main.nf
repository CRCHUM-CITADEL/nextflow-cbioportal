process MERGE_MUTATIONS {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(mutation_files)

    output:
    path("data_mutations.txt")

    script:
    """
    # Header-aware merge: per-sample MAFs can have different column sets (e.g. older
    # vcf2maf runs vs mafsmith). Pass 1 builds the union of columns in first-seen order,
    # pass 2 writes each row mapped onto it, leaving missing columns blank.
    files=( *_mutations.txt )
    awk -F'\\t' -v OFS='\\t' -v nf=\${#files[@]} '
        FNR==1 { fileno++ }
        fileno <= nf {
            if (FNR==1) for (i = 1; i <= NF; i++) if (!(\$i in seen)) { seen[\$i] = 1; cols[++n] = \$i }
            next
        }
        FNR==1 {
            if (!printed) { out = cols[1]; for (j = 2; j <= n; j++) out = out OFS cols[j]; print out; printed = 1 }
            split("", idx); for (i = 1; i <= NF; i++) idx[\$i] = i
            next
        }
        {
            out = ""
            for (j = 1; j <= n; j++) out = out (j > 1 ? OFS : "") ((cols[j] in idx) ? \$(idx[cols[j]]) : "")
            print out
        }
    ' "\${files[@]}" "\${files[@]}" > data_mutations.txt
    """
}
