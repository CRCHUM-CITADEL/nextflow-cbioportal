process PREPARE_MAFSMITH {
    label "mafsmith"

    input:
    path mafsmith_data, stageAs: 'mafsmith_bundle' // mafsmith home: <data>/GRCh37/{reference.fa,genes.gff3.gz}

    output:
    path 'mafsmith_bundle', includeInputs: true, emit: data_dir

    script:
    """
    # Runs once per pipeline run, before any MAFSMITH task, and completes the bundle in place:
    #  - mafsmith needs GRCh37/reference.fa.fai, which `mafsmith fetch` never writes
    #    ("Cannot read FASTA index"); no samtools in the image, so build it with awk.
    #  - mafsmith decompresses genes.gff3.gz to genes.gff3 on first use, non-atomically and
    #    only if genes.gff3 is absent; parallel MAFSMITH tasks would otherwise read a
    #    half-written GFF and silently lose annotation.
    # Both are written to a temp name and renamed, so an interrupted run leaves no partial file.
    G="${mafsmith_data}/GRCh37"
    [ -s "\$G/reference.fa" ] || { echo "ERROR: \$G/reference.fa missing or empty; delete the bundle (e.g. assets/mafsmith) to re-fetch, or fix --mafsmith_data" >&2; exit 1; }
    [ -s "\$G/genes.gff3.gz" ] || [ -s "\$G/genes.gff3" ] || { echo "ERROR: \$G/genes.gff3.gz missing or empty" >&2; exit 1; }

    if [ ! -s "\$G/reference.fa.fai" ]; then
        # samtools-faidx-compatible (uniform line length): NAME LENGTH OFFSET LINEBASES LINEWIDTH
        LC_ALL=C awk '
            function out() { printf "%s\\t%.0f\\t%.0f\\t%d\\t%d\\n", name, len, off, lb, lw }
            /^>/ {
                if (name != "") out()
                name = substr(\$1, 2); pos += length(\$0) + 1; off = pos; len = 0; lb = 0; lw = 0; next
            }
            { if (lb == 0) { lb = length(\$0); lw = lb + 1 } len += length(\$0); pos += length(\$0) + 1 }
            END { if (name != "") out() }
        ' "\$G/reference.fa" > "\$G/reference.fa.fai.tmp.\$\$"
        mv "\$G/reference.fa.fai.tmp.\$\$" "\$G/reference.fa.fai"
    fi

    # A genes.gff3 whose size differs from the decompressed .gz was left half-written (or is
    # stale): redo it. The new mtime also makes fastVEP rebuild a cache built from the bad copy.
    if [ -s "\$G/genes.gff3.gz" ]; then
        WANT=\$(zcat "\$G/genes.gff3.gz" | wc -c)
        HAVE=\$(wc -c < "\$G/genes.gff3" 2>/dev/null || echo 0)
        if [ "\$WANT" != "\$HAVE" ]; then
            echo "Decompressing genes.gff3 (have \$HAVE bytes, want \$WANT)" >&2
            zcat "\$G/genes.gff3.gz" > "\$G/genes.gff3.tmp.\$\$"
            mv "\$G/genes.gff3.tmp.\$\$" "\$G/genes.gff3"
        fi
    fi
    """

    stub:
    """
    true
    """
}
