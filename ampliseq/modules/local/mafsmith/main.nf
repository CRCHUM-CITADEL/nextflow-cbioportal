process MAFSMITH {
    tag "${meta.sample_id}"
    label "mafsmith"

    input:
    tuple val(meta), path(sample_folder)
    path(mafsmith_data) // mafsmith home: <data>/GRCh37/{reference.fa,genes.gff3.gz}

    output:
    tuple val(meta), path("${meta.sample_id}.maf")

    when:
    task.ext.when == null || task.ext.when

    script:
    def retain_ann = params.mafsmith_retain_ann ? "--retain-ann ${params.mafsmith_retain_ann}" : ''
    def args       = task.ext.args ?: ''
    """
    VCF=\$(find -L "${sample_folder}" -maxdepth 1 \\( -name '*-basespace-pisces.final.vcf.gz' -o -name '*-basespace-pisces.final.vcf' \\) | head -1)
    [ -n "\$VCF" ] || { echo "ERROR: No VCF found in ${sample_folder}" >&2; exit 1; }

    # An empty bundle (e.g. left in assets/mafsmith by a -stub-run of DOWNLOAD_MAFSMITH) would
    # otherwise be reused by every later run; fail loudly unless annotation is skipped.
    if [[ " ${args} " != *" --skip-annotation "* ]]; then
        for f in GRCh37/reference.fa GRCh37/genes.gff3.gz; do
            [ -s "${mafsmith_data}/\$f" ] || { echo "ERROR: ${mafsmith_data}/\$f missing or empty; delete the bundle (e.g. assets/mafsmith) to re-fetch, or fix --mafsmith_data" >&2; exit 1; }
        done
    fi

    # The bundle uses mafsmith fetch (Ensembl) naming: 1, 2, X, MT. mafsmith refuses a VCF whose
    # contigs are not in the FASTA, and Pisces VCFs say chr1, so strip chr (chrM -> MT) from the
    # records and ##contig lines. Keep the header and PASS records only.
    VCF_HAS_CHR=0
    zcat -f "\$VCF" | grep -v '^#' | head -1 | grep -q '^chr' && VCF_HAS_CHR=1 || true

    zcat -f "\$VCF" | awk '
        function no_chr(c) { if (c !~ /^chr/) return c; sub(/^chr/, "", c); return (c == "M" ? "MT" : c) }
        BEGIN { FS = OFS = "\\t" }
        /^##contig=<ID=chrM[,>]/ { sub(/ID=chrM/, "ID=MT"); print; next }
        /^##contig=<ID=chr/      { sub(/ID=chr/, "ID="); print; next }
        /^#/ { print; next }
        \$7 == "PASS" { \$1 = no_chr(\$1); print }
    ' > "${meta.sample_id}.pass.vcf"

    # mafsmith looks under \$HOME/.mafsmith for its reference data
    export HOME=./
    ln -s ${mafsmith_data} .mafsmith

    # Ampliseq VCFs are tumor-only: the single sample column is the tumor
    mafsmith vcf2maf \\
        --genome grch37 \\
        --tumor-id "${meta.sample_id}" \\
        ${retain_ann} \\
        ${args} \\
        --input-vcf "${meta.sample_id}.pass.vcf" \\
        --output-maf "tmp.${meta.sample_id}.maf"

    # Downstream expects the column header on line 1 (drop mafsmith's #version line) and the
    # original VCF's contig naming, so put chr back (MT -> chrM) for chr-prefixed VCFs.
    grep -v '^#' "tmp.${meta.sample_id}.maf" \\
        | awk -v addchr="\$VCF_HAS_CHR" 'BEGIN { FS = OFS = "\\t" }
            NR == 1 {
                for (i = 1; i <= NF; i++) if (\$i == "Chromosome") col = i
                if (!col) { print "ERROR: no Chromosome column in mafsmith output" > "/dev/stderr"; exit 1 }
                print; next
            }
            addchr == 1 && \$col !~ /^chr/ { \$col = (\$col == "MT" ? "chrM" : "chr" \$col) }
            { print }' > "${meta.sample_id}.maf"
    """

    stub:
    """
    printf 'Hugo_Symbol\\tEntrez_Gene_Id\\tCenter\\tNCBI_Build\\tChromosome\\tStart_Position\\tEnd_Position\\tStrand\\tVariant_Classification\\tVariant_Type\\tReference_Allele\\tTumor_Seq_Allele1\\tTumor_Seq_Allele2\\tdbSNP_RS\\tdbSNP_Val_Status\\tTumor_Sample_Barcode\\n' > "${meta.sample_id}.maf"
    """
}
