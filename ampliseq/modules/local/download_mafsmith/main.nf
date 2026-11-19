process DOWNLOAD_MAFSMITH {
    label "mafsmith"
    storeDir "${projectDir}/assets"

    input:
    val _ready // gate: only run when there is work to do

    output:
    path "mafsmith", emit: data_dir, type: 'dir'

    script:
    """
    # Lays out mafsmith/GRCh37/{genes.gff3.gz,reference.fa}. For GRCh37 mafsmith 0.1.0 ignores
    # --ensembl-release: it always fetches Homo_sapiens.GRCh37.87 (frozen, no canonical tags).
    # fastVEP is built into the container.
    # NOTE: Ensembl's primary assembly names contigs 1, 2, X (no chr prefix); MAFSMITH
    # restores the VCF's naming on output.
    mafsmith fetch \\
        --data-dir mafsmith \\
        --genome grch37 \\
        --ensembl-release 113 \\
        --skip-fastvep
    """

    stub:
    """
    mkdir -p mafsmith
    """
}
