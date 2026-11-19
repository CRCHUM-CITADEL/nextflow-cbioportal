process DOWNLOAD_MAFSMITH {
    // Release in the path: an Ensembl 113 download in assets/ must never be reused.
    storeDir "${projectDir}/assets/mafsmith/ensembl_116"
    container "${params.container_mafsmith}"

    input:
        val _ready   // gate: only run when there is work to do

    output:
        path "mafsmith", emit: data_dir, type: 'dir'

    script:
    """
    # Test-only bundle: Ensembl contigs (1, not chr1); fastvep comes from the container's PATH. See docs/usage.txt
    mafsmith fetch \\
        --data-dir mafsmith \\
        --genome grch38 \\
        --ensembl-release 116 \\
        --skip-fastvep
    """

    stub:
    """
    mkdir -p mafsmith
    """
}
