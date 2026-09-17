process DOWNLOAD_MAFSMITH {
    storeDir "${projectDir}/assets"
    container "${params.container_mafsmith}"

    input:
        val _ready   // gate: only run when there is work to do

    output:
        path "mafsmith", emit: data_dir, type: 'dir'

    script:
    """
    # `mafsmith fetch` lays out <data-dir>/GRCh38/{genes.gff3.gz,reference.fa}, which is
    # what `mafsmith vcf2maf` reads out of \$HOME/.mafsmith. Ensembl release 113 matches
    # the BioMart TSV expected in params.ensembl_annotations.
    #
    # --skip-fastvep: fastvep is built into the container and on PATH, and mafsmith falls
    # back to `which fastvep` when <data-dir>/bin/fastvep is absent. Building it here
    # instead would need cargo plus a full toolchain in the task.
    #
    # NOTE: this downloads Ensembl's primary assembly, whose contigs are named 1, 2, X —
    # not chr1, chr2, chrX as in the GATK/hg38 reference oncoanalyser aligns against.
    # Pre-stage a bundle via --mafsmith_data for real cohorts; see docs/usage.txt.
    mafsmith fetch \\
        --data-dir mafsmith \\
        --genome grch38 \\
        --ensembl-release 113 \\
        --skip-fastvep
    """

    stub:
    """
    mkdir -p mafsmith
    """
}
