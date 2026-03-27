process INDEX_FASTA {

    label 'process_low'

    // take apptainer container from nf-core https://github.com/nf-core/modules/blob/master/modules/nf-core/samtools/view/main.nf
    container 'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0'

    input:
    path(fasta)

    output:
    path("${fasta}.fai")

    script:
    """
    samtools faidx ${fasta}
    """
}
