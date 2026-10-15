process ANNOTATE_GNOMAD {
    label 'python'

    input:
    path(mutations_file)

    output:
    path("gnomad/data_mutations.txt")

    script:
    """
    # Adds gnomADe_AF / gnomADg_AF from Genome Nexus; needs outbound HTTPS from this node.
    # Only chromosome, positions and alleles are sent, never sample IDs.
    mkdir -p gnomad
    annotate_gnomad.py ${mutations_file} --out gnomad/data_mutations.txt --url '${params.genome_nexus_url}'
    """
}
