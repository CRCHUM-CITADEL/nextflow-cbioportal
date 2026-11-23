process PASSTHROUGH_MUTATIONS {
    tag "${meta.sample_id}"
    label 'python'
    publishDir { "${params.outdir}/samples/${meta.sample_id}" }, mode: 'copy'

    input:
    tuple val(meta), path(maf)
    path(gene_panel)  // gene panel TSV, or [] when unset

    output:
    path("${meta.sample_id}_mutations.txt")

    script:
    """
    cp "${maf}" "${meta.sample_id}_mutations.txt"
    if [ -n "${gene_panel}" ]; then filter_gene_panel.py "${meta.sample_id}_mutations.txt" --panel "${gene_panel}" --type mutation --out "${meta.sample_id}_mutations.txt.tmp" && mv "${meta.sample_id}_mutations.txt.tmp" "${meta.sample_id}_mutations.txt"; fi
    """
}
