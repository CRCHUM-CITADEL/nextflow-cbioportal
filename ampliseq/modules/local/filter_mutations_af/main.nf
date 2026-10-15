process FILTER_MUTATIONS_AF {
    label 'python'

    input:
    path(mutations_file)

    output:
    path("filtered/data_mutations.txt")

    script:
    def gnomad_args = params.filter_gnomad_af.toString().toBoolean() ? "--max-gnomad-af ${params.max_gnomad_af} --gnomad-columns '${params.gnomad_af_columns}'" : ''
    def vaf_args    = params.filter_vaf.toString().toBoolean() ? "--max-vaf ${params.max_vaf}" : ''
    """
    mkdir -p filtered
    filter_mutations_af.py ${mutations_file} --out filtered/data_mutations.txt ${gnomad_args} ${vaf_args}
    """
}
