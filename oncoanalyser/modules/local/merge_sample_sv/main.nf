process MERGE_SAMPLE_SV {
    tag "${meta.sample}"
    label 'process_single'

    publishDir "${params.outdir}/${meta.group}/${meta.subject}", mode: 'copy'

    input:
        tuple val(meta), path(sv_files)

    output:
        tuple val(meta), path("${meta.sample}.data_sv.txt"), emit: sv

    when:
        task.ext.when == null || task.ext.when

    script:
    def files = sv_files instanceof List ? sv_files : [sv_files]
    """
    head -n 1 ${files[0]} > "${meta.sample}.data_sv.txt"
    ${files.collect { "tail -n +2 ${it} >> '${meta.sample}.data_sv.txt'" }.join('\n')}
    """

    stub:
    """
    touch "${meta.sample}.data_sv.txt"
    """
}
