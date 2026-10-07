process MERGE_SAMPLE_SV {
    tag "${meta.sample}"
    label 'process_single'

    publishDir { "${params.outdir}/${meta.group}/${meta.subject}" }, mode: 'copy'

    input:
        tuple val(meta), path(sv_files)

    output:
        tuple val(meta), path("${meta.sample}.data_sv.txt"), emit: sv

    when:
        task.ext.when == null || task.ext.when

    script:
    // Sort by name: the input channel's order varies between runs, which would
    // otherwise make this file (and the group-level data_sv.txt) non-reproducible.
    def files = (sv_files instanceof List ? sv_files : [sv_files]).sort { it.name }
    def data_cmds = files.collect { f -> "tail -n +2 ${f}" }.join('; ')
    """
    head -n 1 ${files[0]} > "tmp_sv_merge.txt"
    { ${data_cmds}; } | awk '!seen[\$0]++' >> "tmp_sv_merge.txt"
    mv "tmp_sv_merge.txt" "${meta.sample}.data_sv.txt"
    """

    stub:
    """
    touch "${meta.sample}.data_sv.txt"
    """
}
