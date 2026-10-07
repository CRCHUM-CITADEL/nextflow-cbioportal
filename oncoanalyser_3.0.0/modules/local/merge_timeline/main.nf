process MERGE_TIMELINE {
    publishDir { "${params.outdir}/${group}" }, mode: 'copy'

    container params.container_r

    tag { group }

    input:
        tuple val(group), path(timeline_parts, stageAs: "part_*.txt")

    output:
        tuple val(group), path("data_timeline.txt"), emit: ch_timeline, optional: true

    script:
    def inputs_arg = timeline_parts instanceof List ? timeline_parts.join(",") : timeline_parts
    """
    merge_timeline.R --inputs ${inputs_arg}
    """

    stub:
    """
    touch data_timeline.txt
    """
}
