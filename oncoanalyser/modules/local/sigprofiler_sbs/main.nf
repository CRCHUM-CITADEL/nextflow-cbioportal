process SIGPROFILER_SBS {
    tag "$meta.sample"
    label 'process_low'

    container params.container_sigprofiler

    publishDir "${params.outdir}/${meta.group}/${meta.subject}", mode: 'copy'

    input:
        tuple val(meta), path(snv_counts)
        path cosmic_zip
        path sbs_metadata

    output:
        tuple val(meta), path("${meta.sample}.data_mutational_signatures_contribution_SBS.txt"), emit: sigs

    when:
        task.ext.when == null || task.ext.when

    script:
    """
    python3 -c "import zipfile, sys; zipfile.ZipFile(sys.argv[1]).extract('COSMIC_Human_SBS-96_GRCh38_v3.6.csv')" ${cosmic_zip}

    python3 ${projectDir}/bin/run_sigprofiler_sbs.py \\
        --snv_counts    ${snv_counts} \\
        --signatures_db COSMIC_Human_SBS-96_GRCh38_v3.6.csv \\
        --metadata      ${sbs_metadata} \\
        --sample        ${meta.sample} \\
        --output        ${meta.sample}.data_mutational_signatures_contribution_SBS.txt
    """

    stub:
    """
    touch ${meta.sample}.data_mutational_signatures_contribution_SBS.txt
    """
}
