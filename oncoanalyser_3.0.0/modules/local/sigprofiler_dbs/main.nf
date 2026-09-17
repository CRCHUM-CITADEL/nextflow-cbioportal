process SIGPROFILER_DBS {
    tag "$meta.sample"
    label 'process_low'

    container params.container_sigprofiler

    publishDir { "${params.outdir}/${meta.group}/${meta.subject}" }, mode: 'copy'

    input:
        tuple val(meta), path(somatic_vcf)
        path cosmic_zip
        path dbs_metadata

    output:
        tuple val(meta), path("${meta.sample}.data_mutational_signatures_contribution_DBS.txt"), emit: sigs_dbs
        tuple val(meta), path("${meta.sample}.data_mutational_signatures_counts_DBS.txt"), emit: sigs_counts_dbs

    when:
        task.ext.when == null || task.ext.when

    script:
    """
    python3 -c "import zipfile, sys; zipfile.ZipFile(sys.argv[1]).extract('COSMIC_Human_DBS-78_GRCh38_v3.6.csv')" ${cosmic_zip}

    python3 ${projectDir}/bin/run_sigprofiler_dbs.py \\
        --vcf             ${somatic_vcf} \\
        --signatures_db   COSMIC_Human_DBS-78_GRCh38_v3.6.csv \\
        --metadata        ${dbs_metadata} \\
        --sample          ${meta.sample} \\
        --output_contrib  ${meta.sample}.data_mutational_signatures_contribution_DBS.txt \\
        --output_counts   ${meta.sample}.data_mutational_signatures_counts_DBS.txt
    """

    stub:
    """
    touch ${meta.sample}.data_mutational_signatures_contribution_DBS.txt
    touch ${meta.sample}.data_mutational_signatures_counts_DBS.txt
    """
}
