process SIGPROFILER_SV {
    tag "$meta.sample"
    label 'process_low'

    container params.container_sigprofiler

    publishDir "${params.outdir}/${meta.group}/${meta.subject}", mode: 'copy'

    input:
        tuple val(meta), path(esvee_vcf)
        path cosmic_zip
        path sv_metadata

    output:
        tuple val(meta), path("${meta.sample}.data_mutational_signatures_contribution_SV.txt"), emit: sigs_sv
        tuple val(meta), path("${meta.sample}.data_mutational_signatures_counts_SV.txt"), emit: sigs_counts_sv

    when:
        task.ext.when == null || task.ext.when

    script:
    """
    python3 -c "import zipfile, sys; zipfile.ZipFile(sys.argv[1]).extract('COSMIC_Human_SV-32_GRCh38_v3.6.csv')" ${cosmic_zip}

    python3 ${projectDir}/bin/run_sigprofiler_sv.py \\
        --vcf             ${esvee_vcf} \\
        --signatures_db   COSMIC_Human_SV-32_GRCh38_v3.6.csv \\
        --metadata        ${sv_metadata} \\
        --sample          ${meta.sample} \\
        --output_contrib  ${meta.sample}.data_mutational_signatures_contribution_SV.txt \\
        --output_counts   ${meta.sample}.data_mutational_signatures_counts_SV.txt
    """

    stub:
    """
    touch ${meta.sample}.data_mutational_signatures_contribution_SV.txt
    touch ${meta.sample}.data_mutational_signatures_counts_SV.txt
    """
}
