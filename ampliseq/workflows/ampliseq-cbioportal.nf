/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PER_SAMPLE_FORMAT } from '../subworkflows/local/per_sample_format/main'
include { MERGE_DEANON      } from '../subworkflows/local/merge_deanon/main'
include { STUDY_METADATA    } from '../subworkflows/local/study_metadata/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow AMPLISEQ_CBIOPORTAL {

    take:
    ch_samplesheet // channel: rows from samplesheet CSV

    main:

    // -------------------------------------------------------------------------
    // Build per-sample channel: tuple(meta, tsv, sample_folder)
    // -------------------------------------------------------------------------
    ch_samples = ch_samplesheet
        .map { row ->
            def meta = [
                group     : row.group,
                subject_id: row.subject_id,
                sample_id : row.sample_id,
            ]
            def folder = file(row.folder_location)
            if (!folder.exists()) {
                error "folder_location does not exist: ${folder}"
            }
            def tsvList = file("${row.folder_location}/analysis_*_export.tsv", glob: true)
            if (!tsvList) {
                error "No analysis_*_export.tsv found in: ${folder}"
            }
            def tsv = tsvList[0]
            tuple(meta, tsv, folder)
        }

    ch_tsv       = ch_samples.map { meta, tsv, folder -> tuple(meta, tsv) }
    ch_vcf_input = ch_samples.map { meta, tsv, folder -> tuple(meta, folder) }

    // -------------------------------------------------------------------------
    // Per-sample: format SV, CNA, and mutations
    // -------------------------------------------------------------------------
    PER_SAMPLE_FORMAT(ch_tsv, ch_vcf_input)

    // -------------------------------------------------------------------------
    // Collect: merge and deanonymise all data types
    // -------------------------------------------------------------------------
    ch_linking = Channel.value(file(params.linking_file))

    MERGE_DEANON(
        PER_SAMPLE_FORMAT.out.sv.collect(),
        PER_SAMPLE_FORMAT.out.cna.collect(),
        PER_SAMPLE_FORMAT.out.mutations.collect(),
        ch_linking
    )

    // -------------------------------------------------------------------------
    // Once: clinical files, case lists, and meta files
    // -------------------------------------------------------------------------
    STUDY_METADATA(
        Channel.fromPath(params.patient_file),
        Channel.fromPath(params.sample_file),
        ch_linking,
        params.study_id
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
