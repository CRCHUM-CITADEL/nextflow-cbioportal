/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PER_SAMPLE_FORMAT  } from '../subworkflows/local/per_sample_format/main'
include { MERGE_DEANON       } from '../subworkflows/local/merge_deanon/main'
include { STUDY_METADATA     } from '../subworkflows/local/study_metadata/main'
include { FILTER_LINKING     } from '../modules/local/filter_linking/main'
include { BUILD_ANON_LINKING } from '../modules/local/build_anon_linking/main'
include { PACKAGE_CBIOPORTAL } from '../modules/local/package_cbioportal/main'

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

    // -------------------------------------------------------------------------
    // Incremental skip: branch samples by whether all 4 per-sample outputs exist
    // -------------------------------------------------------------------------
    ch_samples_branched = ch_samples.branch { meta, tsv, folder ->
        existing: ['_sv.txt', '_cna.txt', '_seg.txt', '_mutations.txt'].every { suffix ->
            file("${params.outdir}/samples/${meta.sample_id}/${meta.sample_id}${suffix}").exists()
        }
        new_sample: true
    }

    ch_samples_branched.existing.subscribe { meta, tsv, folder ->
        log.info "Skipping already-processed sample: ${meta.sample_id}"
    }

    // -------------------------------------------------------------------------
    // Per-sample: format SV, CNA, mutations, and seg (new samples only)
    // -------------------------------------------------------------------------
    ch_tsv_new       = ch_samples_branched.new_sample.map { meta, tsv, folder -> tuple(meta, tsv) }
    ch_vcf_input_new = ch_samples_branched.new_sample.map { meta, tsv, folder -> tuple(meta, folder) }

    PER_SAMPLE_FORMAT(ch_tsv_new, ch_vcf_input_new)

    // -------------------------------------------------------------------------
    // Read existing per-sample outputs from disk (already-processed samples)
    // -------------------------------------------------------------------------
    ch_existing_sv        = ch_samples_branched.existing.map { meta, tsv, folder ->
        file("${params.outdir}/samples/${meta.sample_id}/${meta.sample_id}_sv.txt")
    }
    ch_existing_cna       = ch_samples_branched.existing.map { meta, tsv, folder ->
        file("${params.outdir}/samples/${meta.sample_id}/${meta.sample_id}_cna.txt")
    }
    ch_existing_mutations = ch_samples_branched.existing.map { meta, tsv, folder ->
        file("${params.outdir}/samples/${meta.sample_id}/${meta.sample_id}_mutations.txt")
    }
    ch_existing_seg       = ch_samples_branched.existing.map { meta, tsv, folder ->
        file("${params.outdir}/samples/${meta.sample_id}/${meta.sample_id}_seg.txt")
    }

    // -------------------------------------------------------------------------
    // Filter linking file to only samples in the samplesheet
    // -------------------------------------------------------------------------
    ch_linking = Channel.value(file(params.linking_file))

    ch_samplesheet_ids = ch_samplesheet
        .map { row -> row.subject_id.toUpperCase() }
        .collectFile(name: 'samplesheet_ids.txt', newLine: true)

    FILTER_LINKING(ch_samplesheet_ids, ch_linking)
    ch_filtered_linking = FILTER_LINKING.out

    // -------------------------------------------------------------------------
    // When anonymize=true, build an anonymised linking file where
    // deanon_sample_id and deanon_patient_id = lowercase of sample_id.
    // The existing DEANON scripts then map to anonymised IDs instead of real IDs.
    // -------------------------------------------------------------------------
    if (params.anonymize) {
        BUILD_ANON_LINKING(ch_filtered_linking)
        ch_output_linking = BUILD_ANON_LINKING.out
    } else {
        ch_output_linking = ch_filtered_linking
    }

    // -------------------------------------------------------------------------
    // Collect: merge new + existing outputs, then remap sample IDs
    // -------------------------------------------------------------------------
    MERGE_DEANON(
        PER_SAMPLE_FORMAT.out.sv.mix(ch_existing_sv).collect(),
        PER_SAMPLE_FORMAT.out.cna.mix(ch_existing_cna).collect(),
        PER_SAMPLE_FORMAT.out.mutations.mix(ch_existing_mutations).collect(),
        PER_SAMPLE_FORMAT.out.seg.mix(ch_existing_seg).collect(),
        ch_output_linking
    )

    // -------------------------------------------------------------------------
    // Once: clinical files, case lists, and meta files
    // -------------------------------------------------------------------------
    STUDY_METADATA(
        Channel.fromPath(params.patient_file),
        Channel.fromPath(params.sample_file),
        ch_filtered_linking,
        ch_output_linking,
        params.study_id
    )

    // -------------------------------------------------------------------------
    // Package all cBioPortal files into a tar.gz for sharing between instances
    // -------------------------------------------------------------------------
    all_data_files = MERGE_DEANON.out.mutations
        .mix(MERGE_DEANON.out.sv)
        .mix(MERGE_DEANON.out.cna)
        .mix(MERGE_DEANON.out.seg)
        .mix(STUDY_METADATA.out.clinical_patient)
        .mix(STUDY_METADATA.out.clinical_sample)
        .collect()

    PACKAGE_CBIOPORTAL(
        params.study_id,
        all_data_files,
        STUDY_METADATA.out.meta_files,
        STUDY_METADATA.out.case_lists
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
