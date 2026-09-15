include { BUILD_CLINICAL_TABLE }                    from '../../../modules/local/build_clinical_table'
include { WRITE_CLINICAL_SAMPLE }                    from '../../../modules/local/write_clinical_sample'
include { WRITE_CLINICAL_PATIENT }                   from '../../../modules/local/write_clinical_patient'
include { SPLIT_CLINICAL }                           from '../../../modules/local/split_clinical'
include { GENERATE_META_FILE }                       from '../../../modules/local/generate_meta_file'
include { GENERATE_META_FILE as GENERATE_META_FILE_TIMELINE } from '../../../modules/local/generate_meta_file'
include { GENERATE_TIMELINE_SURGERY }                from '../../../modules/local/generate_timeline_surgery'
include { GENERATE_TIMELINE_TREATMENT }              from '../../../modules/local/generate_timeline_treatment'
include { GENERATE_TIMELINE_STATUS }                 from '../../../modules/local/generate_timeline_status'
include { GENERATE_TIMELINE_SPECIMEN }               from '../../../modules/local/generate_timeline_specimen'
include { GENERATE_TIMELINE_LAB_TEST }               from '../../../modules/local/generate_timeline_lab_test'
include { MERGE_TIMELINE }                           from '../../../modules/local/merge_timeline'

workflow CLINICAL_AGGREGATE {
    take:
        filelist          // channel: [meta(group, pipeline, extraction_date), csv_path]
        genomic_subjects  // val: path to genomic subjects TSV, or "" to skip filtering

    main:

        ch_clinical_common = Channel.fromPath("${projectDir}/bin/clinical_common.R").first()

        csvs = filelist
            .map { meta, csv ->
                tuple(meta.group, [(meta.pipeline): csv])
            }
            .groupTuple()
            .map { group, data_list ->
                tuple(group, data_list.collectEntries())
            }
            .combine(genomic_subjects)
            .map { group, csv_map, gs ->
                if (gs) csv_map.genomic_subjects = gs
                tuple(group, csv_map)
            }

        all_groups = csvs.map { group, csv_map -> group }.unique()

        // ── Clinical sample/patient files ────────────────────────────────────
        // BUILD_CLINICAL_TABLE runs the loaders + merge cascade once per group;
        // the two writers below are thin column-selection passes over its output,
        // so no work is duplicated the way the old two-mode FORMAT_CLINICAL was.
        csvs
            .map { group, csv_map ->
                return tuple(
                    [group: group],
                    csv_map.donors               ? file(csv_map.donors)               : [],
                    csv_map.primary_diagnoses    ? file(csv_map.primary_diagnoses)    : [],
                    csv_map.specimens            ? file(csv_map.specimens)            : [],
                    csv_map.sample_registrations ? file(csv_map.sample_registrations) : [],
                    csv_map.treatments           ? file(csv_map.treatments)           : [],
                    csv_map.surgeries            ? file(csv_map.surgeries)            : [],
                    csv_map.systemic_therapies   ? file(csv_map.systemic_therapies)   : [],
                    csv_map.radiations           ? file(csv_map.radiations)           : [],
                    csv_map.follow_ups           ? file(csv_map.follow_ups)           : [],
                    csv_map.biomarkers           ? file(csv_map.biomarkers)           : [],
                    csv_map.genomic_subjects     ? file(csv_map.genomic_subjects)     : [],
                    file(params.mohccn_primary_site_map),
                    file(params.mohccn_specimen_tissue_source_map),
                    file(params.mohccn_treatment_intent_map)
                )
            }
            .set { ch_build_input }

        BUILD_CLINICAL_TABLE(ch_build_input, ch_clinical_common)

        WRITE_CLINICAL_SAMPLE(BUILD_CLINICAL_TABLE.out.ch_merged)
        WRITE_CLINICAL_PATIENT(BUILD_CLINICAL_TABLE.out.ch_merged)

        // ── Per-sample clinical split (both mode only) ───────────────────────
        ch_group_clinical = WRITE_CLINICAL_SAMPLE.out.ch_clinical_sample
            .join(WRITE_CLINICAL_PATIENT.out.ch_clinical_patient)

        // Parse linking file for per-sample tuples (empty channel in clinical-only mode)
        ch_per_sample = genomic_subjects
            .filter { it && it != "" }
            .flatMap { gs_path ->
                file(gs_path).readLines().drop(1).findAll { it.trim() }.collect { line ->
                    def fields = line.split('\t')
                    [fields[0], fields[1]]
                }
            }

        ch_split_input = ch_per_sample
            .combine(ch_group_clinical)
            .map { subject, sample, group, clin_sample, clin_patient ->
                tuple(
                    [group: group, subject: subject, sample: sample],
                    clin_sample,
                    clin_patient
                )
            }

        SPLIT_CLINICAL(ch_split_input)

        meta_text = Channel.of("""cancer_study_identifier: add_text
genetic_alteration_type: CLINICAL
datatype: SAMPLE_ATTRIBUTES
data_filename: data_clinical_sample.txt
        """,
        """cancer_study_identifier: add_text
genetic_alteration_type: CLINICAL
datatype: PATIENT_ATTRIBUTES
data_filename: data_clinical_patient.txt
        """)

        file_names = Channel.of("clinical_sample", "clinical_patient")

        all_groups_times_two = all_groups.combine(file_names).map { g, f -> g }

        GENERATE_META_FILE(
            all_groups_times_two,
            file_names,
            meta_text
        )

        // ── Timeline files ────────────────────────────────────────────────────
        // Five per-EVENT_TYPE processes (each optional-output, only surfacing a
        // part when its source CSV produced rows), unioned by MERGE_TIMELINE.
        // Column ordering there is a hard contract with combine_cbioportal_outputs.py
        // — see merge_timeline.R.
        csvs
            .map { group, csv_map ->
                tuple(
                    [group: group],
                    csv_map.sample_registrations ? file(csv_map.sample_registrations) : [],
                    csv_map.treatments           ? file(csv_map.treatments)           : [],
                    csv_map.surgeries            ? file(csv_map.surgeries)            : [],
                    csv_map.genomic_subjects     ? file(csv_map.genomic_subjects)     : [],
                    file(params.mohccn_primary_site_map),
                    file(params.mohccn_treatment_intent_map)
                )
            }
            .set { ch_surgery_input }
        GENERATE_TIMELINE_SURGERY(ch_surgery_input, ch_clinical_common)

        csvs
            .map { group, csv_map ->
                tuple(
                    [group: group],
                    csv_map.sample_registrations ? file(csv_map.sample_registrations) : [],
                    csv_map.treatments           ? file(csv_map.treatments)           : [],
                    csv_map.systemic_therapies   ? file(csv_map.systemic_therapies)   : [],
                    csv_map.genomic_subjects     ? file(csv_map.genomic_subjects)     : [],
                    file(params.mohccn_treatment_intent_map)
                )
            }
            .set { ch_treatment_input }
        GENERATE_TIMELINE_TREATMENT(ch_treatment_input, ch_clinical_common)

        csvs
            .map { group, csv_map ->
                tuple(
                    [group: group],
                    csv_map.sample_registrations ? file(csv_map.sample_registrations) : [],
                    csv_map.follow_ups           ? file(csv_map.follow_ups)           : [],
                    csv_map.genomic_subjects     ? file(csv_map.genomic_subjects)     : [],
                    file(params.mohccn_primary_site_map)
                )
            }
            .set { ch_status_input }
        GENERATE_TIMELINE_STATUS(ch_status_input, ch_clinical_common)

        csvs
            .map { group, csv_map ->
                tuple(
                    [group: group],
                    csv_map.sample_registrations ? file(csv_map.sample_registrations) : [],
                    csv_map.specimens            ? file(csv_map.specimens)            : [],
                    csv_map.follow_ups           ? file(csv_map.follow_ups)           : [],
                    csv_map.genomic_subjects     ? file(csv_map.genomic_subjects)     : [],
                    file(params.mohccn_primary_site_map),
                    file(params.mohccn_specimen_tissue_source_map)
                )
            }
            .set { ch_specimen_input }
        GENERATE_TIMELINE_SPECIMEN(ch_specimen_input, ch_clinical_common)

        csvs
            .map { group, csv_map ->
                tuple(
                    [group: group],
                    csv_map.sample_registrations ? file(csv_map.sample_registrations) : [],
                    csv_map.biomarkers           ? file(csv_map.biomarkers)           : [],
                    csv_map.genomic_subjects     ? file(csv_map.genomic_subjects)     : []
                )
            }
            .set { ch_lab_test_input }
        GENERATE_TIMELINE_LAB_TEST(ch_lab_test_input, ch_clinical_common)

        ch_timeline_parts = GENERATE_TIMELINE_SURGERY.out.ch_timeline_part
            .mix(GENERATE_TIMELINE_TREATMENT.out.ch_timeline_part)
            .mix(GENERATE_TIMELINE_STATUS.out.ch_timeline_part)
            .mix(GENERATE_TIMELINE_SPECIMEN.out.ch_timeline_part)
            .mix(GENERATE_TIMELINE_LAB_TEST.out.ch_timeline_part)
            .groupTuple()

        MERGE_TIMELINE(ch_timeline_parts)

        // MERGE_TIMELINE's output is declared non-optional (see merge_timeline.R):
        // it always creates data_timeline.txt, writing a zero-byte file when none
        // of the five event types produced any rows for a group. Downstream, that
        // is exactly the "no timeline data" case the original single-process
        // GENERATE_TIMELINE represented by not emitting a file at all — so treat
        // an empty file the same way here: skip the meta file and keep it out of
        // the study package.
        ch_timeline_nonempty = MERGE_TIMELINE.out.ch_timeline
            .filter { _group, f -> f.size() > 0 }

        // Generate meta file for the combined timeline data file
        ch_timeline_nonempty
            .map { group, f -> group }
            .set { ch_timeline_groups }

        GENERATE_META_FILE_TIMELINE(
            ch_timeline_groups,
            ch_timeline_groups.map { "timeline" },
            ch_timeline_groups.map { """cancer_study_identifier: add_text
genetic_alteration_type: CLINICAL
datatype: TIMELINE
data_filename: data_timeline.txt
""" }
        )

        // ── Files that make up the clinical part of the study package ────────
        // SPLIT_CLINICAL output is deliberately excluded: per-subject slices are
        // not loadable by cBioPortal.

        ch_package_files = WRITE_CLINICAL_SAMPLE.out.ch_clinical_sample
            .mix(WRITE_CLINICAL_PATIENT.out.ch_clinical_patient)
            .mix(GENERATE_META_FILE.out)
            .mix(ch_timeline_nonempty)
            .mix(GENERATE_META_FILE_TIMELINE.out)

    emit:
        csvs
        package_files = ch_package_files
}
