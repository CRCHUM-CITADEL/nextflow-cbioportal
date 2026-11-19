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
include { resource }                                 from '../utils'

// An optional clinical CSV: its file, or [] when the group has none
def optCsv(csv_map, key) {
    return csv_map[key] ? file(csv_map[key]) : []
}

workflow CLINICAL_AGGREGATE {
    take:
        filelist          // channel: [meta(group, pipeline, extraction_date), csv_path]
        genomic_subjects  // val: path to genomic subjects TSV, or "" to skip filtering
        genomic_attributes // val: merged TMB/HRD sample attributes, or [] — joined onto data_clinical_sample.txt

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

        // ── Clinical sample/patient files: BUILD_CLINICAL_TABLE once per group; the writers only select columns ──
        ch_build_input = csvs
            .map { group, csv_map ->
                return tuple(
                    [group: group],
                    optCsv(csv_map, 'donors'),
                    optCsv(csv_map, 'primary_diagnoses'),
                    optCsv(csv_map, 'specimens'),
                    optCsv(csv_map, 'sample_registrations'),
                    optCsv(csv_map, 'treatments'),
                    optCsv(csv_map, 'surgeries'),
                    optCsv(csv_map, 'systemic_therapies'),
                    optCsv(csv_map, 'radiations'),
                    optCsv(csv_map, 'follow_ups'),
                    optCsv(csv_map, 'biomarkers'),
                    optCsv(csv_map, 'genomic_subjects'),
                    file(resource('mohccn_primary_site_map')),
                    file(resource('mohccn_specimen_tissue_source_map')),
                    file(resource('mohccn_treatment_intent_map'))
                )
            }

        BUILD_CLINICAL_TABLE(ch_build_input, ch_clinical_common)

        WRITE_CLINICAL_SAMPLE(BUILD_CLINICAL_TABLE.out.ch_merged, genomic_attributes)
        WRITE_CLINICAL_PATIENT(BUILD_CLINICAL_TABLE.out.ch_merged)

        // ── Per-sample clinical split (whenever genomic subjects are given) ──
        ch_group_clinical = WRITE_CLINICAL_SAMPLE.out.ch_clinical_sample
            .join(WRITE_CLINICAL_PATIENT.out.ch_clinical_patient)

        // Parse the genomic subjects for per-sample tuples (empty channel when unset)
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

        // ── Timeline: one optional part per EVENT_TYPE, unioned by MERGE_TIMELINE (column order: merge_timeline.R) ──
        ch_surgery_input = csvs
            .map { group, csv_map ->
                tuple(
                    [group: group],
                    optCsv(csv_map, 'sample_registrations'),
                    optCsv(csv_map, 'treatments'),
                    optCsv(csv_map, 'surgeries'),
                    optCsv(csv_map, 'genomic_subjects'),
                    file(resource('mohccn_primary_site_map')),
                    file(resource('mohccn_treatment_intent_map'))
                )
            }
        GENERATE_TIMELINE_SURGERY(ch_surgery_input, ch_clinical_common)

        ch_treatment_input = csvs
            .map { group, csv_map ->
                tuple(
                    [group: group],
                    optCsv(csv_map, 'sample_registrations'),
                    optCsv(csv_map, 'treatments'),
                    optCsv(csv_map, 'systemic_therapies'),
                    optCsv(csv_map, 'genomic_subjects'),
                    file(resource('mohccn_treatment_intent_map'))
                )
            }
        GENERATE_TIMELINE_TREATMENT(ch_treatment_input, ch_clinical_common)

        ch_status_input = csvs
            .map { group, csv_map ->
                tuple(
                    [group: group],
                    optCsv(csv_map, 'sample_registrations'),
                    optCsv(csv_map, 'follow_ups'),
                    optCsv(csv_map, 'genomic_subjects'),
                    file(resource('mohccn_primary_site_map'))
                )
            }
        GENERATE_TIMELINE_STATUS(ch_status_input, ch_clinical_common)

        ch_specimen_input = csvs
            .map { group, csv_map ->
                tuple(
                    [group: group],
                    optCsv(csv_map, 'sample_registrations'),
                    optCsv(csv_map, 'specimens'),
                    optCsv(csv_map, 'follow_ups'),
                    optCsv(csv_map, 'genomic_subjects'),
                    file(resource('mohccn_primary_site_map')),
                    file(resource('mohccn_specimen_tissue_source_map'))
                )
            }
        GENERATE_TIMELINE_SPECIMEN(ch_specimen_input, ch_clinical_common)

        ch_lab_test_input = csvs
            .map { group, csv_map ->
                tuple(
                    [group: group],
                    optCsv(csv_map, 'sample_registrations'),
                    optCsv(csv_map, 'biomarkers'),
                    optCsv(csv_map, 'genomic_subjects')
                )
            }
        GENERATE_TIMELINE_LAB_TEST(ch_lab_test_input, ch_clinical_common)

        ch_timeline_parts = GENERATE_TIMELINE_SURGERY.out.ch_timeline_part
            .mix(GENERATE_TIMELINE_TREATMENT.out.ch_timeline_part)
            .mix(GENERATE_TIMELINE_STATUS.out.ch_timeline_part)
            .mix(GENERATE_TIMELINE_SPECIMEN.out.ch_timeline_part)
            .mix(GENERATE_TIMELINE_LAB_TEST.out.ch_timeline_part)
            .groupTuple()

        MERGE_TIMELINE(ch_timeline_parts)

        // No timeline rows in a group → no data_timeline.txt or meta_timeline.txt (MERGE_TIMELINE output is optional)

        // Generate meta file for the combined timeline data file
        MERGE_TIMELINE.out.ch_timeline
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

        // ── Clinical package files (per-subject SPLIT_CLINICAL slices are not loadable, so excluded) ──

        ch_package_files = WRITE_CLINICAL_SAMPLE.out.ch_clinical_sample
            .mix(WRITE_CLINICAL_PATIENT.out.ch_clinical_patient)
            .mix(GENERATE_META_FILE.out)
            .mix(MERGE_TIMELINE.out.ch_timeline)
            .mix(GENERATE_META_FILE_TIMELINE.out)

    emit:
        csvs
        package_files = ch_package_files
}
