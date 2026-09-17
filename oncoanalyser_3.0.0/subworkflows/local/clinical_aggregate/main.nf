include { FORMAT_CLINICAL }                          from '../../../modules/local/format_clinical'
include { SPLIT_CLINICAL }                           from '../../../modules/local/split_clinical'
include { GENERATE_META_FILE }                       from '../../../modules/local/generate_meta_file'
include { GENERATE_META_FILE as GENERATE_META_FILE_TIMELINE } from '../../../modules/local/generate_meta_file'
include { GENERATE_TIMELINE }                        from '../../../modules/local/generate_timeline'

workflow CLINICAL_AGGREGATE {
    take:
        filelist          // channel: [meta(group, pipeline, extraction_date), csv_path]
        genomic_subjects  // val: path to genomic subjects TSV, or "" to skip filtering

    main:

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

        mode_ch = channel.of("sample", "patient")

        mode_ch
            .combine(csvs)
            .map { mode, group, csv_map ->
                return tuple(
                    [group: group, mode: mode],
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
            .set { ch_formatted_input }

        FORMAT_CLINICAL(ch_formatted_input)

        // ── Per-sample clinical split (both mode only) ───────────────────────
        // Pair up group-level sample + patient files
        ch_sample_file = FORMAT_CLINICAL.out
            .filter { meta, f -> meta.mode == "sample" }
            .map { meta, f -> tuple(meta.group, f) }

        ch_patient_file = FORMAT_CLINICAL.out
            .filter { meta, f -> meta.mode == "patient" }
            .map { meta, f -> tuple(meta.group, f) }

        ch_group_clinical = ch_sample_file.join(ch_patient_file)

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
        csvs
            .map { group, csv_map ->
                tuple(
                    [group: group],
                    csv_map.sample_registrations ? file(csv_map.sample_registrations) : [],
                    csv_map.treatments           ? file(csv_map.treatments)           : [],
                    csv_map.surgeries            ? file(csv_map.surgeries)            : [],
                    csv_map.systemic_therapies   ? file(csv_map.systemic_therapies)   : [],
                    csv_map.follow_ups           ? file(csv_map.follow_ups)           : [],
                    csv_map.specimens            ? file(csv_map.specimens)            : [],
                    csv_map.biomarkers           ? file(csv_map.biomarkers)           : [],
                    csv_map.genomic_subjects     ? file(csv_map.genomic_subjects)     : [],
                    file(params.mohccn_primary_site_map),
                    file(params.mohccn_specimen_tissue_source_map),
                    file(params.mohccn_treatment_intent_map)
                )
            }
            .set { ch_timeline_input }

        GENERATE_TIMELINE(ch_timeline_input)

        // Generate meta file for the combined timeline data file
        GENERATE_TIMELINE.out.ch_timeline
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

        ch_package_files = FORMAT_CLINICAL.out
            .map { meta, f -> tuple(meta.group, f) }
            .mix(GENERATE_META_FILE.out)
            .mix(GENERATE_TIMELINE.out.ch_timeline.filter { _group, f -> f })
            .mix(GENERATE_META_FILE_TIMELINE.out)

    emit:
        csvs
        package_files = ch_package_files
}
