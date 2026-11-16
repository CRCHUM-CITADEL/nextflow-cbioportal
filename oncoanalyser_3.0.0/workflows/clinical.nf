/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { CLINICAL_AGGREGATE } from '../subworkflows/local/clinical_aggregate'
include { GENERATE_META_FILE } from '../modules/local/generate_meta_file'
include { GENERATE_CLINICAL_TEMPLATE } from '../modules/local/generate_clinical_template'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CLINICAL {

    take:
        file_list            // clinical samplesheet rows (may be empty channel)
        sample_registrations // path or channel: sample_registrations.csv (template mode only)
        genomic_subjects     // val: path to genomic subjects TSV (both mode, or clinical mode with a genomic samplesheet), or "" to skip filtering
        genomic_attributes   // val: merged TMB/HRD sample attributes (both mode), or [] — joined onto data_clinical_sample.txt

    main:
        ch_versions = Channel.empty()

        if (params.clinical_samplesheet) {
            // ── Full clinical processing from clinical CSVs ──────────────────
            ch_file_list = file_list
                .map { row ->
                    def group = params.study_id
                    def sub_file = row[0].file
                    def file = sub_file.startsWith("/") ? sub_file : "${projectDir}/${sub_file}"
                    def pipeline = row[0].pipeline
                    def extraction_date = row[0].date
                    return tuple([group: group, pipeline: pipeline, extraction_date: extraction_date], file)
                }

            CLINICAL_AGGREGATE(ch_file_list, genomic_subjects, genomic_attributes)

            ch_package_files = CLINICAL_AGGREGATE.out.package_files
        } else {
            // ── Template fallback: generate minimal clinical files ──
            // Clinical mode: sample_registrations is a path to sample_registrations.csv.
            // Both mode: it is GENOMIC's util_linking_file.txt (subject_id\tsample_id), which already is the cohort.
            if (sample_registrations instanceof CharSequence) {
                log.info "No clinical samplesheet provided. Generating template clinical files from sample_registrations."

                // subject_id → genomic sample_id, or [:] when there is no genomic filter
                ch_genomic_map = genomic_subjects.map { gs ->
                    gs ? file(gs).readLines().drop(1).findAll { it.trim() }
                            .collectEntries { line -> def f = line.split('\t'); [(f[0]): f[1]] }
                       : [:]
                }

                ch_linking = Channel.fromPath(sample_registrations)
                    .splitCsv(header: true, sep: ',')
                    .filter { row -> row.sample_type == "Total DNA" &&
                        row.tumour_normal_designation == "Tumour" }
                    .toSortedList { a, b -> a.submitter_sample_id <=> b.submitter_sample_id }
                    .map { rows -> [rows.unique { [it.submitter_donor_id, it.submitter_specimen_id] }] }
                    .combine(ch_genomic_map)
                    // Same rule as build_clinical_table.R: keep genomic subjects; single-sample ones adopt the genomic sample_id.
                    .flatMap { rows, gmap ->
                        def kept   = gmap ? rows.findAll { gmap.containsKey(it.submitter_donor_id) } : rows
                        def counts = kept.countBy { it.submitter_donor_id }
                        kept.collect { row ->
                            def remap = gmap && counts[row.submitter_donor_id] == 1
                            [subject     : row.submitter_donor_id,
                             sample      : remap ? gmap[row.submitter_donor_id] : row.submitter_sample_id,
                             sample_type : "Primary"]
                        }
                    }
            } else {
                log.info "No clinical samplesheet provided. Generating template clinical files from the genomic linking file."

                // One Primary row per genomic sample, sorted for a reproducible file
                ch_linking = sample_registrations
                    .flatMap { linking ->
                        file(linking).readLines().drop(1).findAll { it.trim() }
                            .collect { line -> line.split('\t') }
                            .sort { f -> f[1] }
                            .collect { f -> [subject: f[0], sample: f[1], sample_type: "Primary"] }
                    }
            }

            // Build newline-separated "subject\tsample\tsample_type" lines
            sample_lines = ch_linking
                .map { meta -> "${meta.subject}\t${meta.sample}\t${meta.sample_type}" }
                .collect()
                .map { lines -> lines.join('\n') }

            // Use study_id as the group
            group = params.study_id

            GENERATE_CLINICAL_TEMPLATE(
                group,
                sample_lines,
                genomic_attributes
            )

            // Generate meta file for clinical sample
            meta_text = Channel.of(
                """cancer_study_identifier: add_text
genetic_alteration_type: CLINICAL
datatype: SAMPLE_ATTRIBUTES
data_filename: data_clinical_sample.txt
                """)

            file_names = Channel.of("clinical_sample")

            all_groups = channel.of(group)

            GENERATE_META_FILE(
                all_groups,
                file_names,
                meta_text
            )

            ch_package_files = GENERATE_CLINICAL_TEMPLATE.out
                .mix(GENERATE_META_FILE.out)
        }

        //
        // TASK: Aggregate software versions
        //
        softwareVersionsToYAML(ch_versions)
            .collectFile(
                storeDir: "${params.outdir}/pipeline_info",
                name: 'software_versions.yml',
                sort: true,
                newLine: true
            )

    emit:
        package_files = ch_package_files  // channel<tuple(group, file)> — files to put in the study archive
}
