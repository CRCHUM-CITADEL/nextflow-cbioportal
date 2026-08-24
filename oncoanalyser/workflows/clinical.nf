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
        genomic_subjects     // val: path to genomic subjects TSV (both mode), or "" to skip filtering

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

            CLINICAL_AGGREGATE(ch_file_list, genomic_subjects)

            ch_package_files = CLINICAL_AGGREGATE.out.package_files
        } else {
            // ── Template fallback: generate minimal clinical files from sample_registrations ──
            log.info "No clinical samplesheet provided. Generating template clinical files from sample_registrations."

            // sample_registrations is a channel (from GENOMIC output) in "both" mode,
            // or a string path in "clinical" mode
            ch_linking_path = (sample_registrations instanceof CharSequence) ?
                Channel.fromPath(sample_registrations) : sample_registrations

            ch_linking = ch_linking_path
                .splitCsv(header: true, sep: ',')
                .filter { row -> row.tumour_normal_designation == "Tumour" && row.sample_type == "Total DNA" }
                .toSortedList { a, b -> a.submitter_sample_id <=> b.submitter_sample_id }
                .flatMap()
                .unique { [it.submitter_donor_id, it.submitter_specimen_id] }
                .map { row -> [subject: row.submitter_donor_id, sample: row.submitter_sample_id] }

            // Build newline-separated "subject\tsample" lines
            sample_lines = ch_linking
                .map { meta -> "${meta.subject}\t${meta.sample}" }
                .collect()
                .map { lines -> lines.join('\n') }

            // Use study_id as the group
            group = params.study_id

            GENERATE_CLINICAL_TEMPLATE(
                group,
                sample_lines
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
