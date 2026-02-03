/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { GENOMIC_CNV } from '../subworkflows/local/genomic_cnv'
include { GENOMIC_SV } from '../subworkflows/local/genomic_sv'
include { GENOMIC_EXPRESSION } from '../subworkflows/local/genomic_expression'
include { GENOMIC_MUTATIONS } from '../subworkflows/local/genomic_mutations'
include { GENOMIC_AGGREGATE_OUTPUT } from '../subworkflows/local/genomic_aggregate_output'
include { GENERATE_META_FILE } from '../modules/local/generate_meta_file'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GENOMIC {

    take:
        samplesheet_list
        ensembl_annotations
        ensembl_annotations_expr
        vep_data
        pcgr_data
        needs_vep
        needs_pcgr
        fasta
        cosmic_data

    main:

        ch_versions = channel.empty()

        // Create a channel where each record has: sample, filepath, germinal or somatic, pipeline label, and dna or rna
        ch_files_all = samplesheet_list
            .map { rec ->
                def group = rec[0].group
                def subject = "${rec[0].subject}"
                def sample = "${rec[0].sample}" // need to wrap it because if it's just number it will become integer and we need strings
		// TODO: fix this...
                def sub_file = "${rec[0].file}"
                def filepath = sub_file.startsWith("/") ? sub_file : "${projectDir}/${sub_file}"
                def type = rec[0].type
                def pipeline = rec[0].pipeline  // e.g. "cnv", "hard_filtered", etc.
                def sequence = rec[0].sequence  // e.g. "dna", "rna"
                return tuple([group: group, subject : subject, sample: sample, type: type, pipeline : pipeline, sequence: sequence],filepath)
            }

       // create a channel using meta of files already ran (Channel : [meta, file]
       ch_files_ran = ch_files_all
            .filter {meta, filepath ->
                def sample_dir = file("${params.outdir}/${meta.group}/${meta.subject}")
                sample_dir.exists() && sample_dir.isDirectory()
            }
            .map { meta, filepath ->
                def baseDir = file("${params.outdir}/${meta.group}/${meta.subject}")

                if (meta.pipeline == 'cnv' ) {
                    def seg = file("${baseDir}/${meta.sample}_data_cna_hg38.seg")
                    def longfile = file("${baseDir}/${meta.sample}_data_cna_long.txt")
                    return [meta, [seg : seg, longfile: longfile]]
                }

                if (meta.pipeline == 'sv' ) {
                    def sv = file("${baseDir}/${meta.sample}.data_sv.txt")
                    return tuple(meta, sv)
                }

                if (meta.pipeline == 'expression') {
                    def tpm = file("${baseDir}/${meta.sample}.tpm.tsv")
                    return tuple(meta, tpm)
                }

                if (meta.pipeline == 'hard_filtered') {
                    def maf = file("${baseDir}/${meta.subject}.somatic_rna_germline.maf")
                    return tuple(meta, maf)
                }
           }

        // get subject names of that have not yet been run
        existing_subject_names = ch_files_ran
            .map { meta, filepath -> meta.subject }
            .collect()
            .map { it.toSet() }
            .ifEmpty([] as Set)

        ch_files_not_ran = ch_files_all
            .combine(existing_subject_names)
            .filter { meta, filepath, sample_set -> meta.subject !in sample_set }
            .map { meta, filepath, sample_set -> tuple(meta, filepath) }


        ch_files_not_ran
            .collect()
            .filter { list ->
                if (list.isEmpty()) {
                    error "According to current output directory, not samples are left to run."
                }
                return true
            }
        ch_file_to_run = ch_files_not_ran
            .branch { meta, filepath ->
                cnv             : meta.pipeline == 'cnv' && meta.type == 'somatic' && meta.sequence == 'dna'
                sv              : meta.pipeline == 'sv'
                expression      : meta.pipeline == 'expression'
                germinal_dna    : meta.pipeline == 'hard_filtered' && meta.type == 'germinal' && meta.sequence == 'dna'
                somatic_dna     : meta.pipeline == 'hard_filtered' && meta.type == 'somatic' && meta.sequence == 'dna'
                somatic_rna     : meta.pipeline == 'hard_filtered' && meta.type == 'somatic' && meta.sequence == 'rna'
            }

        GENOMIC_CNV(
            ch_file_to_run.cnv,
            ensembl_annotations
        )

        // merge new CNV results with those already ra
        all_cnv_seg_results = GENOMIC_CNV.out.segfile
            .mix(ch_files_ran
                    .filter{meta, filepath -> meta.pipeline == 'cnv'}
                    .map{meta, filepath -> [meta, filepath.seg] }
                )

        all_cnv_long_results = GENOMIC_CNV.out.longfile
            .mix(ch_files_ran
                    .filter{meta, filepath -> meta.pipeline == 'cnv'}
                    .map{meta, filepath -> [meta, filepath.longfile]}
                )

        GENOMIC_SV(
            ch_file_to_run.sv
        )

        all_sv_results = GENOMIC_SV.out
            .mix(ch_files_ran
                    .filter{meta, filepath -> meta.pipeline == 'sv'}
                )

        GENOMIC_EXPRESSION(
           ch_file_to_run.expression,
           ensembl_annotations_expr
        )

        all_expression_results = GENOMIC_EXPRESSION.out
            .mix(ch_files_ran
                    .filter{meta, filepath -> meta.pipeline == 'expression'}
                )

        GENOMIC_MUTATIONS(
            ch_file_to_run.germinal_dna,
            ch_file_to_run.somatic_dna,
            ch_file_to_run.somatic_rna,
            fasta,
            vep_data,
            pcgr_data,
            needs_vep,
            needs_pcgr
        )

        all_mutations_results = GENOMIC_MUTATIONS.out
            .mix(ch_files_ran
                    .filter{meta, filepath -> meta.pipeline == 'hard_filtered'}
                )

        // if results already existed, try to merge -----------
        GENOMIC_AGGREGATE_OUTPUT(
            all_cnv_seg_results,
            all_cnv_long_results,
            all_sv_results,
            all_expression_results,
            all_mutations_results
        )

<<<<<<< Updated upstream
        // generate meta files and linking file -------------
        all_groups = ch_files_all.map {meta, filepath -> meta.group}.unique()
=======
        
        // get output ML ready ---------------
        GENOMIC_ML(
            GENOMIC_AGGREGATE_OUTPUT.out.ml_cnv,
            GENOMIC_AGGREGATE_OUTPUT.out.ml_expression,
            GENOMIC_AGGREGATE_OUTPUT.out.ml_mutation,
            GENOMIC_AGGRAGATE_OUTPUT.out.ml_sv,
            cosmic_data,
        )


        // generate meta files and linking file -------------
        all_groups = ch_meta_all.map {meta -> meta.group}.unique()
>>>>>>> Stashed changes

        meta_text = """type_of_cancer: add_text
cancer_study_identifier: add_text
name: add_text
description: add_text
add_global_case_list: true
reference_genome: hg38
        """

        GENERATE_META_FILE(
            all_groups,
            "study",
            meta_text
         )

    ch_files_all = samplesheet_list
        .filter { rec -> rec[0].type != "germinal" && rec[0].sequence == "dna"}
        .map { rec ->
            def full_name = "${rec[0].subject}"
            def sample = "${rec[0].sample}"
            def group = rec[0].group
            return tuple(group, "${full_name}\t${sample}")
        }
        .unique()
        .groupTuple()
        .map { group, lines ->
            def file_content = "subject_id\tsample_id\n" + lines.join("\n")
            def output_file = file("${params.outdir}/${group}/util_linking_file.txt")
            output_file.parent.mkdirs()
            output_file.text = file_content
            return tuple(group, output_file)
        }

        //
        // TASK: Aggregate software versions
        //
        // TODO : add versions of software
        softwareVersionsToYAML(ch_versions)
            .collectFile(
                storeDir: "${params.outdir}/pipeline_info",
                name: 'software_versions.yml',
                sort: true,
                newLine: true,
            )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
