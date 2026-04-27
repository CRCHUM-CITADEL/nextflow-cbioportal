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
include { GENOMIC_ML } from '../subworkflows/local/genomic_ml'
include { GENOMIC_AGGREGATE_OUTPUT } from '../subworkflows/local/genomic_aggregate_output'
include { GENERATE_META_FILE } from '../modules/local/generate_meta_file'


// finds a file given a pattern and sends a warning to console if it's not found
// returns tuple(meta, filepath)
def findFile(meta, pattern, pipeline, quiet = false) {
    def file_path = file(pattern, checkIfExists: false)

    // If glob pattern, file() may return a list
    if (file_path instanceof List) {
        if (file_path.isEmpty()) {
            if (!quiet) log.warn "File not found for ${meta.sample} (${pipeline}): ${pattern}"
            return null
        }
        file_path = file_path[0]
    }

    // Now check if single file exists
    if (!file_path.exists() || file_path.isEmpty()) {
        if (!quiet) log.warn "File not found for ${meta.sample} (${pipeline}): ${pattern}"
        return null
    }

    def meta_with_subworkflow = meta + [pipeline: pipeline]
    return [meta_with_subworkflow, file_path]
}


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
        ch_meta_all = samplesheet_list
            .map { rec ->
                def group = "${rec[0].group}"
                def subject = "${rec[0].subject}"
                def sample = "${rec[0].sample}"
                def sub_folder = "${rec[0].folder}"
                def input_path = sub_folder.startsWith("/") ? sub_folder : "${projectDir}/${sub_folder}"
                def type = "${rec[0].type}" // e.g. germinal, somatic
                def sequence = "${rec[0].sequence}"  // e.g. "dna", "rna"
                return [group: group, subject : subject, sample: sample, type: type, sequence: sequence, input_path: input_path]
            }


        // create a channel of already-processed rows (Channel : [meta, file])
        ch_files_ran = ch_meta_all
            .filter { meta ->
                def baseDir = file("${params.outdir}/${meta.group}/${meta.subject}")
                baseDir.exists() && baseDir.isDirectory()
            }
            .flatMap { meta ->
                def baseDir = file("${params.outdir}/${meta.group}/${meta.subject}")
                def files = []

                if (meta.sequence == "dna") {
                    def cnv_seg  = findFile(meta, "${baseDir}/${meta.sample}_data_cna_hg38.seg", "cnv", true)
                    def cnv_long = findFile(meta, "${baseDir}/${meta.sample}_data_cna_long.txt", "cnv", true)
                    if (cnv_seg && cnv_long) {
                        def meta_cnv  = cnv_seg[0]
                        def seg_file  = cnv_seg[1]
                        def long_file = cnv_long[1]
                        files.add(tuple(meta_cnv, [seg: seg_file, longfile: long_file]))
                    }
                }

                if (meta.sequence == "rna") {
                    def expression = findFile(meta, "${baseDir}/${meta.sample}.tpm.tsv", "expression", true)
                    if (expression) files.add(expression)

                    def sv = findFile(meta, "${baseDir}/${meta.sample}.data_sv.txt", "sv", true)
                    if (sv) files.add(sv)
                }

                def mutation = findFile(meta, "${baseDir}/${meta.subject}.somatic_rna_germline.maf", "mutation", true)
                if (mutation) files.add(mutation)

                return files
            }

        // Filter to rows whose specific expected outputs are all present (per-row check).
        // A germline-DNA row is done when the merged MAF exists.
        // A somatic-DNA row is done when both CNV files and the merged MAF exist.
        // A somatic-RNA row is done when TPM and SV files exist.
        ch_files_not_ran = ch_meta_all
            .filter { meta ->
                def baseDir = file("${params.outdir}/${meta.group}/${meta.subject}")
                if (!baseDir.exists()) return true

                if (meta.type == 'germinal' && meta.sequence == 'dna') {
                    def maf = file("${baseDir}/${meta.subject}.somatic_rna_germline.maf")
                    return !maf.exists()
                }
                if (meta.type == 'somatic' && meta.sequence == 'dna') {
                    def seg  = file("${baseDir}/${meta.sample}_data_cna_hg38.seg")
                    def lng  = file("${baseDir}/${meta.sample}_data_cna_long.txt")
                    def maf  = file("${baseDir}/${meta.subject}.somatic_rna_germline.maf")
                    return !(seg.exists() && lng.exists() && maf.exists())
                }
                if (meta.type == 'somatic' && meta.sequence == 'rna') {
                    def tpm = file("${baseDir}/${meta.sample}.tpm.tsv")
                    def sv  = file("${baseDir}/${meta.sample}.data_sv.txt")
                    return !(tpm.exists() && sv.exists())
                }
                return true
            }

        // Log skipped rows
        ch_meta_all
            .filter { meta ->
                def baseDir = file("${params.outdir}/${meta.group}/${meta.subject}")
                if (!baseDir.exists()) return false

                if (meta.type == 'germinal' && meta.sequence == 'dna') {
                    return file("${baseDir}/${meta.subject}.somatic_rna_germline.maf").exists()
                }
                if (meta.type == 'somatic' && meta.sequence == 'dna') {
                    def seg = file("${baseDir}/${meta.sample}_data_cna_hg38.seg")
                    def lng = file("${baseDir}/${meta.sample}_data_cna_long.txt")
                    def maf = file("${baseDir}/${meta.subject}.somatic_rna_germline.maf")
                    return seg.exists() && lng.exists() && maf.exists()
                }
                if (meta.type == 'somatic' && meta.sequence == 'rna') {
                    def tpm = file("${baseDir}/${meta.sample}.tpm.tsv")
                    def sv  = file("${baseDir}/${meta.sample}.data_sv.txt")
                    return tpm.exists() && sv.exists()
                }
                return false
            }
            .subscribe { meta ->
                log.info "Skipping already-processed row: ${meta.subject} (${meta.type}/${meta.sequence})"
            }

        // Error if there are no rows left to run
        ch_files_not_ran
            .collect()
            .filter { list ->
                if (list.isEmpty()) {
                    error "According to current output directory, no samples are left to run."
                }
                return true
            }

        ch_meta_file_to_run = ch_files_not_ran
            .flatMap { meta ->

                def files = []

                if (meta.type == 'germinal') {
                    def result = findFile(meta, "${meta.input_path}/*WGS_germinal.hard-filtered.vcf.gz", "mutation")
                    if (result) files.add(result)
                }

                if (meta.type == 'somatic' && meta.sequence == 'dna') {
                    def cnv = findFile(meta, "${meta.input_path}/*.WGS_somatic-tumor_normal.cnv.vcf.gz", "cnv")
                    if (cnv) files.add(cnv)

                    def mutation = findFile(meta, "${meta.input_path}/*.WGS_somatic-tumor_normal.hard-filtered.vcf.gz", "mutation")
                    if (mutation) files.add(mutation)
                }

                if (meta.type == 'somatic' && meta.sequence == 'rna') {
                    def expression = findFile(meta, "${meta.input_path}/*.quant.genes.sf", "expression")
                    if (expression) files.add(expression)

                    def sv = findFile(meta, "${meta.input_path}/*.fusion_candidates.final", "sv")
                    if (sv) files.add(sv)

                    def mutation = findFile(meta, "${meta.input_path}/*.RNASeq_somatic.hard-filtered.vcf.gz", "mutation")
                    if (mutation) files.add(mutation)
                }

                return files
            }

        ch_file_to_run = ch_meta_file_to_run
            .branch { meta, filepath ->
                cnv             : meta.pipeline == 'cnv' && meta.type == 'somatic' && meta.sequence == 'dna'
                sv              : meta.pipeline == 'sv' && meta.type == 'somatic' && meta.sequence == 'rna'
                expression      : meta.pipeline == 'expression' && meta.type == 'somatic' && meta.sequence == 'rna'
                germinal_dna    : meta.pipeline == 'mutation' && meta.type == 'germinal' && meta.sequence == 'dna'
                somatic_dna     : meta.pipeline == 'mutation' && meta.type == 'somatic' && meta.sequence == 'dna'
                somatic_rna     : meta.pipeline == 'mutation' && meta.type == 'somatic' && meta.sequence == 'rna'
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

        all_sv_results = GENOMIC_SV.out.sv_out
            .mix(ch_files_ran
                    .filter{meta, filepath -> meta.pipeline == 'sv'}
                )

        GENOMIC_EXPRESSION(
           ch_file_to_run.expression,
           ensembl_annotations_expr
        )

        all_expression_results = GENOMIC_EXPRESSION.out.out
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

        all_mutations_results = GENOMIC_MUTATIONS.out.out
            .mix(ch_files_ran
                    .filter{meta, filepath -> meta.pipeline == 'mutation'}
                )
            .unique { it[0].subject }

        // if results already existed, try to merge -----------
        GENOMIC_AGGREGATE_OUTPUT(
            all_cnv_seg_results,
            all_cnv_long_results,
            all_sv_results,
            all_expression_results,
            all_mutations_results,
        )

        // get output ML ready ---------------
        GENOMIC_ML(
            GENOMIC_AGGREGATE_OUTPUT.out.cnv,
            GENOMIC_AGGREGATE_OUTPUT.out.expression,
            GENOMIC_AGGREGATE_OUTPUT.out.mutation,
            GENOMIC_AGGREGATE_OUTPUT.out.sv,
            cosmic_data,
        )


        // generate meta files and linking file -------------
        all_groups = ch_meta_all.map {meta -> meta.group}.unique()

        meta_text = """type_of_cancer: add_text
cancer_study_identifier: add_text
name: ${params.project_name}
description: ${params.project_description}
add_global_case_list: true
reference_genome: hg38
        """

        GENERATE_META_FILE(
            all_groups,
            "study",
            meta_text
         )

    samplesheet_list
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
