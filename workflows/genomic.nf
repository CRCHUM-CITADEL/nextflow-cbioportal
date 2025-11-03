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
        gencode_annotations
        vep_cache
        pcgr_data
        needs_vep
        needs_pcgr
        fasta

    main:

        ch_versions = Channel.empty()

        // Create a channel where each record has: sample, filepath, germinal or somatic, pipeline label, and dna or rna
        ch_files_all = samplesheet_list
            .map { rec ->
                def sample = "${rec[0].sample}" // need to wrap it because if it's just number it will become integer and we need strings
                def file = "${params.input_dir}/${rec[0].file}"
                def type = rec[0].type
                def pipeline = rec[0].pipeline  // e.g. "cnv", "hard_filtered", etc.
                def sequence = rec[0].sequence  // e.g. "dna", "rna"
                return tuple(sample, file, type, pipeline, sequence)
            }


        // Filter out only the ones for the “cnv” pipeline
        ch_vcf_cnv = ch_files_all
            .filter { sample, file, type, pipeline, sequence ->
                pipeline == 'cnv' && type == 'somatic' && sequence == 'dna'
            }
            // then drop the pipeline field (if GENOMIC_CNV expects only sample + file)
            .map { sample, file, type, pipeline, sequence ->
                tuple(sample, file)
            }


        GENOMIC_CNV(
            ch_vcf_cnv,
            ensembl_annotations
        )

        // Filter out only the ones for the “sv” pipeline
        ch_vcf_sv = ch_files_all
            .filter { sample, file, type, pipeline, sequence ->
                pipeline == 'sv'
            }
            // then drop the pipeline field
            .map { sample, file, type, pipeline, sequence ->
                tuple(sample, file)
            }

        GENOMIC_SV(
            ch_vcf_sv
        )

        // Filter out only the ones for the “expression” pipeline
        ch_vcf_expression = ch_files_all
            .filter { sample, file, type, pipeline, sequence ->
                pipeline == 'expression'
            }
            // then drop the pipeline field
            .map { sample, file, type, pipeline, sequence ->
                tuple(sample, file)
            }

        GENOMIC_EXPRESSION(
            ch_vcf_expression,
            gencode_annotations
        )

        ch_vcf_gen_ger_dna = ch_files_all
            .filter { sample, file, type, pipeline, sequence ->
                pipeline == 'hard_filtered' && type == "germinal" && sequence == "dna"
            }
            // then drop the pipeline field
            .map { sample, file, type, pipeline, sequence ->
                tuple(sample, file)
            }

        // Filter out only the ones for the “expression” pipeline
        ch_vcf_gen_som_dna = ch_files_all
            .filter { sample, file, type, pipeline, sequence ->
                pipeline == 'hard_filtered' && type == 'somatic' && sequence == "dna"
            }
            // then drop the pipeline field
            .map { sample, file, type, pipeline, sequence ->
                tuple(sample, file)
            }

        ch_vcf_gen_som_rna = ch_files_all
            .filter { sample, file, type, pipeline, sequence ->
                pipeline == 'hard_filtered' && sequence == "rna"
            }
            // then drop the pipeline field
            .map { sample, file, type, pipeline, sequence ->
                tuple(sample, file)
            }

        GENOMIC_MUTATIONS(
            ch_vcf_gen_ger_dna,
            ch_vcf_gen_som_dna,
            ch_vcf_gen_som_rna,
            fasta,
            vep_cache,
            pcgr_data,
            needs_vep,
            needs_pcgr
        )

        meta_text = """type_of_cancer: add_text
cancer_study_identifier: add_text
name: add_text
description: add_text
add_global_case_list: true
reference_genome: hg38
        """

        GENERATE_META_FILE(
            "study",
            meta_text
        )

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
