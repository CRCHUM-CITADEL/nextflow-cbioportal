/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { softwareVersionsToYAML    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { GENOMIC_CNV               } from '../subworkflows/local/genomic_cnv'
include { GENOMIC_SV                } from '../subworkflows/local/genomic_sv'
include { GENOMIC_EXPRESSION        } from '../subworkflows/local/genomic_expression'
include { GENOMIC_MUTATIONS         } from '../subworkflows/local/genomic_mutations'
include { GENOMIC_ML                } from '../subworkflows/local/genomic_ml'
include { GENOMIC_AGGREGATE_OUTPUT  } from '../subworkflows/local/genomic_aggregate_output'
include { GENERATE_META_FILE        } from '../modules/local/generate_meta_file'


// Resolve an oncoanalyser output file path; log a warning and return null if absent.
def findOncoFile(meta, path_str, label) {
    def f = file(path_str, checkIfExists: false)
    if (!f.exists() || f.isEmpty()) {
        log.warn "File not found for ${meta.sample} (${label}): ${path_str}"
        return null
    }
    return f
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GENOMIC {

    take:
        samplesheet_list        // from nf-schema — one record per samplesheet row
        ensembl_annotations     // path — BioMart TSV for CNV gene mapping + ESVEE gene overlap
        ensembl_annotations_expr// path — BioMart TSV for Isofox expression Entrez ID mapping
        vep_data                // channel<path> — pre-staged VEP cache (may be empty)
        needs_vep               // boolean — true when vep_data is not supplied
        fasta                   // path — GRCh38 reference FASTA (for vcf2maf)
        cosmic_data             // channel<path> — COSMIC/ChimerKB fusion data for ML step

    main:

        ch_versions = channel.empty()

        // ── Parse samplesheet ────────────────────────────────────────────────
        // One row per subject: group, subject_id, sample_id, folder.
        // All modality files are resolved relative to folder.

        ch_samples = samplesheet_list
            .map { rec ->
                [
                    group  : "${rec[0].group}",
                    subject: "${rec[0].subject}",
                    sample : "${rec[0].sample}",
                    folder : "${rec[0].folder}",
                ]
            }

        // ── Build per-modality input channels ─────────────────────────────────
        // File naming conventions (relative to folder):
        //   sage/somatic/${sample_id}-T.sage.somatic.vcf.gz   — SAGE somatic DNA
        //   sage/germline/${sample_id}-N.sage.germline.vcf.gz — SAGE germline DNA
        //   sage/append/${sample_id}-T.sage.append.vcf.gz     — SAGE somatic RNA append
        //   esvee/caller/${sample_id}-T.esvee.unfiltered.vcf.gz
        //   purple/${sample_id}-T.purple.cnv.somatic.tsv
        //   purple/${sample_id}-T.purple.cnv.gene.tsv
        //   ${sample_id}-T.isf.fusions.csv
        //   ${sample_id}-T.isf.gene_data.csv

        // SAGE somatic VCF → mutations
        ch_sage_vcf = ch_samples
            .map { meta ->
                def vcf = findOncoFile(meta,
                    "${meta.folder}/sage/somatic/${meta.sample}-T.sage.somatic.vcf.gz",
                    'mutation (SAGE somatic)')
                vcf ? [meta + [pipeline: 'mutation'], vcf] : null
            }
            .filter { it != null }

        // SAGE germline VCF → germline mutations
        ch_sage_germline_vcf = ch_samples
            .map { meta ->
                def vcf = findOncoFile(meta,
                    "${meta.folder}/sage/germline/${meta.sample}-N.sage.germline.vcf.gz",
                    'germline mutation (SAGE)')
                vcf ? [meta + [pipeline: 'mutation_germline'], vcf] : null
            }
            .filter { it != null }

        // SAGE RNA-append VCF → somatic RNA mutations
        ch_sage_rna_vcf = ch_samples
            .map { meta ->
                def vcf = findOncoFile(meta,
                    "${meta.folder}/sage/append/${meta.sample}-T.sage.append.vcf.gz",
                    'mutation (SAGE RNA append)')
                vcf ? [meta + [pipeline: 'mutation_rna'], vcf] : null
            }
            .filter { it != null }

        // PURPLE CNV somatic + gene TSV → copy-number
        ch_purple_cnv = ch_samples
            .map { meta ->
                def somatic = findOncoFile(meta,
                    "${meta.folder}/purple/${meta.sample}-T.purple.cnv.somatic.tsv",
                    'cnv (PURPLE somatic)')
                def gene = findOncoFile(meta,
                    "${meta.folder}/purple/${meta.sample}-T.purple.cnv.gene.tsv",
                    'cnv (PURPLE gene)')
                (somatic && gene) ? [meta + [pipeline: 'cnv'], somatic, gene] : null
            }
            .filter { it != null }

        // ESVEE unfiltered VCF → structural variants
        ch_esvee_vcf = ch_samples
            .map { meta ->
                def vcf = findOncoFile(meta,
                    "${meta.folder}/esvee/caller/${meta.sample}-T.esvee.unfiltered.vcf.gz",
                    'sv (ESVEE)')
                vcf ? [meta + [pipeline: 'sv'], vcf] : null
            }
            .filter { it != null }

        // Isofox gene expression CSV → TPM
        ch_isofox_exp = ch_samples
            .map { meta ->
                def exp = findOncoFile(meta,
                    "${meta.folder}/${meta.sample}-T.isf.gene_data.csv",
                    'expression (Isofox)')
                exp ? [meta + [pipeline: 'expression'], exp] : null
            }
            .filter { it != null }

        // Isofox fusions CSV → RNA fusions
        ch_isofox_fusion = ch_samples
            .map { meta ->
                def fusion = findOncoFile(meta,
                    "${meta.folder}/${meta.sample}-T.isf.fusions.csv",
                    'sv (Isofox fusions)')
                fusion ? [meta + [pipeline: 'sv'], fusion] : null
            }
            .filter { it != null }

        // ── Run subworkflows ──────────────────────────────────────────────────

        GENOMIC_CNV(ch_purple_cnv, ensembl_annotations)

        GENOMIC_SV(ch_esvee_vcf, ch_isofox_fusion, ensembl_annotations)

        GENOMIC_EXPRESSION(ch_isofox_exp, ensembl_annotations_expr)

        GENOMIC_MUTATIONS(ch_sage_vcf, ch_sage_germline_vcf, fasta, vep_data, needs_vep)

        // ── Aggregate per-group outputs ───────────────────────────────────────

        GENOMIC_AGGREGATE_OUTPUT(
            GENOMIC_CNV.out.segfile,
            GENOMIC_CNV.out.longfile,
            GENOMIC_SV.out.sv_out,
            GENOMIC_EXPRESSION.out.out,
            GENOMIC_MUTATIONS.out.out,
        )

        // ── ML formatting ─────────────────────────────────────────────────────

        GENOMIC_ML(
            GENOMIC_AGGREGATE_OUTPUT.out.cnv,
            GENOMIC_AGGREGATE_OUTPUT.out.expression,
            GENOMIC_AGGREGATE_OUTPUT.out.mutation,
            GENOMIC_AGGREGATE_OUTPUT.out.sv,
            cosmic_data,
        )

        // ── Study-level metadata ──────────────────────────────────────────────

        all_groups = ch_samples.map { meta -> meta.group }.unique()

        meta_text = """type_of_cancer: add_text
cancer_study_identifier: add_text
name: ${params.project_name}
description: ${params.project_description}
add_global_case_list: true
reference_genome: hg38
        """

        GENERATE_META_FILE(all_groups, "study", meta_text)

        // ── Subject → tumor sample linking file ───────────────────────────────

        ch_samples
            .map    { meta -> tuple(meta.group, "${meta.subject}\t${meta.sample}") }
            .unique()
            .groupTuple()
            .map { group, lines ->
                def file_content = "subject_id\tsample_id\n" + lines.join("\n")
                def output_file  = file("${params.outdir}/${group}/util_linking_file.txt")
                output_file.parent.mkdirs()
                output_file.text = file_content
                return tuple(group, output_file)
            }

        // ── Software versions ─────────────────────────────────────────────────

        softwareVersionsToYAML(ch_versions)
            .collectFile(
                storeDir: "${params.outdir}/pipeline_info",
                name:     'software_versions.yml',
                sort:     true,
                newLine:  true,
            )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
