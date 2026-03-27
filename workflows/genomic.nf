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
        // The oncoanalyser samplesheet has one row per file (BAM, BAI, …).
        // We deduplicate to one entry per unique (group, subject, sample,
        // sample_type, sequence_type) combination.

        ch_samples = samplesheet_list
            .map { rec ->
                [
                    group        : "${rec[0].group}",
                    subject      : "${rec[0].subject}",
                    sample       : "${rec[0].sample}",
                    sample_type  : "${rec[0].sample_type}",    // 'tumor' or 'normal'
                    sequence_type: "${rec[0].sequence_type}"   // 'dna' or 'rna'
                ]
            }
            .unique { meta -> "${meta.group}_${meta.sample}_${meta.sample_type}_${meta.sequence_type}" }

        // Tumor DNA → SAGE mutations, PURPLE CNV, ESVEE SVs
        ch_tumor_dna = ch_samples
            .filter { meta -> meta.sample_type == 'tumor' && meta.sequence_type == 'dna' }

        // Normal DNA → SAGE germline mutations
        ch_normal_dna = ch_samples
            .filter { meta -> meta.sample_type == 'normal' && meta.sequence_type == 'dna' }

        // Tumor RNA → Isofox expression + fusions
        ch_tumor_rna = ch_samples
            .filter { meta -> meta.sample_type == 'tumor' && meta.sequence_type == 'rna' }

        def onco = params.oncoanalyser_outdir

        // ── Build per-modality input channels ─────────────────────────────────

        // SAGE somatic VCF → mutations
        ch_sage_vcf = ch_tumor_dna
            .map { meta ->
                def vcf = findOncoFile(meta,
                    "${onco}/${meta.group}/sage/${meta.sample}.sage.vcf.gz",
                    'mutation (SAGE)')
                vcf ? [meta + [pipeline: 'mutation'], vcf] : null
            }
            .filter { it != null }

        // SAGE germline VCF → germline mutations
        ch_sage_germline_vcf = ch_normal_dna
            .map { meta ->
                def vcf = findOncoFile(meta,
                    "${onco}/${meta.group}/sage/germline/${meta.sample}.sage.germline.vcf.gz",
                    'germline mutation (SAGE)')
                vcf ? [meta + [pipeline: 'mutation_germline'], vcf] : null
            }
            .filter { it != null }

        // PURPLE CNV somatic + gene TSV → copy-number
        ch_purple_cnv = ch_tumor_dna
            .map { meta ->
                def somatic = findOncoFile(meta,
                    "${onco}/${meta.group}/purple/${meta.sample}.purple.cnv.somatic.tsv",
                    'cnv (PURPLE somatic)')
                def gene = findOncoFile(meta,
                    "${onco}/${meta.group}/purple/${meta.sample}.purple.cnv.gene.tsv",
                    'cnv (PURPLE gene)')
                (somatic && gene) ? [meta + [pipeline: 'cnv'], somatic, gene] : null
            }
            .filter { it != null }

        // ESVEE somatic VCF → structural variants
        ch_esvee_vcf = ch_tumor_dna
            .map { meta ->
                def vcf = findOncoFile(meta,
                    "${onco}/${meta.group}/esvee/${meta.sample}.esvee.somatic.vcf.gz",
                    'sv (ESVEE)')
                vcf ? [meta + [pipeline: 'sv'], vcf] : null
            }
            .filter { it != null }

        // Isofox expression TSV → TPM
        ch_isofox_exp = ch_tumor_rna
            .map { meta ->
                def exp = findOncoFile(meta,
                    "${onco}/${meta.group}/isofox/${meta.sample}.isofox.exp.tsv",
                    'expression (Isofox)')
                exp ? [meta + [pipeline: 'expression'], exp] : null
            }
            .filter { it != null }

        // Isofox fusion TSV → RNA fusions
        ch_isofox_fusion = ch_tumor_rna
            .map { meta ->
                def fusion = findOncoFile(meta,
                    "${onco}/${meta.group}/isofox/${meta.sample}.isofox.fusion.tsv",
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
            .filter { meta -> meta.sample_type == 'tumor' && meta.sequence_type == 'dna' }
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
