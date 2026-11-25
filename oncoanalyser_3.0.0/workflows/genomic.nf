// ── Imports ──
include { softwareVersionsToYAML    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { GENOMIC_CNV                  } from '../subworkflows/local/genomic_cnv'
include { GENOMIC_SV                   } from '../subworkflows/local/genomic_sv'
include { GENOMIC_EXPRESSION           } from '../subworkflows/local/genomic_expression'
include { GENOMIC_MUTATIONS            } from '../subworkflows/local/genomic_mutations'
include { GENOMIC_ML                   } from '../subworkflows/local/genomic_ml'
include { GENOMIC_AGGREGATE_OUTPUT     } from '../subworkflows/local/genomic_aggregate_output'
include { GENERATE_META_FILE           } from '../modules/local/generate_meta_file'
include { GENERATE_META_FILE as GENERATE_META_FILE_CLINICAL_SAMPLE } from '../modules/local/generate_meta_file'
include { GENERATE_CANCER_TYPE        } from '../modules/local/generate_cancer_type'
include { ISOFOX_FUSION_TO_CBIOPORTAL  } from '../modules/local/isofox_fusion_to_cbioportal'
include { MERGE_SAMPLE_SV              } from '../modules/local/merge_sample_sv'
include { SIGPROFILER_SBS               } from '../modules/local/sigprofiler_sbs'
include { SIGPROFILER_DBS              } from '../modules/local/sigprofiler_dbs'
include { SIGPROFILER_ID               } from '../modules/local/sigprofiler_id'
include { SIGS_COUNTS_TO_CBIOPORTAL    } from '../modules/local/sigs_counts_to_cbioportal'
include { GENOMIC_SAMPLE_ATTRIBUTES    } from '../modules/local/genomic_sample_attributes'
include { resource                     } from '../subworkflows/local/utils'
include { isDnaOnly                    } from '../subworkflows/local/utils'


// Resolve an oncoanalyser output file path; log a warning and return null if absent.
def findOncoFile(meta, path_str, label) {
    def f = file(path_str, checkIfExists: false)
    if (!f.exists() || f.isEmpty()) {
        log.warn "File not found for ${meta.subject} (${label}): ${path_str}"
        return null
    }
    return f
}


// False only under dna-only when <folder>/<rel> is absent: that step is skipped for the subject.
def hasOncoOutput(meta, rel) {
    return !isDnaOnly() || file("${meta.folder}/${rel}", checkIfExists: false).exists()
}


// [meta + [pipeline: p], <folder>/<rel>] per subject that has the file; `gate` (dna-only) skips subjects lacking it. {s} = subject
def oncoInput(ch_samples, rel, pipeline, label, gate = null) {
    return ch_samples
        .filter { meta -> gate == null || hasOncoOutput(meta, gate.replace('{s}', meta.subject)) }
        .map { meta ->
            def f = findOncoFile(meta, "${meta.folder}/${rel.replace('{s}', meta.subject)}", label)
            f ? [meta + [pipeline: pipeline], f] : null
        }
        .filter { it != null }
}


// dna-only inputs: label, path relative to the subject folder, steps skipped without it.
def dnaOnlyInputs(meta) {
    return [
        [label: 'purple/',            rel: 'purple',                                     steps: 'copy number, SVs, mutations and DBS/ID signatures'],
        [label: 'PURPLE SV VCF',      rel: "purple/${meta.subject}-T.purple.sv.vcf.gz",  steps: 'structural variants'],
        [label: 'sigs/',              rel: 'sigs',                                       steps: 'SBS signatures'],
    ]
}


// Subjects of the previous run, from its linking file; call before ch_linking_file overwrites it
def previousCohort(group) {
    def linking = file("${params.outdir}/${group}/util_linking_file.txt", checkIfExists: false)
    if (!linking.exists()) {
        return [] as Set
    }
    return linking.readLines()
        .drop(1)
        .findAll { line -> line?.trim() }
        .collect { line -> line.split('\t')[0] }
        .toSet()
}


// ── Main workflow ──

workflow GENOMIC {

    take:
        samplesheet_list        // from nf-schema — one record per samplesheet row
        ensembl_annotations     // path — BioMart TSV for CNV/MAF gene mapping + SV gene overlap
        vep_data                // channel<path> — pre-staged VEP cache (may be empty)
        pcgr_data               // channel<path> — pre-staged PCGR reference data (may be empty)
        mafsmith_data           // channel<path> — pre-staged mafsmith home (may be empty)
        needs_vep               // boolean — true when vep_data is not supplied
        needs_pcgr              // boolean — true when pcgr_data is not supplied
        needs_mafsmith          // boolean — true when mafsmith_data is not supplied
        fasta                   // path — GRCh38 reference FASTA (for SigProfiler indel context)
        cosmic_data             // channel<path> — COSMIC/ChimerKB fusion data for ML step
        chimer_data
        hotspots_data           // channel<path> — pre-staged cancerhotspots.org data (may be empty)

    main:

        ch_versions = channel.empty()

        // ── Parse samplesheet: one row per subject; files resolve relative to its folder ──

        ch_samples = samplesheet_list
            .map { rec ->
                [
                    group  : "${params.study_id}",
                    subject: "${rec[0].subject}",
                    sample : "${rec[0].sample}",
                    folder : "${rec[0].folder}",
                ]
            }

        // A subject dropped from the samplesheet leaves the study; read the old cohort before it is overwritten
        previous_cohort = previousCohort("${params.study_id}")

        ch_samples
            .map { meta -> meta.subject }
            .toList()
            .subscribe { subjects ->
                def dropped = (previous_cohort - (subjects as Set)).sort()
                if (dropped) {
                    log.warn "Subject(s) ${dropped.join(', ')} were in the previous run of study " +
                             "${params.study_id} but are not in this samplesheet. The samplesheet defines the " +
                             "cohort, so they will NOT appear in the merged study files, even though their " +
                             "output folders remain on disk."
                }
            }

        // ── Per-modality inputs (oncoanalyser 3.0 layout, paths relative to folder) ──

        // PURPLE somatic VCF → mutations
        ch_somatic_vcf = oncoInput(ch_samples, 'purple/{s}-T.purple.somatic.vcf.gz', 'mutation', 'mutation (PURPLE somatic)', 'purple')

        // PURPLE germline VCF → germline mutations (to CPSR unfiltered, normal column only)
        ch_germline_vcf = oncoInput(ch_samples, 'purple/{s}-T.purple.germline.vcf.gz', 'mutation_germline', 'germline mutation (PURPLE)', 'purple')

        // SAGE RNA-append VCF → somatic RNA mutations (none under dna-only)
        ch_sage_rna_vcf = isDnaOnly() ? channel.empty()
            : oncoInput(ch_samples, 'sage_append/{s}-T/{s}-T.sage.append.vcf.gz', 'mutation_rna', 'mutation (SAGE RNA append)')

        // PURPLE CNV somatic + gene TSV → copy-number; the optional purity file only sets sex (else diploid baseline)
        ch_purple_cnv = ch_samples
            .filter { meta -> hasOncoOutput(meta, 'purple') }
            .map { meta ->
                def somatic = findOncoFile(meta,
                    "${meta.folder}/purple/${meta.subject}-T.purple.cnv.somatic.tsv",
                    'cnv (PURPLE somatic)')
                def gene = findOncoFile(meta,
                    "${meta.folder}/purple/${meta.subject}-T.purple.cnv.gene.tsv",
                    'cnv (PURPLE gene)')
                def purity = file("${meta.folder}/purple/${meta.subject}-T.purple.purity.tsv", checkIfExists: false)
                if (!purity.exists() || purity.isEmpty()) {
                    log.warn "No PURPLE purity file for ${meta.subject}; chrX/chrY will use a diploid baseline"
                    purity = []
                }
                (somatic && gene) ? [meta + [pipeline: 'cnv'], somatic, gene, purity] : null
            }
            .filter { it != null }

        // PURPLE somatic SV VCF (ESVEE's records + PURPLE copy number) → structural variants
        ch_purple_sv = oncoInput(ch_samples, 'purple/{s}-T.purple.sv.vcf.gz', 'sv', 'sv (PURPLE)', 'purple/{s}-T.purple.sv.vcf.gz')

        // Isofox gene expression → TPM (none under dna-only)
        ch_isofox_exp = isDnaOnly() ? channel.empty()
            : oncoInput(ch_samples, 'isofox/{s}-T.isf.gene_data.tsv', 'expression', 'expression (Isofox)')

        // Isofox pass_fusions (tumor RNA) → RNA fusions for data_sv.txt (none under dna-only)
        ch_isofox_fusion = isDnaOnly() ? channel.empty()
            : oncoInput(ch_samples, 'isofox/{s}-T.isf.pass_fusions.tsv', 'sv_rna_fusion', 'rna fusion (Isofox)')

        // Signatures: SBS fit + trinucleotide counts from snv_counts.csv; DBS and ID extracted from the PURPLE somatic VCF
        ch_sigs_for_assignment = oncoInput(ch_samples, 'sigs/{s}-T.sig.snv_counts.csv', 'sigs', 'SNV counts (SigProfiler SBS)', 'sigs')
        ch_sigs_counts = oncoInput(ch_samples, 'sigs/{s}-T.sig.snv_counts.csv', 'sigs_counts', 'mutational signature SNV counts', 'sigs')
        ch_sigs_dbs    = oncoInput(ch_samples, 'purple/{s}-T.purple.somatic.vcf.gz', 'sigs_dbs', 'somatic VCF (SigProfiler DBS)', 'purple')
        ch_sigs_id     = oncoInput(ch_samples, 'purple/{s}-T.purple.somatic.vcf.gz', 'sigs_id', 'somatic VCF (SigProfiler ID)', 'purple')

        // PURPLE purity + CHORD → TMB/TML/HRD sample attributes. Absent file → NA.
        ch_sample_attributes = ch_samples
            .map { meta ->
                def purity = file("${meta.folder}/purple/${meta.subject}-T.purple.purity.tsv", checkIfExists: false)
                def chord  = file("${meta.folder}/chord/${meta.subject}-T.chord.prediction.tsv", checkIfExists: false)
                if (!purity.exists() || purity.isEmpty()) {
                    log.warn "No PURPLE purity file for ${meta.subject}; TMB/TML will be NA"
                    purity = []
                }
                if (!chord.exists() || chord.isEmpty()) {
                    log.warn "No CHORD prediction file for ${meta.subject}; HRD will be NA"
                    chord = []
                }
                [meta + [pipeline: 'sample_attributes'], purity, chord]
            }

        // dna-only: one summary line per oncoanalyser input absent for some subjects
        if (isDnaOnly()) {
            ch_samples
                .flatMap { meta -> dnaOnlyInputs(meta).collect { input -> [input.label, input.steps, hasOncoOutput(meta, input.rel) ? 0 : 1] } }
                .groupTuple(by: [0, 1])
                .subscribe { label, steps, missing ->
                    def n = missing.sum()
                    if (n > 0) {
                        log.warn "dna-only: ${label} absent for ${n}/${missing.size()} subject(s); ${steps} skipped for those"
                    }
                }
        }

        // ── Run subworkflows ──────────────────────────────────────────────────

        GENOMIC_CNV(ch_purple_cnv, ensembl_annotations)

        GENOMIC_SV(ch_purple_sv, ensembl_annotations)

        ISOFOX_FUSION_TO_CBIOPORTAL(ch_isofox_fusion)

        SIGPROFILER_SBS(ch_sigs_for_assignment, file(resource('sbs_signatures')), file(resource('sbs_metadata')))

        SIGPROFILER_DBS(ch_sigs_dbs, file(resource('dbs_signatures')), file(resource('dbs_metadata')))

        SIGPROFILER_ID(ch_sigs_id, file(resource('id_signatures')), file(resource('id_metadata')), fasta)

        SIGS_COUNTS_TO_CBIOPORTAL(ch_sigs_counts)

        GENOMIC_EXPRESSION(ch_isofox_exp, ensembl_annotations)

        GENOMIC_SAMPLE_ATTRIBUTES(ch_sample_attributes)

        GENOMIC_MUTATIONS(
            ch_germline_vcf,
            ch_somatic_vcf,
            ch_sage_rna_vcf,
            vep_data,
            pcgr_data,
            mafsmith_data,
            needs_vep,
            needs_pcgr,
            needs_mafsmith,
            ensembl_annotations
        )

        // ── Per-sample SV files (DNA + RNA fusion) merged; sorted files + rebuilt meta keep -resume stable ──

        ch_sv_per_sample = GENOMIC_SV.out.sv_out
            .mix(ISOFOX_FUSION_TO_CBIOPORTAL.out.sv)
            .map { meta, sv_file -> [meta.sample, meta, sv_file] }
            .groupTuple()
            .map { _sample, metas, files ->
                def meta = metas[0].subMap(['group', 'subject', 'sample', 'folder']) + [pipeline: 'sv']
                [meta, (files instanceof List ? files : [files]).sort { f -> f.name }]
            }

        MERGE_SAMPLE_SV(ch_sv_per_sample)

        // ── Aggregate per-group outputs ───────────────────────────────────────

        GENOMIC_AGGREGATE_OUTPUT(
            GENOMIC_CNV.out.segfile,
            GENOMIC_CNV.out.longfile,
            MERGE_SAMPLE_SV.out.sv,
            GENOMIC_EXPRESSION.out.out,
            GENOMIC_MUTATIONS.out.out,
            SIGPROFILER_SBS.out.sigs,
            SIGS_COUNTS_TO_CBIOPORTAL.out.sigs_counts,
            SIGPROFILER_DBS.out.sigs_dbs,
            SIGPROFILER_DBS.out.sigs_counts_dbs,
            SIGPROFILER_ID.out.sigs_id,
            SIGPROFILER_ID.out.sigs_counts_id,
            GENOMIC_SAMPLE_ATTRIBUTES.out.attributes,
        )

        // ── ML formatting (opt-in: --run_ml) ──────────────────────────────────

        if (params.run_ml) {
            GENOMIC_ML(
                GENOMIC_AGGREGATE_OUTPUT.out.cnv,
                GENOMIC_AGGREGATE_OUTPUT.out.expression,
                GENOMIC_AGGREGATE_OUTPUT.out.mutation,
                GENOMIC_AGGREGATE_OUTPUT.out.sv,
                cosmic_data,
                chimer_data,
                hotspots_data,
            )
        }

        // ── Study meta + cancer type (skipped on --incremental: cBioPortal already has them) ──

        all_groups = ch_samples.map { meta -> meta.group }.unique()

        ch_study_files = Channel.empty()
        if (!params.incremental) {
            meta_text = """type_of_cancer: ${params.cancer_type}
cancer_study_identifier: add_text
name: ${params.study_id}
description: ${params.project_description}
add_global_case_list: true
reference_genome: hg38
        """

            GENERATE_META_FILE(all_groups, "study", meta_text)
            GENERATE_CANCER_TYPE(all_groups)
            ch_study_files = GENERATE_META_FILE.out
                .mix(GENERATE_CANCER_TYPE.out.flatMap { group, ct, meta_ct ->
                    [tuple(group, ct), tuple(group, meta_ct)]
                })
        }

        // ── TMB/HRD sample attributes (cBioPortal allows only one SAMPLE_ATTRIBUTES file) ──

        ch_sample_attributes = GENOMIC_AGGREGATE_OUTPUT.out.sample_attributes

        ch_genomic_clinical_files = Channel.empty()
        if (params.mode in ['genomic', 'genomic-dna-only']) {
            ch_genomic_clinical_sample = ch_sample_attributes
                .collectFile(storeDir: "${params.outdir}") { group, f ->
                    ["${group}/data_clinical_sample.txt", f.text]
                }
                .map { f -> tuple(f.parent.name, f) }

            meta_text_clinical_sample = """cancer_study_identifier: add_text
genetic_alteration_type: CLINICAL
datatype: SAMPLE_ATTRIBUTES
data_filename: data_clinical_sample.txt
        """
            GENERATE_META_FILE_CLINICAL_SAMPLE(all_groups, "clinical_sample", meta_text_clinical_sample)
            ch_genomic_clinical_files = ch_genomic_clinical_sample.mix(GENERATE_META_FILE_CLINICAL_SAMPLE.out)
        }

        // ── Subject → tumor sample linking file ───────────────────────────────

        ch_linking_file = ch_samples
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

        // ── Study package files (packaged in main.nf so "both" mode adds the clinical files) ──

        all_package_files = GENOMIC_AGGREGATE_OUTPUT.out.cnv
            .mix(GENOMIC_AGGREGATE_OUTPUT.out.cnv_seg)
            .mix(GENOMIC_AGGREGATE_OUTPUT.out.sv)
            .mix(GENOMIC_AGGREGATE_OUTPUT.out.expression)
            .mix(GENOMIC_AGGREGATE_OUTPUT.out.mutation)
            .mix(GENOMIC_AGGREGATE_OUTPUT.out.sigs)
            .mix(GENOMIC_AGGREGATE_OUTPUT.out.sigs_counts)
            .mix(GENOMIC_AGGREGATE_OUTPUT.out.sigs_dbs)
            .mix(GENOMIC_AGGREGATE_OUTPUT.out.sigs_counts_dbs)
            .mix(GENOMIC_AGGREGATE_OUTPUT.out.sigs_id)
            .mix(GENOMIC_AGGREGATE_OUTPUT.out.sigs_counts_id)
            .mix(ch_genomic_clinical_files)
            .mix(GENOMIC_AGGREGATE_OUTPUT.out.meta_files)
            .mix(GENOMIC_AGGREGATE_OUTPUT.out.case_files)
            .mix(ch_study_files)
            .mix(ch_linking_file)

        // ── Software versions ─────────────────────────────────────────────────

        softwareVersionsToYAML(ch_versions)
            .collectFile(
                storeDir: "${params.outdir}/pipeline_info",
                name:     'software_versions.yml',
                sort:     true,
                newLine:  true,
            )

    emit:
        linking_file  = ch_linking_file     // channel<tuple(group, file)> — subject→sample linking file per group
        package_files = all_package_files   // channel<tuple(group, file)> — files to put in the study archive
        sample_attributes = ch_sample_attributes // channel<tuple(group, file)> — merged TMB/HRD, for CLINICAL to join (both mode)
}
