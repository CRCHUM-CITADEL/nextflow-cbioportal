include { GENERATE_CASE_LIST } from '../../../modules/local/generate_case_list'
include { GENERATE_META_FILE } from '../../../modules/local/generate_meta_file'
include { MERGE_EXPRESSION_FILES_TO_CBIOPORTAL } from '../../../modules/local/merge_expression_files_to_cbioportal'
include { MERGE_SIGS_TO_CBIOPORTAL } from '../../../modules/local/merge_sigs_to_cbioportal'
include { MERGE_SIGS_COUNTS_TO_CBIOPORTAL } from '../../../modules/local/merge_sigs_counts_to_cbioportal'
include { MERGE_SIGS_DBS_TO_CBIOPORTAL } from '../../../modules/local/merge_sigs_dbs_to_cbioportal'
include { MERGE_SIGS_COUNTS_DBS_TO_CBIOPORTAL } from '../../../modules/local/merge_sigs_counts_dbs_to_cbioportal'
include { MERGE_SIGS_ID_TO_CBIOPORTAL } from '../../../modules/local/merge_sigs_id_to_cbioportal'
include { MERGE_SIGS_COUNTS_ID_TO_CBIOPORTAL } from '../../../modules/local/merge_sigs_counts_id_to_cbioportal'

include { isDnaOnly } from '../utils'

// Per-sample [meta, file] → <outdir>/<group>/<name>, one header kept (`skip` = header lines per part)
def mergeToStudyFile(ch, name, skip) {
    return ch
        .map { meta, f -> [meta.group, f] }
        .collectFile(storeDir: "${params.outdir}", keepHeader: true, skip: skip, sort: 'deep') { group, f ->
            ["${group}/${name}", f.text]
        }
        .map { f -> tuple(f.parent.name, f) }
}

// Per-sample [meta, file] → [[group: g], files sorted by name]: arrival order never reaches the task hash
def groupFilesByName(ch) {
    return ch
        .map { meta, f -> [meta.group, f] }
        .groupTuple()
        .map { group, files -> tuple([group: group], files.sort(false) { f -> f.name }) }
}

    workflow GENOMIC_AGGREGATE_OUTPUT {

        take:
            cnv_results_seg // channel [long: [meta, seg_file], seg [meta, seg_file]]
            cnv_results_long
            sv_results
            expression_results
            mutation_results
            sigs_results
            sigs_counts_results
            sigs_dbs_results
            sigs_counts_dbs_results
            sigs_id_results
            sigs_counts_id_results
            sample_attributes_results // [meta, <sample>.genomic_sample_attributes.tsv], every subject

        main:

            ch_versions = channel.empty()

            // Groups come from the .seg files; dna-only may lack purple/, so any DNA data type counts
            all_groups = isDnaOnly()
                ? cnv_results_seg.mix(sv_results).mix(mutation_results).map { meta, _f -> meta.group }.unique()
                : cnv_results_seg.map {meta, sample -> meta.group}.unique()

            cnv_seg_output  = mergeToStudyFile(cnv_results_seg, 'data_cna_hg38.seg', 1)
            cnv_long_output = mergeToStudyFile(cnv_results_long, 'data_cna_long.txt', 1)
            sv_output       = mergeToStudyFile(sv_results, 'data_sv.txt', 1)
            mutation_output = mergeToStudyFile(mutation_results, 'data_mutations_dna_rna_germline.txt', 2)

            // Expression merge needs a module (not doable in pure Nextflow); dna-only writes none
            tpm_file_list = groupFilesByName(isDnaOnly() ? channel.empty() : expression_results)
                .filter { meta, files ->
                    if (files.size() < 2) {
                        log.warn "GENOMIC_EXPRESSION: Found ${files.size()} TPM file(s) for group ${meta.group}. Need at least 2 files to merge. Skipping merge step."
                    }
                    files.size() >= 2
                }
            expression_output = MERGE_EXPRESSION_FILES_TO_CBIOPORTAL(tpm_file_list)

            // TMB/HRD parts (4 '#' lines + header) → one intermediate; GENOMIC or CLINICAL places it (one SAMPLE_ATTRIBUTES file allowed)
            sample_attributes_output = sample_attributes_results
                .map { meta, f -> [meta.group, f] }
                .collectFile(keepHeader: true, skip: 5, sort: 'deep') { group, f ->
                    ["${group}.genomic_sample_attributes.txt", f.text]
                }
                .map { f -> tuple(f.name - '.genomic_sample_attributes.txt', f) }

            sigs_output            = MERGE_SIGS_TO_CBIOPORTAL(groupFilesByName(sigs_results))
            sigs_counts_output     = MERGE_SIGS_COUNTS_TO_CBIOPORTAL(groupFilesByName(sigs_counts_results))
            sigs_dbs_output        = MERGE_SIGS_DBS_TO_CBIOPORTAL(groupFilesByName(sigs_dbs_results))
            sigs_counts_dbs_output = MERGE_SIGS_COUNTS_DBS_TO_CBIOPORTAL(groupFilesByName(sigs_counts_dbs_results))
            sigs_id_output         = MERGE_SIGS_ID_TO_CBIOPORTAL(groupFilesByName(sigs_id_results))
            sigs_counts_id_output  = MERGE_SIGS_COUNTS_ID_TO_CBIOPORTAL(groupFilesByName(sigs_counts_id_results))

            // create meta files and case lists ---------------------------------------------------------

            // Standard cBioPortal case lists; `needs` = profiles a sample must have (`_all` comes from add_global_case_list)
            case_list_defs = [
                [suffix: 'sequenced', name: 'Samples with mutation data',
                 description: 'Samples profiled for somatic and germline mutations',
                 category: 'all_cases_with_mutation_data', needs: ['mut']],
                [suffix: 'cna', name: 'Samples with CNA data',
                 description: 'Samples profiled for copy-number alterations',
                 category: 'all_cases_with_cna_data', needs: ['cna']],
                [suffix: 'cnaseq', name: 'Samples with mutation and CNA data',
                 description: 'Samples profiled for both mutations and copy-number alterations',
                 category: 'all_cases_with_mutation_and_cna_data', needs: ['mut', 'cna']],
                [suffix: 'sv', name: 'Samples with structural variant data',
                 description: 'Samples profiled for structural variants and fusions',
                 category: 'all_cases_with_sv_data', needs: ['sv']],
                [suffix: 'rna_seq_v2_mrna', name: 'Samples with mRNA data (RNA-Seq)',
                 description: 'Samples with RNA-Seq gene expression (TPM)',
                 category: 'all_cases_with_mrna_rnaseq_data', needs: ['rna']],
                [suffix: '3way_complete', name: 'Samples with mutation, CNA and mRNA data',
                 description: 'Samples profiled for mutations, copy-number alterations and mRNA expression',
                 category: 'all_cases_with_mutation_and_cna_and_mrna_data', needs: ['mut', 'cna', 'rna']],
            ]

            // RNA: only samples in data_expression.txt (no merge under dna-only or with < 2 TPM files)
            rna_profiled = expression_output
                .map { g, _f -> g }
                .combine(expression_results.map { meta, _f -> meta.sample })
                .map { _g, sample -> ['rna', sample] }

            // toList(), not collect(): it still emits once when every input is empty
            ch_case_lists = cnv_results_seg.map { meta, _f -> ['cna', meta.sample] }
                .mix(mutation_results.map { meta, _f -> ['mut', meta.sample] })
                .mix(sv_results.map { meta, _f -> ['sv', meta.sample] })
                .mix(rna_profiled)
                .toList()
                .flatMap { pairs ->
                    def profiled = pairs.groupBy { p -> p[0] }
                        .collectEntries { k, v -> [k, v.collect { p -> p[1] } as Set] }
                    case_list_defs.collect { d ->
                        def ids = d.needs.collect { n -> profiled.get(n, [] as Set) }
                            .inject { a, b -> a.intersect(b) }
                        [d.suffix, d.name, d.description, d.category, ids.toList().sort().join('\t')]
                    }
                }
                .filter { cl -> cl[4] }  // an empty case_list_ids fails cBioPortal validation

            GENERATE_CASE_LIST(all_groups.combine(ch_case_lists))

            // Meta file templates; GENERATE_META_FILE replaces add_text with the study id
            meta_text_seg = """cancer_study_identifier: add_text
genetic_alteration_type: COPY_NUMBER_ALTERATION
datatype: SEG
reference_genome_id: hg38
description: Somatic CNA data (copy number segment file)
data_filename: data_cna_hg38.seg
            """

            meta_text_long = """cancer_study_identifier: add_text
genetic_alteration_type: COPY_NUMBER_ALTERATION
datatype: DISCRETE_LONG
stable_id: cna
show_profile_in_analysis_tab: TRUE
profile_name: Copy-number alterations
profile_description: ADD TEXT
data_filename: data_cna_long.txt
            """

            meta_text_sv = """cancer_study_identifier: add_text
genetic_alteration_type: STRUCTURAL_VARIANT
datatype: SV
stable_id: structural_variants
show_profile_in_analysis_tab: true
profile_name: Structural variants from DNA
profile_description: Structural Variant Data DNA
data_filename: data_sv.txt
namespaces: HMF
            """

            meta_text_expression = """cancer_study_identifier: add_text
genetic_alteration_type: MRNA_EXPRESSION
datatype: CONTINUOUS
stable_id: rna_seq_v2_mrna
show_profile_in_analysis_tab: true
profile_name: mRNA expression (RNA-Seq TPM)
profile_description: Expression levels (RNA-Seq TPM values)
data_filename: data_expression.txt
            """

            meta_text_mutations = """cancer_study_identifier: add_text
genetic_alteration_type: MUTATION_EXTENDED
stable_id: mutations
datatype: MAF
show_profile_in_analysis_tab: true
profile_description: ADD TEXT
profile_name: Mutations
data_filename: data_mutations_dna_rna_germline.txt
namespaces: HMF,RNA
    """

            meta_text_sigs = """cancer_study_identifier: add_text
genetic_alteration_type: GENERIC_ASSAY
generic_assay_type: MUTATIONAL_SIGNATURE
datatype: LIMIT-VALUE
stable_id: mutational_signatures_contribution_SBS
show_profile_in_analysis_tab: true
profile_name: Mutational Signatures
profile_description: Mutational signatures based on COSMIC v3.6 (WGS)
data_filename: data_mutational_signatures_contribution_SBS.txt
generic_entity_meta_properties: NAME,DESCRIPTION
pivot_threshold_value: 0.0
    """

            meta_text_counts = """cancer_study_identifier: add_text
genetic_alteration_type: GENERIC_ASSAY
generic_assay_type: MUTATIONAL_SIGNATURE
datatype: LIMIT-VALUE
stable_id: mutational_signatures_counts_SBS
show_profile_in_analysis_tab: false
profile_name: Mutational Signatures Counts
profile_description: Mutational signature trinucleotide counts based on COSMIC v3.6 (WGS)
data_filename: data_mutational_signatures_counts_SBS.txt
generic_entity_meta_properties: NAME
    """

            meta_text_sigs_dbs = """cancer_study_identifier: add_text
genetic_alteration_type: GENERIC_ASSAY
generic_assay_type: MUTATIONAL_SIGNATURE
datatype: LIMIT-VALUE
stable_id: mutational_signatures_contribution_DBS
show_profile_in_analysis_tab: true
profile_name: DBS Mutational Signatures
profile_description: Doublet base substitution signatures based on COSMIC v3.6 (WGS)
data_filename: data_mutational_signatures_contribution_DBS.txt
generic_entity_meta_properties: NAME,DESCRIPTION
pivot_threshold_value: 0.0
"""

            meta_text_counts_dbs = """cancer_study_identifier: add_text
genetic_alteration_type: GENERIC_ASSAY
generic_assay_type: MUTATIONAL_SIGNATURE
datatype: LIMIT-VALUE
stable_id: mutational_signatures_counts_DBS
show_profile_in_analysis_tab: false
profile_name: DBS Mutational Signatures Counts
profile_description: Doublet base substitution counts based on COSMIC v3.6 (WGS)
data_filename: data_mutational_signatures_counts_DBS.txt
generic_entity_meta_properties: NAME
    """

            meta_text_sigs_id = """cancer_study_identifier: add_text
genetic_alteration_type: GENERIC_ASSAY
generic_assay_type: MUTATIONAL_SIGNATURE
datatype: LIMIT-VALUE
stable_id: mutational_signatures_contribution_ID
show_profile_in_analysis_tab: true
profile_name: ID Mutational Signatures
profile_description: Indel signatures based on COSMIC v3.6 (WGS)
data_filename: data_mutational_signatures_contribution_ID.txt
generic_entity_meta_properties: NAME,DESCRIPTION
pivot_threshold_value: 0.0
    """

            meta_text_counts_id = """cancer_study_identifier: add_text
genetic_alteration_type: GENERIC_ASSAY
generic_assay_type: MUTATIONAL_SIGNATURE
datatype: LIMIT-VALUE
stable_id: mutational_signatures_counts_ID
show_profile_in_analysis_tab: false
profile_name: ID Mutational Signatures Counts
profile_description: Indel classification counts based on COSMIC v3.6 (WGS)
data_filename: data_mutational_signatures_counts_ID.txt
generic_entity_meta_properties: NAME
    """

        // No expression meta file under dna-only: it would point at a data file that is never written
        meta_entries = [
            ["cna_hg38", meta_text_seg],
            ["cna_long", meta_text_long],
            ["sv", meta_text_sv],
        ]
        if (!isDnaOnly()) {
            meta_entries << ["expression", meta_text_expression]
        }
        meta_entries += [
            ["sequenced", meta_text_mutations],
            ["mutational_signatures_contribution_SBS", meta_text_sigs],
            ["mutational_signatures_counts_SBS", meta_text_counts],
            ["mutational_signatures_contribution_DBS", meta_text_sigs_dbs],
            ["mutational_signatures_counts_DBS", meta_text_counts_dbs],
            ["mutational_signatures_contribution_ID", meta_text_sigs_id],
            ["mutational_signatures_counts_ID", meta_text_counts_id],
        ]

        file_name_all = channel.fromList(meta_entries.collect { e -> e[0] })
        meta_text_all  = channel.fromList(meta_entries.collect { e -> e[1] })
        all_groups_meta = all_groups.combine(file_name_all).map { g, _f -> g }

        // dna-only: a meta file only for data files produced (e.g. no sigs/ on WES → no SBS meta)
        if (isDnaOnly()) {
            ch_present = cnv_seg_output.map { g, _f -> [g, 'cna_hg38'] }
                .mix(cnv_long_output.map        { g, _f -> [g, 'cna_long'] })
                .mix(sv_output.map              { g, _f -> [g, 'sv'] })
                .mix(mutation_output.map        { g, _f -> [g, 'sequenced'] })
                .mix(sigs_output.map            { g, _f -> [g, 'mutational_signatures_contribution_SBS'] })
                .mix(sigs_counts_output.map     { g, _f -> [g, 'mutational_signatures_counts_SBS'] })
                .mix(sigs_dbs_output.map        { g, _f -> [g, 'mutational_signatures_contribution_DBS'] })
                .mix(sigs_counts_dbs_output.map { g, _f -> [g, 'mutational_signatures_counts_DBS'] })
                .mix(sigs_id_output.map         { g, _f -> [g, 'mutational_signatures_contribution_ID'] })
                .mix(sigs_counts_id_output.map  { g, _f -> [g, 'mutational_signatures_counts_ID'] })
                .map { g, name -> "${g}/${name}".toString() }
                .collect()
                .map { names -> names.toSet() }
                .ifEmpty([] as Set)

            ch_meta_in = all_groups
                .combine(channel.fromList(meta_entries))
                .combine(ch_present)
                .filter { g, name, _text, present -> "${g}/${name}".toString() in present }
                .multiMap { g, name, text, _present -> group: g; label: name; text: text }
            meta_groups = ch_meta_in.group
            meta_labels = ch_meta_in.label
            meta_texts  = ch_meta_in.text
        } else {
            meta_groups = all_groups_meta
            meta_labels = file_name_all
            meta_texts  = meta_text_all
        }

        GENERATE_META_FILE(
            meta_groups,
            meta_labels,
            meta_texts
        )

    emit:
        cnv         = cnv_long_output
        cnv_seg     = cnv_seg_output
        expression  = expression_output
        mutation    = mutation_output
        sv          = sv_output
        sigs            = sigs_output
        sigs_counts     = sigs_counts_output
        sigs_dbs        = sigs_dbs_output
        sigs_counts_dbs = sigs_counts_dbs_output
        sigs_id         = sigs_id_output
        sigs_counts_id  = sigs_counts_id_output
        sample_attributes = sample_attributes_output
        meta_files      = GENERATE_META_FILE.out
        case_files  = GENERATE_CASE_LIST.out
}
