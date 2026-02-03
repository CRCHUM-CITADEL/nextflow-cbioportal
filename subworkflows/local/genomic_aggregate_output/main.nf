include { GENERATE_CASE_LIST } from '../../../modules/local/generate_case_list'
include { GENERATE_META_FILE } from '../../../modules/local/generate_meta_file'
include { MERGE_EXPRESSION_FILES_TO_CBIOPORTAL } from '../../../modules/local/merge_expression_files_to_cbioportal'

workflow GENOMIC_AGGREGATE_OUTPUT {

    take:
        cnv_results_seg // channel [long: [meta, seg_file], seg [meta, seg_file]]
        cnv_results_long
        sv_results
        expression_results
        mutation_results

    main:

        ch_versions = channel.empty()

        // to get all groups, just take .seg files (we assume seg and long are the same)
        all_groups = cnv_results_seg.map {meta, sample -> meta.group}.unique()

        // merge cnv ----------------------------------------
        cnv_results_seg
             .map {meta, filepath -> [meta.group, filepath]}
             .groupTuple()
             .flatMap { group, files ->
                files.collect { filepath -> [group, filepath]}
             }
             .collectFile(storeDir: "${params.outdir}",
                         keepHeader : true,
                         skip: 1,
                         sort: 'deep') { group, filepath ->
                             ["${group}/data_cna_hg38.seg", filepath.text]
                         }

         cnv_results_long
             .map {meta, filepath -> [meta.group, filepath]}
             .groupTuple()
             .flatMap {group, files ->
                 files.collect { filepath -> [group, filepath]}
             }
             .collectFile(storeDir : "${params.outdir}",
                          keepHeader : true,
                          skip : 1,
                          sort: 'deep') { group, filepath ->
                             ["${group}/data_cna_long.txt", filepath.text]
                          }

        // merge sv -------------------------------------------------------------
        sv_output = sv_results
            .map {meta, filepath -> [meta.group, filepath]}
            .groupTuple()
            .flatMap { group, files ->
                files.collect { filepath -> [group, filepath]}
            }
            .collectFile(storeDir: "${params.outdir}",
                        keepHeader : true,
                        skip: 1,
                        sort : 'deep') { group, filepath ->
                            ["${group}/data_sv.txt", filepath.text]
                        }

        // merge expression with merging not possible within pure nextflow. This will publish the file to output via the module
        tpm_file_list = expression_results
            .map { meta, filepath -> tuple(meta.group, meta, filepath) }
            .groupTuple()
            .map { group, metas, files ->
                def meta = metas[0]  // Take first meta since they share the same group
                def sortedFiles = files.sort { a, b ->
                    def na = a.toString().split(/[\/\\]/).last()
                    def nb = b.toString().split(/[\/\\]/).last()
                    na <=> nb
                }
                if (sortedFiles.size() < 2) {
                    log.warn "GENOMIC_EXPRESSION: Found ${sortedFiles.size()} TPM file(s) for group ${group}. Need at least 2 files to merge. Skipping merge step."
                    return null
                }
                return tuple(meta, sortedFiles)
            }
            .filter { it != null }

<<<<<<< Updated upstream
        MERGE_EXPRESSION_FILES_TO_CBIOPORTAL(
=======
        // TODO: add meta (or group, in the very least)
        expression_output = MERGE_EXPRESSION_FILES_TO_CBIOPORTAL(
>>>>>>> Stashed changes
            tpm_file_list
        )

        // merge mutations
        mutation_results
            .map {meta, filepath -> [meta.group, filepath]}
            .groupTuple()
            .flatMap { group, files ->
                files.collect { filepath -> [group, filepath]}
            }
            .collectFile(storeDir: "${params.outdir}",
                       keepHeader : true,
                       skip: 2,
                       sort: 'deep') { group, filepath ->
                          ["${group}/data_mutations_dna_rna_germline.txt", filepath.text]
                       }

        // create meta files and case lists ---------------------------------------------------------

        case_name_all = channel.of("cnv","sequenced")
        
        cnv_sample_list = cnv_results
            .map {meta, filepath -> meta.sample} 
            .collect()
            .map { it.sort(false).join('\t') }

        mutation_sample_list = mutation_results
            .map {meta, filepath -> meta.sample}
            .collect()
            .map { it.sort(false).join('\t') }

        case_simple_lists = cnv_sample_list.concat(mutation_sample_list)
        all_groups_cases = all_groups.combine(case_name_all).map{all_groups, case_name_all -> all_groups }

        GENERATE_CASE_LIST(
            all_groups_cases,
            case_name_all,
            case_simple_lists
        )

        // to get all groups, just take .seg files (we assume seg and long are the same)

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
"""

        meta_text_all = channel.of(meta_text_seg, meta_text_long, meta_text_sv, meta_text_expression, meta_text_mutations)
        file_name_all = channel.of("cna_hg38", "cna_long", "sv", "expression", "sequenced")
        all_groups_meta = all_groups.combine(file_name_all).map{all_groups, file_name_all -> all_groups }

        GENERATE_META_FILE(
            all_groups_meta,
            file_name_all,
            meta_text_all
        )

<<<<<<< Updated upstream

=======
    emit:
        cnv         = cnv_long_output
        expression  = expression_output
        mutation    = mutation_output
        sv          = sv_output
>>>>>>> Stashed changes
}
