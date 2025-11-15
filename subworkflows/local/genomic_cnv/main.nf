include { EXTRACT_GENE_CNV_FOLD_CHANGES } from '../../../modules/local/extract_gene_cnv_fold_changes'
include { GENE_CNV_FOLD_CHANGES_TO_CBIOPORTAL } from '../../../modules/local/gene_cnv_fold_changes_to_cbioportal'
include { GENERATE_CASE_LIST } from '../../../modules/local/generate_case_list'
include { GENERATE_META_FILE } from '../../../modules/local/generate_meta_file'

workflow GENOMIC_CNV {
    take:
        cnv_vcf // tuple (meta, filepath)
        ensembl_annotations
    main:

        all_groups = cnv_vcf.map {meta, sample -> meta.group}.unique()

        cna_case_list = GENERATE_CASE_LIST(
            all_groups,
            "cnv",
            cnv_vcf.map{ meta, file -> meta.sample}.collect().map{ it.sort(false).join('\t') } // item at index 0 is samplename, join all by tabs in order to send a list
        )

        fold_change_per_gene_cnv = EXTRACT_GENE_CNV_FOLD_CHANGES(
            cnv_vcf,
            ensembl_annotations
            )

        cnv_vcf_with_fold_changes = cnv_vcf.join(fold_change_per_gene_cnv)

        cbioportal_genomic_cnv_files = GENE_CNV_FOLD_CHANGES_TO_CBIOPORTAL(
            cnv_vcf_with_fold_changes
            )

        cbioportal_genomic_cnv_seg_merged = cbioportal_genomic_cnv_files.seg
            .map {meta, file -> [meta.group, file]}
            .groupTuple()
            .flatMap { group, files ->
                files.collect { file -> [group, file]}
            }
            .collectFile(storeDir: "${params.outdir}",
                        keepHeader : true,
                        skip: 1,
                        sort: 'deep') { group, file ->
                            ["${group}/data_cna_hg38.seg", file.text]
                        }

        cbioportal_genomic_cnv_long_merged = cbioportal_genomic_cnv_files.long
            .map {meta, file -> [meta.group, file]}
            .groupTuple()
            .flatMap {group, files ->
                files.collect { file -> [group, file]}
            }
            .collectFile(storeDir : "${params.outdir}",
                        keepHeader : true,
                        skip : 1,
                        sort: 'deep') { group, file ->
                            ["${group}/data_cna_long.txt", file.text]
                        }

        meta_text_cna = """cancer_study_identifier: add_text
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

        meta_text_all = Channel.of(meta_text_cna, meta_text_long)
        file_name_all = Channel.of("cna_hg38", "cna_long")

	all_groups_times_two = all_groups.combine(meta_text_all).map {all_groups, meta_text_all -> all_groups }

        GENERATE_META_FILE(
	    all_groups_times_two,
            file_name_all,
            meta_text_all
        )


    emit:
        cna_case_list
        cbioportal_genomic_cnv_seg_merged
        cbioportal_genomic_cnv_long_merged
}
