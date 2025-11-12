// include modules
include { DRAGEN_FUSION_SV_TO_CBIOPORTAL } from '../../../modules/local/dragen_fusion_sv_to_cbioportal'
include { GENERATE_META_FILE } from '../../../modules/local/generate_meta_file'

workflow GENOMIC_SV {
    take:
        sv_vcf

    main:

        all_groups = sv_vcf.map {meta, sample -> meta.group}.unique()

        cbioportal_genomic_sv_files = DRAGEN_FUSION_SV_TO_CBIOPORTAL(
            sv_vcf
        )

        cbioportal_genomic_sv_merged = cbioportal_genomic_sv_files
            .map {meta, file -> [meta.group, file]}
            .groupTuple()
            .flatMap { group, files ->
                files.collect { file -> [group, file]}
            }
            .collectFile(storeDir: "${params.outdir}",
                        keepHeader : true,
                        skip: 1,
                        sort : 'deep') { group, file ->
                            ["${group}/data_sv.txt", file.text]
                        }

        meta_text = """cancer_study_identifier: add_text
genetic_alteration_type: STRUCTURAL_VARIANT
datatype: SV
stable_id: structural_variants
show_profile_in_analysis_tab: true
profile_name: Structural variants from DNA
profile_description: Structural Variant Data DNA
data_filename: data_sv.txt
        """

        meta_file = GENERATE_META_FILE(
            all_groups,
            "sv",
            meta_text
        )

    emit:
        meta_file
        cbioportal_genomic_sv_merged
}
