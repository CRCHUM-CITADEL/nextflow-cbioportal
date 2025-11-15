include { VCF2MAF } from '../../../modules/nf-core/vcf2maf'
include { INTEGRATE_RNA_VARIANTS } from '../../../modules/local/integrate_rna_variants'
include { PCGR } from '../../../modules/local/pcgr'
include { CONVERT_CPSR_TO_MAF } from '../../../modules/local/convert_cpsr_to_maf'
include { DOWNLOAD_VEP_TEST } from '../../../modules/local/download_vep_test'
include { DOWNLOAD_PCGR } from '../../../modules/local/download_pcgr'
include { BCFTOOLS_INDEX } from '../../../modules/nf-core/bcftools/index'
include { GENERATE_CASE_LIST } from '../../../modules/local/generate_case_list'
include { GENERATE_META_FILE } from '../../../modules/local/generate_meta_file'

workflow GENOMIC_MUTATIONS {
    take:
        ger_dna_vcf // tuple (meta, filepath)
        som_dna_vcf // tuple (meta, filepath)
        som_rna_vcf // tuple (meta, filepath)
        fasta
        vep_cache
        pcgr_data
        needs_vep
        needs_pcgr

    main:


        ch_vep_data = needs_vep ? DOWNLOAD_VEP_TEST().cache_dir.first() : vep_cache.first()
        ch_pcgr_data = needs_pcgr ? DOWNLOAD_PCGR().data_dir.first() : pcgr_data.first()

        ger_dna_index = BCFTOOLS_INDEX(ger_dna_vcf).tbi

        ger_dna_vcf_with_index = ger_dna_vcf
            .join(ger_dna_index)
            .map {meta, file, index -> tuple(meta, file, index)}

        ger_dna_tsv = PCGR(
            ger_dna_vcf_with_index,
            ch_vep_data,
            ch_pcgr_data
        )

        // in order to get meta.tumor_sample and meta.normal_sample,
        // we need to join dna's on the same subject name.
        som_dna_vcf_input = som_dna_vcf
            .map { meta, vcf -> tuple(meta.subject, meta, vcf) }
            .join(
                ger_dna_vcf.map { meta, vcf -> tuple(meta.subject, meta) }
            )
            .map { subject, som_meta, som_vcf, ger_meta ->
                def meta = [*:som_meta, germinal_sample: ger_meta.sample]
                return tuple(meta, som_vcf)
            }

        VCF2MAF(
            som_dna_vcf_input,
            fasta,
            ch_vep_data
        )

        som_dna_maf = VCF2MAF.out.maf.map { meta, vcf ->
            return tuple(meta, vcf)
        }

        // join on ID to create tuple(subject, dna, rna)
        som_rna_dna_tuple = som_rna_vcf
           .map { meta, file -> tuple(meta.subject, meta, file) }
           .join(
                som_dna_maf.map { meta, file -> tuple(meta.subject, meta, file) }
            )

        som_dna_rna_maf = INTEGRATE_RNA_VARIANTS(
            som_rna_dna_tuple
        )

        som_dna_maf_tsv = som_dna_rna_maf.join(ger_dna_tsv) 

        cbioportal_genomic_mutation_files = CONVERT_CPSR_TO_MAF(
            som_dna_maf_tsv
        )

        cbioportal_genomic_mutations_merged = cbioportal_genomic_mutation_files
            .map {meta, file -> [meta.group, file]}
           .groupTuple()
            .flatMap { group, files ->
                files.collect { file -> [group, file]}
            }
            .collectFile(storeDir: "${params.outdir}",
                       keepHeader : true,
                       skip: 2,
                        sort: 'deep') { group, file ->
                            ["${group}/data_mutations_dna_rna_germline.txt", file.text]
                        }

        all_groups = cbioportal_genomic_mutation_files.map {meta, sample -> meta.group}.unique()

        sequenced_case_list = GENERATE_CASE_LIST(
            all_groups,
            "sequenced",
            som_dna_vcf.map { meta, file -> meta.sample}.collect().map{ it.sort(false).join('\t') } // item at index 1 is sample_id, join by tabs in order to send a list
        )

        meta_text = """cancer_study_identifier: add_text
genetic_alteration_type: MUTATION_EXTENDED
stable_id: mutations
datatype: MAF
show_profile_in_analysis_tab: true
profile_description: ADD TEXT
profile_name: Mutations
data_filename: data_mutations_dna_rna_germline.txt
"""

        meta_file = GENERATE_META_FILE(
            all_groups,
            "mutations",
            meta_text
        )

    emit:
        meta_file
        sequenced_case_list
        cbioportal_genomic_mutations_merged

}
