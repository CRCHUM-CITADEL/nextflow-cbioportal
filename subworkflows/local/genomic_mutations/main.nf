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
        ger_dna_vcf // tuple (sample_id, filepath)
        som_dna_vcf // tuple (sample_id, filepath)
        som_rna_vcf // tuple (sample_id, filepath)
        fasta
        vep_cache
        pcgr_data
        needs_vep
        needs_pcgr

    main:

        ch_vep_data = needs_vep ? DOWNLOAD_VEP_TEST().cache_dir.first() : vep_cache.first()
        ch_pcgr_data = needs_pcgr ? DOWNLOAD_PCGR().data_dir.first() : pcgr_data.first()


        vcf_index = BCFTOOLS_INDEX(ger_dna_vcf)

        // Join them
        ger_dna_vcf_with_index = ger_dna_vcf.join(vcf_index.tbi)

        ger_dna_tsv = PCGR(
            ger_dna_vcf_with_index,
            ch_vep_data,
            ch_pcgr_data
        )

        som_dna_vcf_input = som_dna_vcf.map { id, vcf ->
            def meta = [ id: id ]
            return tuple(meta, vcf)
        }

        VCF2MAF(
            som_dna_vcf_input,
            fasta,
            ch_vep_data
        )

        som_dna_maf = VCF2MAF.out.maf.map { meta, vcf ->
            return tuple(meta.id, vcf)
        }

        // join on ID to create tuple(id, dna, rna)
        som_rna_dna_tuple = som_rna_vcf
            .join(som_dna_maf)

        som_dna_rna_maf = INTEGRATE_RNA_VARIANTS(
            som_rna_dna_tuple
        )

        cbioportal_genomic_mutation_files = CONVERT_CPSR_TO_MAF (
            som_dna_rna_maf,
            ger_dna_tsv
        )

        sequenced_case_list = GENERATE_CASE_LIST(
            Channel.of("sequenced"),
            som_dna_rna_maf.map { it[0]}.collect().map{it.sort(false).join('\t') } // item at index 0 is sample_id, join by tabs in order to send a list
        )

        cbioportal_genomic_mutations_merged = cbioportal_genomic_mutation_files
            .collectFile( name : 'data_mutations_dna_rna_germline.txt', storeDir: "${params.outdir}", keepHeader : true, skip: 2, sort: 'deep')
        
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
            "mutations",
            meta_text
        )

    emit:
        meta_file
        sequenced_case_list
        cbioportal_genomic_mutations_merged

}
