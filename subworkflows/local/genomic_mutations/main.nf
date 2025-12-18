include { VCF2MAF } from '../../../modules/local/vcf2maf'
include { INTEGRATE_RNA_VARIANTS } from '../../../modules/local/integrate_rna_variants'
include { PCGR } from '../../../modules/local/pcgr'
include { CONVERT_CPSR_TO_MAF } from '../../../modules/local/convert_cpsr_to_maf'
include { DOWNLOAD_VEP_TEST } from '../../../modules/local/download_vep_test'
include { DOWNLOAD_PCGR } from '../../../modules/local/download_pcgr'
include { BCFTOOLS_INDEX } from '../../../modules/nf-core/bcftools/index'
include { FILTER_GERMLINE_DNA } from '../../../modules/local/filter_germline_dna'
include { GENERATE_CASE_LIST } from '../../../modules/local/generate_case_list'
include { GENERATE_META_FILE } from '../../../modules/local/generate_meta_file'
include { INDEX_FASTA } from '../../../modules/local/index_fasta'
include { SPLIT_VCF_BY_INTERVAL } from '../../../modules/local/split_vcf_by_interval'
include { GENERATE_INTERVALS } from '../../../modules/local/generate_intervals'

workflow GENOMIC_MUTATIONS {
    take:
        ger_dna_vcf // tuple (meta, filepath)
        som_dna_vcf // tuple (meta, filepath)
        som_rna_vcf // tuple (meta, filepath)
        fasta
        vep_data
        pcgr_data
        needs_vep
        needs_pcgr

    main:

        ch_vep_data = needs_vep ? DOWNLOAD_VEP_TEST().cache_dir.first() : vep_data.first()
        ch_pcgr_data = needs_pcgr ? DOWNLOAD_PCGR().data_dir.first() : pcgr_data.first()

        // index fasta
        fasta_fai = INDEX_FASTA(fasta)

	    ger_dna_filtered = FILTER_GERMLINE_DNA(ger_dna_vcf)

        ger_dna_index = BCFTOOLS_INDEX(ger_dna_filtered).tbi

        ger_dna_vcf_with_index = ger_dna_vcf
            .join(ger_dna_index)
            .map {meta, filepath, index -> tuple(meta, filepath, index)}

        // Generate intervals from reference index
        def chunk_size = 1000000  // 50000000 = 50 Mb

        intervals_ch = GENERATE_INTERVALS(fasta_fai, chunk_size)
            .flatten()
            .map { bed -> [[id: bed.baseName], bed] }
        
        // Combine VCF with each interval
        vcf_chunked = ger_dna_vcf_with_index.combine(intervals_ch)

        vcf_chunked.view()
        
        // Split VCF by intervals
        split_vcf_with_index = SPLIT_VCF_BY_INTERVAL(vcf_chunked)

        split_vcf_pcgr = split_vcf_with_index
            .map { meta, interval_id, vcf, vcf_index -> 
                def meta_interval = meta + [interval_id : interval_id ]
                return [meta_interval, vcf, vcf_index]
            }     

        ger_dna_tsv = PCGR(
            split_vcf_pcgr,
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

        som_dna_maf_tsv = som_dna_rna_maf
		.map {meta, file -> return tuple(meta.subject, meta, file)}
		.join(
			ger_dna_tsv.map {meta, file -> return tuple(meta.subject, meta, file)}
		)

        cbioportal_genomic_mutation_files = CONVERT_CPSR_TO_MAF(
            som_dna_maf_tsv
        )

    emit:
        out = cbioportal_genomic_mutation_files
}
