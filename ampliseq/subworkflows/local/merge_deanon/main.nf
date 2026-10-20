include { MERGE_SV            } from '../../../modules/local/merge_sv/main'
include { MERGE_CNA           } from '../../../modules/local/merge_cna/main'
include { MERGE_MUTATIONS     } from '../../../modules/local/merge_mutations/main'
include { MERGE_SEG           } from '../../../modules/local/merge_seg/main'
include { ANNOTATE_GNOMAD     } from '../../../modules/local/annotate_gnomad/main'
include { FILTER_MUTATIONS_AF } from '../../../modules/local/filter_mutations_af/main'
include { FILTER_CNV_SAMPLES as FILTER_CNV_SAMPLES_CNA } from '../../../modules/local/filter_cnv_samples/main'
include { FILTER_CNV_SAMPLES as FILTER_CNV_SAMPLES_SEG } from '../../../modules/local/filter_cnv_samples/main'
include { DEANON_MUTATIONS    } from '../../../modules/local/deanon_mutations/main'
include { DEANON_SV           } from '../../../modules/local/deanon_sv/main'
include { DEANON_CNA          } from '../../../modules/local/deanon_cna/main'
include { DEANON_SEG          } from '../../../modules/local/deanon_seg/main'

workflow MERGE_DEANON {

    take:
    ch_sv_files       // channel: collected per-sample sv files
    ch_cna_files      // channel: collected per-sample cna files
    ch_mutation_files // channel: collected per-sample mutation files
    ch_seg_files      // channel: collected per-sample seg files
    ch_linking        // value channel: linking file
    ch_cnv_failed     // channel: list of sample IDs that failed CNV depth QC

    main:
    MERGE_SV(ch_sv_files)
    MERGE_CNA(ch_cna_files)
    MERGE_MUTATIONS(ch_mutation_files)

    // -------------------------------------------------------------------------
    // Optional filters on the merged data (before ID remapping), so they also
    // apply to samples cached from earlier incremental runs
    // -------------------------------------------------------------------------
    ch_mutations = MERGE_MUTATIONS.out
    if (params.filter_gnomad_af.toString().toBoolean()) {
        // gnomAD frequencies from Genome Nexus, fetched for every run so cached samples get them too
        ANNOTATE_GNOMAD(ch_mutations)
        ch_mutations = ANNOTATE_GNOMAD.out
    }
    if (params.filter_gnomad_af.toString().toBoolean() || params.filter_vaf.toString().toBoolean()) {
        FILTER_MUTATIONS_AF(ch_mutations)
        ch_mutations = FILTER_MUTATIONS_AF.out
    }

    ch_cna = MERGE_CNA.out
    if (params.filter_cnv_depth.toString().toBoolean()) {
        FILTER_CNV_SAMPLES_CNA(ch_cna, ch_cnv_failed, 'Sample_Id')
        ch_cna = FILTER_CNV_SAMPLES_CNA.out
    }

    DEANON_MUTATIONS(ch_mutations, ch_linking)
    DEANON_SV(MERGE_SV.out, ch_linking)
    DEANON_CNA(ch_cna, ch_linking)

    // -------------------------------------------------------------------------
    // Segments: only exported with include_seg (cBioPortal derives Fraction
    // Genome Altered from them)
    // -------------------------------------------------------------------------
    ch_seg_out = Channel.empty()
    if (params.include_seg.toString().toBoolean()) {
        MERGE_SEG(ch_seg_files)
        ch_seg = MERGE_SEG.out
        if (params.filter_cnv_depth.toString().toBoolean()) {
            FILTER_CNV_SAMPLES_SEG(ch_seg, ch_cnv_failed, 'ID')
            ch_seg = FILTER_CNV_SAMPLES_SEG.out
        }
        DEANON_SEG(ch_seg, ch_linking)
        ch_seg_out = DEANON_SEG.out
    }

    emit:
    mutations = DEANON_MUTATIONS.out
    sv        = DEANON_SV.out
    cna       = DEANON_CNA.out
    seg       = ch_seg_out
}
