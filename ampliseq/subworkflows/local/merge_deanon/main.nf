include { MERGE_SV         } from '../../../modules/local/merge_sv/main'
include { MERGE_CNA        } from '../../../modules/local/merge_cna/main'
include { MERGE_MUTATIONS  } from '../../../modules/local/merge_mutations/main'
include { MERGE_SEG        } from '../../../modules/local/merge_seg/main'
include { DEANON_MUTATIONS } from '../../../modules/local/deanon_mutations/main'
include { DEANON_SV        } from '../../../modules/local/deanon_sv/main'
include { DEANON_CNA       } from '../../../modules/local/deanon_cna/main'
include { DEANON_SEG       } from '../../../modules/local/deanon_seg/main'

workflow MERGE_DEANON {

    take:
    ch_sv_files       // channel: collected per-sample sv files
    ch_cna_files      // channel: collected per-sample cna files
    ch_mutation_files // channel: collected per-sample mutation files
    ch_seg_files      // channel: collected per-sample seg files
    ch_linking        // value channel: linking file

    main:
    MERGE_SV(ch_sv_files)
    MERGE_CNA(ch_cna_files)
    MERGE_MUTATIONS(ch_mutation_files)
    MERGE_SEG(ch_seg_files)

    if (params.anonymize) {
        ch_mutations = MERGE_MUTATIONS.out
        ch_sv        = MERGE_SV.out
        ch_cna       = MERGE_CNA.out
        ch_seg       = MERGE_SEG.out
    } else {
        DEANON_MUTATIONS(MERGE_MUTATIONS.out, ch_linking)
        DEANON_SV(MERGE_SV.out, ch_linking)
        DEANON_CNA(MERGE_CNA.out, ch_linking)
        DEANON_SEG(MERGE_SEG.out, ch_linking)
        ch_mutations = DEANON_MUTATIONS.out
        ch_sv        = DEANON_SV.out
        ch_cna       = DEANON_CNA.out
        ch_seg       = DEANON_SEG.out
    }

    emit:
    mutations = ch_mutations
    sv        = ch_sv
    cna       = ch_cna
    seg       = ch_seg
}
