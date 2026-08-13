/*
 * AL_AWGS top workflow.
 *
 *   MinKNOW BAM ─┬─► T2T track  : SV ensemble, fusion annotation, KMT2A-PTD, QC
 *                └─► hg38 track : SNV/indel, CNV, phased germline, FLT3-ITD
 *
 * Each track owns its realignment, so the two run concurrently.
 * Panel selection is per sample via meta.panel; a mixed AML/ALL batch is a
 * single run.
 */

include { PREPARE_INPUT } from '../subworkflows/local/prepare_input.nf'
include { MERGE_MINKNOW } from '../modules/local/merge_minknow.nf'
include { T2T_TRACK     } from '../subworkflows/local/t2t_track.nf'
include { HG38_TRACK    } from '../subworkflows/local/hg38_track.nf'

workflow AL_AWGS {

    // 1. Sample sheet -> [meta, minknow_input]
    PREPARE_INPUT()

    // 2. Concatenate each sample's chunk BAMs once, so both tracks share it.
    MERGE_MINKNOW(PREPARE_INPUT.out.minknow_bams)
    minknow_bams = MERGE_MINKNOW.out.merged_bam

    if (!params.skip_t2t_track)  { T2T_TRACK(minknow_bams)  }
    if (!params.skip_hg38_track) { HG38_TRACK(minknow_bams) }
}
