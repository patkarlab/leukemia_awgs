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
include { IGV_SNAPSHOTS } from '../modules/local/igv_snapshots.nf'
include { REPORT_BUNDLE        } from '../modules/local/report_bundle.nf'
include { DASHBOARD            } from '../modules/local/dashboard.nf'
include { EMBED_REPORT_ASSETS  } from '../modules/local/embed_report_assets.nf'
include { REPORT_ZIP           } from '../modules/local/report_zip.nf'

workflow AL_AWGS {

    // 1. Sample sheet -> [meta, minknow_input]
    PREPARE_INPUT()

    // 2. Concatenate each sample's chunk BAMs once, so both tracks share it.
    MERGE_MINKNOW(PREPARE_INPUT.out.minknow_bams)
    minknow_bams = MERGE_MINKNOW.out.merged_bam

    if (!params.skip_t2t_track)  { T2T_TRACK(minknow_bams)  }
    if (!params.skip_hg38_track) { HG38_TRACK(minknow_bams) }

    // 3. IGV pages. Needs both alignments plus the tables that name the loci,
    //    so it can only run once both tracks have produced their calls.
    if (!params.skip_igv && !params.skip_t2t_track && !params.skip_hg38_track) {
        // FOCAL_TANDEM_DUP emits one tuple per (sample, gene), so a plain
        // join takes whichever arrived first and the other gene's calls never
        // reach the IGV page. Grouped so every duplication table for a sample
        // travels together.
        ch_dups = HG38_TRACK.out.focal_dup
            .mix(T2T_TRACK.out.focal_dup)
            .groupTuple(by: 0)

        ch_igv = T2T_TRACK.out.fusions
            .join(HG38_TRACK.out.al_report, by: 0, remainder: true)
            .join(ch_dups, by: 0, remainder: true)
            .join(T2T_TRACK.out.t2t_bam_bai, by: 0)
            .join(HG38_TRACK.out.hg38_bam_bai, by: 0)
            .map { meta, fus, clin, dup, tbam, tbai, hbam, hbai ->
                tuple(meta, fus,
                      clin ?: file("${projectDir}/assets/NO_FILE"),
                      dup  ?: [file("${projectDir}/assets/NO_FILE")],
                      tbam, tbai, hbam, hbai)
            }
        IGV_SNAPSHOTS(
            ch_igv,
            file(params.t2t_fasta,  checkIfExists: true),
            file(params.t2t_fai,    checkIfExists: true),
            file(params.hg38_fasta, checkIfExists: true),
            file(params.hg38_fai,   checkIfExists: true)
        )
    }

    // 4. Report. The bundle is assembled by scanning the published output
    //    tree, so it must not start until publishDir has run for everything
    //    that feeds it. Nextflow has no "after publish" dependency, so the
    //    completion signals of the producing processes are collected and
    //    handed in as an ordering input.
    if (!params.skip_report && !params.skip_t2t_track && !params.skip_hg38_track) {
        ch_ready = T2T_TRACK.out.qc_summary.collect(flat: false).ifEmpty([])
            .mix(T2T_TRACK.out.translocations.collect(flat: false).ifEmpty([]))
            .mix(T2T_TRACK.out.gene_cnv.collect(flat: false).ifEmpty([]))
            .mix(HG38_TRACK.out.al_report.collect(flat: false).ifEmpty([]))
            .mix(HG38_TRACK.out.focal_dup.collect(flat: false).ifEmpty([]))
            .collect()

        REPORT_BUNDLE(
            ch_ready,
            file("${projectDir}/bin/build_report_bundle.sh", checkIfExists: true),
            file(params.panels[params.panel.toUpperCase()].t2t_chr, checkIfExists: true),
            file(params.panels[params.panel.toUpperCase()].hg38,    checkIfExists: true),
            params.outdir,
            params.report_bundle_name
        )
        DASHBOARD(REPORT_BUNDLE.out.bundle,
                  file("${projectDir}/bin/dashboard_builder", checkIfExists: true))
        EMBED_REPORT_ASSETS(DASHBOARD.out.bundle)
        REPORT_ZIP(EMBED_REPORT_ASSETS.out.bundle)
    }
}
