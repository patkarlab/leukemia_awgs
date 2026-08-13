/*
 * T2T_TRACK
 *
 * Realign to T2T-CHM13v2.0, call SVs with three callers, merge with SURVIVOR,
 * and annotate against the acute-leukaemia fusion dictionary and anchor table.
 * KMT2A-PTD is called here rather than on hg38 because it depends on
 * supplementary alignments inside 11q23, which hg38 alt-contigs disturb.
 *
 * Input:  [meta, minknow_bam]
 */

include { REALIGN_T2T          } from '../../modules/local/realign_t2t.nf'
include { SNIFFLES             } from '../../modules/local/sniffles.nf'
include { CUTESV               } from '../../modules/local/cutesv.nf'
include { SEVERUS              } from '../../modules/local/severus.nf'
include { SURVIVOR_MERGE       } from '../../modules/local/survivor_merge.nf'
include { ANNOTATE_AL_FUSIONS  } from '../../modules/local/annotate_al_fusions.nf'
include { AUGMENT_SV_SUPPORT   } from '../../modules/local/augment_sv_support.nf'
include { MERGE_TRANSLOCATIONS } from '../../modules/local/merge_translocations.nf'
include { QC_ONTARGET          } from '../../modules/local/qc_ontarget.nf'
include { FOCAL_INTRAGENIC_DUP } from '../../modules/local/focal_intragenic_dup.nf'
include { panelPath            } from './utils.nf'
include { focalTargets         } from './utils.nf'

workflow T2T_TRACK {

    take:
    minknow_bams

    main:
    REALIGN_T2T(minknow_bams)
    t2t_bam_bai = REALIGN_T2T.out.bam_bai            // [meta, bam, bai]
    // Attach each sample's own panel BED, so an AML and an ALL sample in the
    // same run each get theirs.
    with_bed = t2t_bam_bai.map { meta, bam, bai ->
        tuple(meta, bam, bai, file(panelPath(meta, 't2t_chr'), checkIfExists: true))
    }

    if (!params.skip_qc) {
        QC_ONTARGET(with_bed)
    }

    // Intragenic duplications. Each sample is paired with every target in
    // scope for its panel, so an AML and an ALL sample in the same run scan
    // different gene sets from the same table.
    if (!params.skip_focal_dup) {
        intragenic_in = with_bed
            .combine(focalTargets('t2t', 'split_read'))
            .filter { meta, bam, bai, bed, target ->
                target.disease == 'BOTH' || target.disease == meta.panel
            }
        FOCAL_INTRAGENIC_DUP(
            intragenic_in,
            file(params.focal_dup_exon_bed_t2t ?: "${projectDir}/assets/NO_FILE")
        )
    }

    if (!params.skip_sv_calling) {
        SNIFFLES(t2t_bam_bai)
        CUTESV(t2t_bam_bai)
        SEVERUS(t2t_bam_bai)

        per_sample_for_merge = SNIFFLES.out.vcf
            .join(CUTESV.out.vcf,  by: 0)
            .join(SEVERUS.out.vcf, by: 0)
            .map { meta, sn_vcf, sn_tbi, cu_vcf, cu_tbi, sv_vcf ->
                tuple(meta, sn_vcf, cu_vcf, sv_vcf)
            }

        SURVIVOR_MERGE(per_sample_for_merge)

        if (!params.skip_fusion_annotation) {
            ANNOTATE_AL_FUSIONS(
                SURVIVOR_MERGE.out.merged_vcf.map { meta, vcf, tbi ->
                    tuple(meta, vcf, tbi, file(panelPath(meta, 't2t_chr'), checkIfExists: true))
                },
                file(params.cytoband_bed_t2t,   checkIfExists: true),
                file(params.al_fusion_dict,     checkIfExists: true),
                file(params.al_fusion_anchors,  checkIfExists: true)
            )

            // Layer real per-caller read support onto the annotated table,
            // then unite near-identical junctions. Drop the .tbi from
            // Sniffles/CuteSV; Severus already emits a plain .vcf.
            ch_augment_in = ANNOTATE_AL_FUSIONS.out.tsv
                .join(SNIFFLES.out.vcf, by: 0)
                .join(CUTESV.out.vcf,   by: 0)
                .join(SEVERUS.out.vcf,  by: 0)
                .map { meta, annotated, sn_vcf, sn_tbi, cu_vcf, cu_tbi, sv_vcf ->
                    tuple(meta, annotated, sn_vcf, cu_vcf, sv_vcf)
                }

            AUGMENT_SV_SUPPORT(ch_augment_in)
            MERGE_TRANSLOCATIONS(AUGMENT_SV_SUPPORT.out.annotated)
        }
    }

    emit:
    t2t_bam_bai   = t2t_bam_bai
    qc_coverage   = params.skip_qc         ? Channel.empty() : QC_ONTARGET.out.coverage
    qc_summary    = params.skip_qc         ? Channel.empty() : QC_ONTARGET.out.summary
    focal_dup     = params.skip_focal_dup   ? Channel.empty() : FOCAL_INTRAGENIC_DUP.out.tsv
    sniffles_vcf  = params.skip_sv_calling ? Channel.empty() : SNIFFLES.out.vcf
    cutesv_vcf    = params.skip_sv_calling ? Channel.empty() : CUTESV.out.vcf
    severus_outdir= params.skip_sv_calling ? Channel.empty() : SEVERUS.out.outdir
    merged_vcf    = params.skip_sv_calling ? Channel.empty() : SURVIVOR_MERGE.out.merged_vcf
    fusions       = (params.skip_sv_calling || params.skip_fusion_annotation) ? Channel.empty() : AUGMENT_SV_SUPPORT.out.annotated
    translocations= (params.skip_sv_calling || params.skip_fusion_annotation) ? Channel.empty() : MERGE_TRANSLOCATIONS.out.translocations
}
