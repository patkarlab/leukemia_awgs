/*
 * HG38_TRACK
 *
 * Realign to hg38, then in parallel:
 *   ClairS-TO   somatic SNV/indel, restricted to the panel BED
 *   ichorCNA    large CNV from off-target reads
 *   Clair3      phased germline, feeding VEP and the panel filter
 *   FLT3-ITD    dedicated caller; general-purpose callers under-call these
 *
 * Input:  [meta, minknow_bam]
 */

include { REALIGN_HG38        } from '../../modules/local/realign_hg38.nf'
include { CLAIRS_TO           } from '../../modules/local/clairs_to.nf'
include { ICHORCNA            } from '../../modules/local/ichorcna.nf'
include { CLAIR3_PHASED       } from '../../modules/local/clair3_phased.nf'
include { VEP_ANNOTATE_CLAIR3 } from '../../modules/local/vep_annotate_clair3.nf'
include { FILTER_AL_REPORT    } from '../../modules/local/filter_al_report.nf'
include { FOCAL_TANDEM_DUP    } from '../../modules/local/focal_tandem_dup.nf'
include { panelPath           } from './utils.nf'
include { focalTargets        } from './utils.nf'

workflow HG38_TRACK {

    take:
    minknow_bams

    main:
    REALIGN_HG38(minknow_bams)
    hg38_bam_bai = REALIGN_HG38.out.bam_bai
    with_bed = hg38_bam_bai.map { meta, bam, bai ->
        tuple(meta, bam, bai, file(panelPath(meta, 'hg38'), checkIfExists: true))
    }

    if (!params.skip_clairs_to) { CLAIRS_TO(with_bed) }
    if (!params.skip_ichorcna)  { ICHORCNA(with_bed)  }

    // Short internal tandem duplications. Targets come from the table; the
    // gene interval comes from each sample's own panel BED.
    if (!params.skip_focal_dup) {
        tandem_in = with_bed
            .combine(focalTargets('hg38', 'insertion'))
            .filter { meta, bam, bai, bed, target ->
                target.disease == 'BOTH' || target.disease == meta.panel
            }
        FOCAL_TANDEM_DUP(
            tandem_in,
            file(params.hg38_fasta, checkIfExists: true),
            file(params.hg38_fai,   checkIfExists: true),
            file(params.focal_dup_hotspot_bed_hg38 ?: "${projectDir}/assets/NO_FILE")
        )
    }

    if (!params.skip_clair3_phased) {
        CLAIR3_PHASED(with_bed)
        if (!params.skip_vep_annotate) {
            VEP_ANNOTATE_CLAIR3(CLAIR3_PHASED.out.merge_output)
            if (!params.skip_al_filter) {
                FILTER_AL_REPORT(
                    VEP_ANNOTATE_CLAIR3.out.candidates_tsv.map { meta, tsv ->
                        tuple(meta, tsv, file(panelPath(meta, 'hg38'), checkIfExists: true))
                    }
                )
            }
        }
    }

    emit:
    hg38_bam_bai         = hg38_bam_bai
    clairs_to_outdir     = params.skip_clairs_to     ? Channel.empty() : CLAIRS_TO.out.outdir
    ichorcna_outdir      = params.skip_ichorcna      ? Channel.empty() : ICHORCNA.out.outdir
    focal_dup            = params.skip_focal_dup     ? Channel.empty() : FOCAL_TANDEM_DUP.out.tsv
    clair3_phased_outdir = params.skip_clair3_phased ? Channel.empty() : CLAIR3_PHASED.out.outdir
    al_report            = (params.skip_clair3_phased || params.skip_vep_annotate || params.skip_al_filter) ? Channel.empty() : FILTER_AL_REPORT.out.clinical
}
