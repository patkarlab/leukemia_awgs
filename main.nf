#!/usr/bin/env nextflow
/*
 * leukemia_awgs
 * ----------------
 * Dual-reference adaptive WGS pipeline for acute leukaemia (AML and B/T-ALL).
 *
 * T2T track:  realign -> Sniffles2 + CuteSV + Severus -> SURVIVOR
 *                     -> AL fusion annotation; KMT2A-PTD; on-target QC
 * hg38 track: realign -> ClairS-TO (somatic SNV/indel) + ichorCNA (large CNV)
 *                     + Clair3 phased germline -> VEP -> panel filter;
 *                       FLT3-ITD
 *
 * Ported from patkarlab/mm-awgs-nextflow.
 */

nextflow.enable.dsl = 2

def printHelp() {
    log.info """
    leukemia_awgs v${workflow.manifest.version}

    USAGE
        nextflow run main.nf \\
            -profile conda,docker,gandalf \\
            --sample_sheet samples.csv \\
            --panel        AML \\
            --outdir       results/

    REQUIRED
        --sample_sheet   CSV: sample_id,minknow_bam,panel,timepoint,notes
                         (see assets/sample_sheet_template.csv)

    PANEL
        --panel          AML | ALL. Default panel for sample-sheet rows with
                         no `panel` column [${params.panel}]. A sheet may mix
                         both; each sample resolves its own BEDs.

    REFERENCES (defaults point at gandalf paths; override elsewhere)
        --t2t_fasta      T2T-CHM13v2.0 chr-named FASTA
        --t2t_mmi        Pre-built minimap2 index (optional, recommended)
        --hg38_fasta     hg38 FASTA
        --severus_vntr_bed / --severus_pon

    FOCAL DUPLICATION CALLERS
        --focal_dup_targets       TSV of ITD/PTD targets [assets/focal_duplication_targets.tsv]
        --focal_dup_hotspot_bed_hg38  Optional BED of named annotation windows
        --focal_dup_exon_bed_t2t      Optional exon BED for the PTD coverage signal

    TRACK TOGGLES
        --skip_t2t_track --skip_hg38_track
        --skip_sv_calling --skip_fusion_annotation --skip_qc
        --skip_focal_dup
        --skip_clairs_to --skip_ichorcna --skip_clair3_phased
        --skip_vep_annotate --skip_al_filter

    PROFILES
        -profile conda,docker,gandalf   recommended on the project server
        -profile test                   small synthetic inputs
        -profile stub                   stub-only mode (CI smoke test)

    OUTPUT
        --outdir         Output directory [${params.outdir}]
    """
}

include { AL_AWGS } from './workflows/al_awgs.nf'

/*
 * Help, version and input validation live inside the workflow body rather
 * than as top-level statements. Nextflow 24.10 and later reject statements
 * mixed with script declarations at file scope; this form parses on every
 * version from 23.10 onward.
 */
/*
 * Completion handlers.
 *
 * Registered from a file-scope function rather than written inline in the
 * workflow body. Inside the workflow block the name `workflow` is shadowed by
 * the block itself, so a handler closure nested there sees null instead of the
 * run metadata. At file scope the function body resolves it correctly, and
 * unlike a bare top-level `workflow.onComplete { }` statement it still parses
 * on Nextflow 24.10 and later.
 */
def registerHandlers() {
    // The metadata object is captured here, at registration time. Resolving
    // the name `workflow` from inside the handler closure returns null,
    // because by then the closure's owner chain no longer reaches the script
    // binding. WorkflowMetadata is live, so the captured reference still
    // reports the final status and duration when the handler fires.
    // params is shadowed inside the handler for the same reason as workflow,
    // so the output directory is captured here too.
    def wf  = workflow
    def out = params.outdir

    workflow.onComplete {
        log.info(
            "\n" + ("=" * 60) +
            "\nleukemia_awgs run complete\n" + ("=" * 60) +
            "\nstatus:   " + (wf.success ? "SUCCESS" : "FAILED") +
            "\nrunName:  " + wf.runName +
            "\nduration: " + wf.duration +
            "\noutdir:   " + out +
            "\nworkDir:  " + wf.workDir + "\n"
        )
    }
    workflow.onError {
        log.error "Pipeline failed: ${wf.errorMessage}"
    }
}

workflow {
    if (params.help) {
        printHelp()
    }
    else if (params.version) {
        log.info "leukemia_awgs v${workflow.manifest.version}"
    }
    else {
        if (!params.sample_sheet) {
            error "Missing required parameter: --sample_sheet. Run with --help for details."
        }
        registerHandlers()
        AL_AWGS()
    }
}
