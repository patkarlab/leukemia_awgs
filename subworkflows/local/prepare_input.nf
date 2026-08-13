/*
 * PREPARE_INPUT
 *
 * Reads the sample sheet and emits [meta, minknow_input].
 *
 * Columns (header row required):
 *   sample_id    Canonical sequencing ID. No PHI.
 *   minknow_bam  Either a FOLDER of per-chunk MinKNOW BAMs (real P2i output,
 *                one BAM per ~1.5 h) or a single BAM (validation sheets).
 *                MERGE_MINKNOW concatenates a folder; a single file passes
 *                through. Point a folder column at the directory that
 *                directly contains the pass BAMs, e.g. .../bam_pass.
 *   panel        Optional. AML or ALL. Defaults to --panel. This is what
 *                allows one run to cover a mixed batch: adaptive sampling
 *                bakes the panel into the flowcell, so samples sequenced on
 *                different flowcells legitimately carry different panels.
 *   timepoint    Optional. Free text; '18h' marks an early snapshot.
 *   notes        Optional. Free text, ignored by the pipeline.
 */

workflow PREPARE_INPUT {

    main:
    if (!params.sample_sheet) {
        error "sample_sheet not set; provide --sample_sheet path/to/samples.csv"
    }

    def known = params.panels?.keySet()?.collect { it.toString().toUpperCase() } ?: []

    minknow_bams = Channel
        .fromPath(params.sample_sheet, checkIfExists: true)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            if (!row.sample_id)   { error "sample sheet row missing sample_id: ${row}" }
            if (!row.minknow_bam) { error "sample sheet row missing minknow_bam for ${row.sample_id}" }

            def panel = ((row.panel ?: '').trim() ?: params.panel).toUpperCase()
            if (!(panel in known)) {
                error "Sample '${row.sample_id}' requests panel '${panel}', " +
                      "which is not defined in conf/panels.config. " +
                      "Known panels: ${known.join(', ')}"
            }

            def meta = [
                id        : row.sample_id.trim(),
                panel     : panel,
                timepoint : (row.timepoint ?: '').trim() ?: null,
                notes     : (row.notes ?: '').trim() ?: null,
            ]
            // checkIfExists so a bad path fails legibly at launch rather than
            // deep inside MERGE_MINKNOW.
            tuple(meta, file(row.minknow_bam.trim(), checkIfExists: true))
        }

    emit:
    minknow_bams = minknow_bams
}
