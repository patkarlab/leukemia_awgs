process IGV_SNAPSHOTS {
    tag      "${meta.id}"
    label    'process_medium'

    // Published directly rather than via modules.config, because
    // REPORT_BUNDLE scans the published tree: if these pages are not on disk
    // under outdir before the bundle runs, the report has no IGV to link to.
    publishDir path: { "${params.outdir}/igv" }, mode: 'copy', overwrite: true

    // Renders self-contained igv-reports pages for one sample, for both
    // evidence classes:
    //
    //   translocations  two standalone pages per rearrangement, one per
    //                   breakpoint, against the T2T BAM, plus a manifest
    //                   linking them. The dashboard loads the pair side by
    //                   side, which is why they are separate documents.
    //
    //   somatic         one page covering the clinical SNV table, against the
    //                   hg38 BAM. It is published as <sample>_igv_report.html
    //                   because that is the exact filename the dashboard
    //                   builder resolves, and it also feeds the per-variant
    //                   IGV links via the row lookup the builder extracts.
    //
    // Both classes are optional per sample. A sample with no clinical SNVs, or
    // none of the selected SV types, produces a placeholder page and exits 0;
    // that is an expected outcome and must not fail the run.

    input:
    // meta.id is the sequencing identifier only. No patient identifier ever
    // reaches this process.
    tuple val(meta), path(leukemia_annotated), path(clinical_tsv), path(dup_tsvs), path(t2t_bam), path(t2t_bai), path(hg38_bam), path(hg38_bai)
    path(t2t_fasta,  stageAs: 't2t_ref/*')
    path(t2t_fai,    stageAs: 't2t_ref/*')
    path(hg38_fasta, stageAs: 'hg38_ref/*')
    path(hg38_fai,   stageAs: 'hg38_ref/*')

    output:
    tuple val(meta), path("${meta.id}"),                        emit: igv
    tuple val(meta), path("${meta.id}/somatic/*.html"),         emit: somatic,  optional: true
    tuple val(meta), path("${meta.id}/focal_dup/*.html"),        emit: focaldup, optional: true
    tuple val(meta), path("${meta.id}/translocations/*.json"),  emit: manifest, optional: true
    path "versions.yml",                                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def flanking     = params.igv_flanking     ?: 5000
    def sv_types     = params.igv_sv_types     ?: 'TRA'
    def min_callers  = params.igv_min_callers  ?: 1
    def max_events   = params.igv_max_events   ?: 200
    """
    set -euo pipefail
    mkdir -p ${meta.id}/translocations ${meta.id}/somatic ${meta.id}/focal_dup

    # Translocations. The annotated table carries the whole merged callset,
    # so selection by sv_type happens inside the script and is reported to
    # stderr before any page is rendered.
    if [ -s "${leukemia_annotated}" ]; then
        igv_snapshots.py \\
            --mode translocations \\
            --sample ${meta.id} \\
            --sites-tsv ${leukemia_annotated} \\
            --bam ${t2t_bam} \\
            --fasta ${t2t_fasta} \\
            --out-html ${meta.id}/translocations/${meta.id}.translocations.html \\
            --out-dir  ${meta.id}/translocations \\
            --flanking ${flanking} \\
            --sv-types '${sv_types}' \\
            --min-callers ${min_callers} \\
            --max-events ${max_events}
    else
        echo "No annotated SV table for ${meta.id}; skipping translocation snapshots" >&2
    fi

    # Clinical SNVs. An empty clinical table is a legitimate result for a
    # sample with no on-panel variants, so this is guarded rather than
    # assumed. The output filename is the one the dashboard builder resolves.
    if [ -n "${clinical_tsv}" ] && [ -s "${clinical_tsv}" ] && [ \$(wc -l < "${clinical_tsv}") -gt 1 ]; then
        igv_snapshots.py \\
            --mode somatic \\
            --sample ${meta.id} \\
            --sites-tsv ${clinical_tsv} \\
            --bam ${hg38_bam} \\
            --fasta ${hg38_fasta} \\
            --out-html ${meta.id}/somatic/${meta.id}_igv_report.html \\
            --flanking ${flanking}
    else
        echo "No clinical SNVs for ${meta.id}; skipping somatic snapshots" >&2
    fi

    # Focal duplications. FLT3-ITD and UBTF-TD are frequently the whole finding
    # in AML and are the calls most in need of eyeballing, since a duplication
    # at an unexpected position is the artefact signature. The duplication
    # tables carry chrom and ref_pos, which is all somatic mode needs, so the
    # rows are reshaped into a sites table rather than adding a third mode.
    #
    # hg38 alignment: both duplication callers that produce ref_pos run there.
    make_focal_dup_sites.py ${meta.id}.focal_dup_sites.tsv ${dup_tsvs}

    if [ -s "${meta.id}.focal_dup_sites.tsv" ] && [ \$(wc -l < "${meta.id}.focal_dup_sites.tsv") -gt 1 ]; then
        igv_snapshots.py \\
            --mode somatic \\
            --sample ${meta.id} \\
            --sites-tsv ${meta.id}.focal_dup_sites.tsv \\
            --bam ${hg38_bam} \\
            --fasta ${hg38_fasta} \\
            --out-html ${meta.id}/focal_dup/${meta.id}_focal_dup_igv.html \\
            --flanking ${flanking}
    else
        echo "No focal duplication calls for ${meta.id}; skipping" >&2
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        igv_snapshots: v1.1
        igv_reports: \$(create_report --version 2>/dev/null || echo "unknown")
        samtools: \$(samtools --version | head -1 | awk '{print \$NF}')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}/translocations ${meta.id}/somatic ${meta.id}/focal_dup
    touch ${meta.id}/translocations/${meta.id}.translocations.html
    touch ${meta.id}/translocations/${meta.id}.translocations.manifest.json
    touch ${meta.id}/somatic/${meta.id}_igv_report.html
    echo '"${task.process}": stub' > versions.yml
    """
}
