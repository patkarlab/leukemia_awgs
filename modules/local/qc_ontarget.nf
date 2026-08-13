process QC_ONTARGET {
    tag      "${meta.id}"
    label    'process_low'

    // conda and publishDir are set in conf/envs.config and
    // conf/modules.config, keeping all per-process wiring in one place.

    input:
    tuple val(meta), path(bam), path(bai), path(panel_bed)

    output:
    tuple val(meta), path("${meta.id}.region_coverage.tsv"),  emit: coverage
    tuple val(meta), path("${meta.id}.readlen_qscore.tsv"),   emit: summary
    tuple val(meta), path("${meta.id}.region_coverage.png"),  emit: coverage_png
    tuple val(meta), path("${meta.id}.readlen_hist.png"),     emit: readlen_png
    tuple val(meta), path("${meta.id}.qscore_hist.png"),      emit: qscore_png
    path "versions.yml",                                       emit: versions

    script:
    """
    set -euo pipefail

    qc_al_sample.py \\
        --bam ${bam} \\
        --bed ${panel_bed} \\
        --sample ${meta.id} \\
        --outdir . \\
        --threshold ${params.qc_depth_threshold} \\
        --threads ${task.cpus}

    # per-region coverage bar chart (per-sample median + threshold reference lines)
    plot_region_coverage.py \\
        --out ${meta.id}.region_coverage.png \\
        --threshold ${params.qc_depth_threshold} \\
        ${meta.id}.region_coverage.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1 | awk '{print \$NF}')
        qc_al_sample.py: v0.1
        plot_region_coverage.py: v0.1
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.region_coverage.tsv
    touch ${meta.id}.readlen_qscore.tsv
    touch ${meta.id}.region_coverage.png
    touch ${meta.id}.readlen_hist.png
    touch ${meta.id}.qscore_hist.png
    echo '"${task.process}": stub' > versions.yml
    """
}
