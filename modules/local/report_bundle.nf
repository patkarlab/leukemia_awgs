process REPORT_BUNDLE {
    tag      "bundle"
    label    'process_low'

    input:
    // These inputs are consumed purely to establish ordering. The bundle is
    // assembled by scanning the published output tree, which only exists once
    // publishDir has run for every producing process, so the process must not
    // start until those have completed. Nextflow provides no direct "after
    // publish" dependency, so the completion signals are taken as inputs and
    // ignored in the script body.
    val  ready_signals
    path bundle_script
    // The panel BEDs are staged rather than resolved relative to the script,
    // because inside a work directory bin/ is staged alone with no assets/
    // beside it. Without them the alignment slicing skips silently and the
    // bundle ships without the reads its calls were made from.
    path panel_bed_t2t
    path panel_bed_hg38
    val  outdir
    val  bundle_name

    output:
    path "${bundle_name}",        emit: bundle
    path "${bundle_name}.zip",    optional: true, emit: archive
    path "${bundle_name}.tar.gz", optional: true, emit: tarball
    path "versions.yml",          emit: versions

    script:
    // The results directory is resolved to an absolute path because the script
    // runs from the task work directory, not the launch directory.
    def results_dir = outdir.toString().startsWith('/')
        ? outdir.toString()
        : "${workflow.launchDir}/${outdir}"
    """
    set -euo pipefail

    export PANEL_BED_T2T="\$(readlink -f ${panel_bed_t2t})"
    export PANEL_BED_HG38="\$(readlink -f ${panel_bed_hg38})"

    if [ ! -d "${results_dir}" ]; then
        echo "ERROR: results directory not found: ${results_dir}" >&2
        echo "The bundle is assembled from published outputs, so it cannot run" >&2
        echo "before publishDir has written them." >&2
        exit 1
    fi

    bash ${bundle_script} "${results_dir}" "${bundle_name}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        build_report_bundle: v1.0
        zip: \$(zip --version 2>/dev/null | sed -n '2s/.*Zip \\([0-9.]*\\).*/\\1/p' || echo "not available")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${bundle_name}
    touch ${bundle_name}/.stub
    echo '"${task.process}": stub' > versions.yml
    """
}
