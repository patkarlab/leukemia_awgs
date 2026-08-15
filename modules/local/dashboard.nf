process DASHBOARD {
    tag      "dashboard"
    label    'process_low'

    input:
    // The bundle directory is staged so that this process runs after
    // REPORT_BUNDLE and operates on a real input rather than on a path handed
    // in as a bare value. The builder writes the report HTML into the bundle
    // tree it is given, so the staged copy is republished with the reports
    // included.
    path bundle
    path builder_dir

    output:
    path "${bundle}",    emit: bundle
    path "versions.yml", emit: versions

    script:
    def subdir_arg = params.dashboard_subdir ? "--subdir ${params.dashboard_subdir}" : ""
    """
    set -euo pipefail

    # The builder is read-only against the run directory apart from writing the
    # report HTML and patching the IGV report's hash router, both of which are
    # idempotent. It is pointed at the staged bundle so that the per-sample
    # layout it expects is the one the bundle already provides.
    python3 ${builder_dir}/build.py "${bundle}" ${subdir_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        jinja2: \$(python3 -c "import jinja2; print(jinja2.__version__)" 2>/dev/null || echo "unknown")
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)" 2>/dev/null || echo "unknown")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${bundle}
    touch ${bundle}/cohort_index.html
    echo '"${task.process}": stub' > versions.yml
    """
}
