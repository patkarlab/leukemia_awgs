process EMBED_REPORT_ASSETS {
    tag      "embed"
    label    'process_low'

    // A relative publishDir resolves against launchDir, which is where
    // REPORT_BUNDLE already writes the bundle.
    publishDir path: { '.' }, mode: 'copy', overwrite: true

    // Rewrites the generated reports so every local dependency is carried
    // inside them. The builder emits relative references to its stylesheets,
    // scripts and figures; those resolve on the server and break the moment a
    // report is copied out of the bundle or extracted to a different depth on
    // another machine, silently and with no indication of what went missing.
    //
    // The IGV pages are deliberately left as separate documents. Each already
    // carries its own alignment slice, so folding them in would multiply the
    // report by the number of events.

    input:
    path bundle

    output:
    path "${bundle}",    emit: bundle
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def embed_igv = params.report_embed_igv ?: 0
    def igv_arg   = embed_igv > 0 ? "--embed-igv ${embed_igv}" : ""
    """
    set -euo pipefail

    # Runs in place on the staged bundle, which is then republished with the
    # rewritten reports. The closing report names any asset that was
    # referenced but not found; that list is the bundle step's to-do.
    embed_report_assets.py "${bundle}" ${igv_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        embed_report_assets: v1.1
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${bundle}
    echo '"${task.process}": stub' > versions.yml
    """
}
