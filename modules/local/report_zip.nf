process REPORT_ZIP {
    tag      "zip"
    label    'process_low'

    // A relative publishDir resolves against launchDir, which is where
    // REPORT_BUNDLE already writes the bundle.
    publishDir path: { '.' }, mode: 'copy', overwrite: true, pattern: '*.zip'

    // Packages the finished bundle for release. Distinct from the bundle
    // itself: the bundle is a working tree the dashboard is built against,
    // this is the artefact that leaves the server.
    //
    // Internal path length is checked before handover, because Windows
    // Explorer still refuses to extract beyond 260 characters and the tree
    // nests <sample>/igv/translocations/<event>.A.html.

    input:
    path bundle

    output:
    path "*.zip",        emit: zip
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def light = params.report_zip_light ? "--light" : ""
    def force = params.report_zip_force ? "--force" : ""
    def nobam = params.report_zip_nobam ? "--no-bam" : ""
    """
    set -euo pipefail

    # --out . keeps the archive in the task directory so publishDir controls
    # where it lands, rather than the script writing outside its workdir.
    make_report_zip.sh "${bundle}" ${light} ${nobam} ${force} --out .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        make_report_zip: v1.1
        zip: \$(zip --version 2>/dev/null | sed -n '2s/.*Zip \\([0-9.]*\\).*/\\1/p' || echo "unknown")
    END_VERSIONS
    """

    stub:
    """
    touch ${bundle}.zip
    echo '"${task.process}": stub' > versions.yml
    """
}
