/*
 * Panel path resolution.
 *
 * meta.panel is set once, in PREPARE_INPUT, and every panel-dependent file is
 * looked up through here. Keeping the lookup in one place is what allows a
 * single run to mix AML and ALL samples: nothing downstream reads a global
 * panel parameter.
 *
 * This deliberately returns a path rather than transforming a channel.
 * A helper that wrapped the .map() would not compile under Nextflow's strict
 * syntax, and the inline form at each call site is easier to follow anyway.
 */

def panelPath(meta, String key) {
    def name = (meta.panel ?: params.panel).toString().toUpperCase()
    def entry = params.panels[name]
    if( entry == null )
        error "Unknown panel '${name}' for sample '${meta.id}'. Known panels: ${params.panels.keySet().join(', ')}"
    def path = entry[key]
    if( path == null )
        error "Panel '${name}' has no '${key}' entry. Check conf/panels.config."
    return path
}

/*
 * Focal-duplication targets.
 *
 * Which genes are scanned for internal tandem or intragenic duplication, on
 * which reference, with which size and support thresholds, is data in
 * assets/focal_duplication_targets.tsv. Nothing here names a gene.
 *
 * Returns rows filtered to one reference and one detection mode; the caller
 * then pairs them with samples and drops rows out of scope for that sample's
 * panel.
 */
def focalTargets(String reference, String mode) {
    return Channel
        .fromPath(params.focal_dup_targets, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .filter { row ->
            row.reference?.toLowerCase() == reference &&
            row.mode?.toLowerCase() == mode
        }
}
