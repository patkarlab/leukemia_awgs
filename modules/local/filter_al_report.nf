process FILTER_AL_REPORT {
    tag   "${meta.id} [${meta.panel}]"
    label 'process_low'

    input:
    tuple val(meta), path(candidates_tsv), path(panel_bed)

    output:
    tuple val(meta), path("al_filtered/*.al_clinical.tsv"),      optional: true, emit: clinical
    tuple val(meta), path("al_filtered/*.al_filtered.tsv"),      optional: true, emit: filtered
    tuple val(meta), path("al_filtered/al_filter_summary.tsv"),  optional: true, emit: summary
    path "versions.yml", emit: versions

    script:
    def include_ig = params.al_include_ig ? '--include-ig' : ''
    """
    filter_al_somatic_candidates.py \\
        --panel-bed  ${panel_bed} \\
        --input      ${candidates_tsv} \\
        --outdir     al_filtered \\
        --max-pop-af ${params.al_max_pop_af} \\
        --aliases       ${params.gene_aliases} \\
        --excluded-loci ${params.excluded_loci} \\
        ${include_ig}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        filter_al_somatic_candidates.py: v0.1
    END_VERSIONS
    """

    stub:
    """
    mkdir -p al_filtered
    touch al_filtered/${meta.id}.al_clinical.tsv
    echo '"${task.process}": stub' > versions.yml
    """
}
