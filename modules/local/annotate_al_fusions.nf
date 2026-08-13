process ANNOTATE_AL_FUSIONS {
    tag   "${meta.id} [${meta.panel}]"
    label 'process_low'

    input:
    tuple val(meta), path(merged_vcf), path(merged_tbi), path(panel_bed)
    path cytoband_bed
    path dictionary
    path anchors

    output:
    tuple val(meta), path("${meta.id}.al_fusions.tsv"), emit: tsv
    path "versions.yml",                                emit: versions

    script:
    """
    annotate_al_fusions.py \\
        --vcf          ${merged_vcf} \\
        --panel-bed    ${panel_bed} \\
        --cytoband-bed ${cytoband_bed} \\
        --dictionary   ${dictionary} \\
        --anchors      ${anchors} \\
        --panel        ${meta.panel} \\
        --sample       ${meta.id} \\
        --output       ${meta.id}.al_fusions.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        annotate_al_fusions.py: \$(annotate_al_fusions.py --version 2>&1 | awk '{print \$NF}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.al_fusions.tsv
    echo '"${task.process}": stub' > versions.yml
    """
}
