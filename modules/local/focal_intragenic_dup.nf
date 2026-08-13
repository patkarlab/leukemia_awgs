process FOCAL_INTRAGENIC_DUP {
    tag   "${meta.id} ${target.label}"
    label 'process_low'

    input:
    tuple val(meta), path(bam), path(bai), path(panel_bed), val(target)
    path exon_bed

    output:
    tuple val(meta), path("${meta.id}.${target.gene}.intragenic_dup.tsv"), emit: tsv
    path "versions.yml",                                                   emit: versions

    script:
    // Large intragenic duplications, from split-read junctions plus large
    // insertions plus exon coverage ratio. Run on the T2T alignment: KMT2A
    // sits in a segmentally duplicated stretch of 11q23 where hg38
    // alt-contigs pull away the supplementary alignments this call depends on.
    def exon_arg = exon_bed.name != 'NO_FILE' ? "--exon-bed ${exon_bed}" : ''
    """
    call_intragenic_dup.py \\
        --bam       ${bam} \\
        --gene      ${target.gene} \\
        --label     ${target.label} \\
        --panel-bed ${panel_bed} \\
        ${exon_arg} \\
        --sample    ${meta.id} \\
        --output    ${meta.id}.${target.gene}.intragenic_dup.tsv \\
        --min-span            ${target.min_len} \\
        --max-span            ${target.max_len} \\
        --min-split-support   ${target.min_support} \\
        --min-mapq            ${params.focal_dup_min_mapq} \\
        --coverage-ratio-threshold ${params.focal_dup_coverage_ratio_threshold}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        call_intragenic_dup.py: \$(call_intragenic_dup.py --version 2>&1 | awk '{print \$NF}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.${target.gene}.intragenic_dup.tsv
    echo '"${task.process}": stub' > versions.yml
    """
}
