process FOCAL_TANDEM_DUP {
    tag   "${meta.id} ${target.label}"
    label 'process_low'

    input:
    tuple val(meta), path(bam), path(bai), path(panel_bed), val(target)
    path fasta
    path fai
    path hotspot_bed

    output:
    tuple val(meta), path("${meta.id}.${target.gene}.tandem_dup.tsv"), emit: tsv
    path "versions.yml",                                               emit: versions

    script:
    // Short internal tandem duplications, detected as CIGAR insertions and
    // confirmed against adjacent reference sequence. The gene, its size range
    // and its support threshold all come from the focal-duplication targets
    // table; the interval comes from the panel BED. No coordinate or gene
    // name is written into the pipeline.
    def hotspot_arg = hotspot_bed.name != 'NO_FILE' ? "--hotspot-bed ${hotspot_bed}" : ''
    """
    call_tandem_dup.py \\
        --bam       ${bam} \\
        --fasta     ${fasta} \\
        --gene      ${target.gene} \\
        --label     ${target.label} \\
        --panel-bed ${panel_bed} \\
        ${hotspot_arg} \\
        --sample    ${meta.id} \\
        --output    ${meta.id}.${target.gene}.tandem_dup.tsv \\
        --min-len       ${target.min_len} \\
        --max-len       ${target.max_len} \\
        --min-support   ${target.min_support} \\
        --min-mapq      ${params.focal_dup_min_mapq} \\
        --min-spanning-depth ${params.focal_dup_min_spanning_depth} \\
        --min-tandem-identity ${params.focal_dup_min_tandem_identity}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        call_tandem_dup.py: \$(call_tandem_dup.py --version 2>&1 | awk '{print \$NF}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.${target.gene}.tandem_dup.tsv
    echo '"${task.process}": stub' > versions.yml
    """
}
