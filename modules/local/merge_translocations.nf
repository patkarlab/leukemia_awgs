process MERGE_TRANSLOCATIONS {
    tag "$meta.id"
    label 'process_low'
    conda params.conda_awgs_sv

    publishDir "${params.outdir}/t2t/calls/mm_annotated", mode: 'copy', pattern: "*.translocations.tsv"

    input:
    tuple val(meta), path(annotated)

    output:
    tuple val(meta), path("*.translocations.tsv"), emit: translocations
    path "versions.yml",                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    merge_translocations.py \\
        --input    ${annotated} \\
        --outdir   . \\
        --max-dist ${params.tra_max_dist}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}
