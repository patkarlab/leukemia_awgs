process AUGMENT_SV_SUPPORT {
    tag "$meta.id"
    label 'process_low'
    conda params.conda_awgs_sv

    publishDir "${params.outdir}/t2t/calls/mm_annotated", mode: 'copy', pattern: "*.mm_annotated.tsv"

    input:
    // annotated staged under a different name to avoid an output-name collision
    tuple val(meta),
          path(annotated, stageAs: 'annotated_in.tsv'),
          path(sniffles_vcf),
          path(cutesv_vcf),
          path(severus_vcf)

    output:
    tuple val(meta), path("${meta.id}.mm_annotated.tsv"), emit: annotated
    path "versions.yml",                                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    augment_sv_support.py \\
        --annotated annotated_in.tsv \\
        --sniffles  ${sniffles_vcf} \\
        --cutesv    ${cutesv_vcf} \\
        --severus   ${severus_vcf} \\
        --output    ${meta.id}.mm_annotated.tsv \\
        --tol       ${params.support_tol}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}
