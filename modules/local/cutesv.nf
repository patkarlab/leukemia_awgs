process CUTESV {
    tag      "${meta.id}"
    label    'process_medium'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.cutesv.t2t.vcf.gz"), path("${meta.id}.cutesv.t2t.vcf.gz.tbi"), emit: vcf
    path "versions.yml",                                                                              emit: versions

    script:
    // CuteSV's positional argument order is: bam ref out.vcf workdir.
    // It writes a plain VCF, so we bgzip+tabix manually afterwards.
    """
    set -euo pipefail

    mkdir -p cutesv_workdir

    cuteSV \\
        ${bam} ${params.t2t_fasta} ${meta.id}.cutesv.t2t.vcf cutesv_workdir \\
        --threads ${params.cutesv_threads} \\
        --sample ${meta.id} \\
        --min_support ${params.cutesv_min_support} \\
        --min_mapq ${params.cutesv_min_mapq} \\
        --min_size ${params.cutesv_min_size} \\
        --min_read_len ${params.cutesv_min_read_len} \\
        --max_cluster_bias_INS 100 \\
        --diff_ratio_merging_INS 0.3 \\
        --max_cluster_bias_DEL 100 \\
        --diff_ratio_merging_DEL 0.3 \\
        --genotype

    bgzip -f ${meta.id}.cutesv.t2t.vcf
    tabix -p vcf ${meta.id}.cutesv.t2t.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cuteSV: \$(cuteSV --version 2>&1 | awk '{print \$NF}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.cutesv.t2t.vcf.gz ${meta.id}.cutesv.t2t.vcf.gz.tbi
    echo '"${task.process}": stub' > versions.yml
    """
}
