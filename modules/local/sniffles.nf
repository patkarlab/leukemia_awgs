process SNIFFLES {
    tag      "${meta.id}"
    label    'process_medium'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.sniffles.t2t.vcf.gz"), path("${meta.id}.sniffles.t2t.vcf.gz.tbi"), emit: vcf
    path "versions.yml",                                                                                  emit: versions

    script:
    """
    set -euo pipefail

    sniffles \\
        --input ${bam} \\
        --reference ${params.t2t_fasta} \\
        --vcf ${meta.id}.sniffles.t2t.vcf.gz \\
        --threads ${params.sniffles_threads} \\
        --minsupport ${params.sniffles_min_support} \\
        --mapq ${params.sniffles_min_mapq} \\
        --min-alignment-length ${params.sniffles_min_aln_len} \\
        --output-rnames

    # Sniffles emits .vcf.gz but doesn't always emit the tabix index — make sure one exists
    if [ ! -f ${meta.id}.sniffles.t2t.vcf.gz.tbi ]; then
        tabix -p vcf ${meta.id}.sniffles.t2t.vcf.gz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sniffles: \$(sniffles --version 2>&1 | awk 'NR==1{print \$NF}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.sniffles.t2t.vcf.gz ${meta.id}.sniffles.t2t.vcf.gz.tbi
    echo '"${task.process}": stub' > versions.yml
    """
}
