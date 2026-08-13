process REALIGN_T2T {
    tag      "${meta.id}"
    label    'process_xhigh'
    label    'process_long'

    input:
    tuple val(meta), path(minknow_bam)

    output:
    tuple val(meta), path("${meta.id}.t2t.bam"), path("${meta.id}.t2t.bam.bai"), emit: bam_bai
    path "versions.yml",                                                          emit: versions

    script:
    // Critical: NO -T on samtools fastq, NO -y on minimap2.
    // Stale NC_-named SA tags from the MinKNOW BAM would otherwise leak through
    // and crash CuteSV/Sniffles when they try to fetch NC_ contigs from the
    // chr-named reference.
    //
    // BAM is indexed here (rather than in a separate SAMTOOLS_INDEX step) so
    // that BAM and BAI both live in this module's work dir. Downstream
    // consumers that do readlink -f on the BAM (Clair3's wrapper does this)
    // look for the BAI next to the resolved BAM path; colocation guarantees
    // it's there.
    def aligner_target = params.t2t_mmi ?: params.t2t_fasta
    """
    set -euo pipefail

    samtools fastq -@ 4 ${minknow_bam} \\
        | minimap2 -ax map-ont -t ${task.cpus} --MD --secondary=no ${aligner_target} - \\
        | samtools sort -@ 4 -m 2G -o ${meta.id}.t2t.bam -

    samtools index -@ 8 ${meta.id}.t2t.bam

    # Post-realign QC: SA tags should not reference NC_ contigs.
    n_nc=\$(samtools view ${meta.id}.t2t.bam 2>/dev/null | grep -c "SA:Z:.*NC_" || true)
    echo "Reads with NC_-referencing SA tags: \$n_nc (must be 0)"
    if [ "\$n_nc" -gt 0 ]; then
        echo "WARNING: stale NC_ contig refs in SA tags; downstream callers may fail." >&2
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(samtools --version | head -1 | awk '{print \$NF}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.t2t.bam ${meta.id}.t2t.bam.bai
    echo '"${task.process}": stub' > versions.yml
    """
}
