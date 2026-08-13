process REALIGN_HG38 {
    tag      "${meta.id}"
    label    'process_xhigh'
    label    'process_long'

    input:
    tuple val(meta), path(minknow_bam)

    output:
    tuple val(meta), path("${meta.id}.hg38.bam"), path("${meta.id}.hg38.bam.bai"), emit: bam_bai
    path "versions.yml",                                                            emit: versions

    script:
    // Critical difference from REALIGN_T2T:
    //   - samtools fastq -T '*' carries ALL aux tags via FASTQ comments
    //   - minimap2 -y re-emits those comment tags into the new BAM
    // Result: methylation tags (MM/ML), read groups, etc. all survive the
    // round-trip. This is intentional for the hg38 track — the SA-tag
    // staleness issue that motivates dropping these flags on T2T does not
    // apply here, since the T2T BAM already has chr-named SA tags.
    //
    // We also samtools-index the BAM here (in the same work dir) instead of
    // in a separate SAMTOOLS_INDEX_HG38 step. Reason: when bam and bai live
    // in different upstream work dirs, downstream consumers that do
    // readlink -f on the BAM (e.g. Clair3's run_clair3.sh wrapper) look for
    // the BAI next to the resolved BAM path and don't find it. Colocating
    // both in this module's work dir fixes that — downstream symlinks
    // resolve to the same upstream dir for both files.
    def aligner_target = params.hg38_mmi && file(params.hg38_mmi).exists() ? params.hg38_mmi : params.hg38_fasta
    """
    set -euo pipefail

    # Capture input primary read count for the post-realign QC check.
    INPUT_PRIMARY=\$(samtools view -c -F 0x900 ${minknow_bam})
    echo "Input primary read count: \$INPUT_PRIMARY"

    samtools fastq -T '*' ${minknow_bam} \\
        | minimap2 -ax map-ont -y --MD --secondary=no -t ${task.cpus} ${aligner_target} - \\
        | samtools sort -@ 8 -o ${meta.id}.hg38.bam -

    # Index in the same work dir so downstream symlinks resolve together.
    samtools index -@ 8 ${meta.id}.hg38.bam

    # Validate output BAM structure.
    samtools quickcheck ${meta.id}.hg38.bam

    # Primary-read-count must match input. A delta means reads were lost
    # (unlikely with map-ont) or the pipe was truncated.
    OUTPUT_PRIMARY=\$(samtools view -c -F 0x900 ${meta.id}.hg38.bam)
    echo "Output primary read count: \$OUTPUT_PRIMARY"
    if [ "\$INPUT_PRIMARY" -ne "\$OUTPUT_PRIMARY" ]; then
        echo "WARNING: primary read count mismatch: delta = \$((OUTPUT_PRIMARY - INPUT_PRIMARY))" >&2
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(samtools --version | head -1 | awk '{print \$NF}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.hg38.bam ${meta.id}.hg38.bam.bai
    echo '"${task.process}": stub' > versions.yml
    """
}
