process MERGE_MINKNOW {
    tag   "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(minknow_input)

    output:
    tuple val(meta), path("${meta.id}.minknow.merged.bam"), emit: merged_bam
    path "versions.yml",                                    emit: versions

    script:
    // The P2i / MinKNOW writes one BAM per ~1.5 h of sequencing, so a sample's
    // reads arrive as a FOLDER of many BAMs. This step concatenates them into
    // a single BAM that both realign steps then consume.
    //
    // We use `samtools cat`, NOT `samtools merge`:
    //   - All chunk BAMs come from the same run with identical @SQ headers, so
    //     concatenation is valid (samtools cat requires matching references).
    //   - cat just glues BGZF blocks: fast, no coordinate re-sort.
    //   - The next pipeline step is `samtools fastq`, which reads every record
    //     regardless of order, so a globally-sorted merge would be wasted work.
    //
    // The input path may be either a DIRECTORY (real P2i output) or a SINGLE
    // BAM file (the validation sample sheets point at one realigned BAM). Both
    // are handled so existing sheets keep working.
    """
    set -euo pipefail

    # Collect the input BAM(s) into a file list.
    bam_list=\$(mktemp)
    if [ -d "${minknow_input}" ]; then
        # Directory: take every *.bam at its top level. Point the sample sheet
        # at the folder that directly contains the pass BAMs (e.g. .../bam_pass).
        find "${minknow_input}/" -maxdepth 1 -name '*.bam' | sort > "\$bam_list"
    else
        # Single file: pass through (cat of one BAM just copies the blocks).
        echo "${minknow_input}" > "\$bam_list"
    fi

    n_bam=\$(wc -l < "\$bam_list")
    if [ "\$n_bam" -eq 0 ]; then
        echo "ERROR: no .bam files found in input '${minknow_input}' for sample ${meta.id}" >&2
        exit 1
    fi
    echo "Merging \$n_bam BAM chunk(s) for sample ${meta.id}"

    # Concatenate. -b takes a file-of-filenames. No sort, no index needed.
    samtools cat -b "\$bam_list" -o ${meta.id}.minknow.merged.bam

    # Sanity: confirm the merged BAM is structurally valid and non-empty.
    samtools quickcheck ${meta.id}.minknow.merged.bam
    n_reads=\$(samtools view -c ${meta.id}.minknow.merged.bam)
    echo "Merged read count for ${meta.id}: \$n_reads"
    if [ "\$n_reads" -eq 0 ]; then
        echo "ERROR: merged BAM for ${meta.id} contains zero reads" >&2
        exit 1
    fi

    rm -f "\$bam_list"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | awk '{print \$NF}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.minknow.merged.bam
    echo '"${task.process}": stub' > versions.yml
    """
}
