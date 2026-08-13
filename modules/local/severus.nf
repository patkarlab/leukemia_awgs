process SEVERUS {
    tag      "${meta.id}"
    label    'process_medium'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("severus_out/${meta.id}.severus.vcf"), emit: vcf
    tuple val(meta), path("severus_out"),                         emit: outdir
    path "versions.yml",                                            emit: versions

    script:
    // Severus writes its output VCF under one of several possible subpaths
    // depending on version and the somatic vs all-SV mode. We normalize by
    // copying whichever exists to a single predictable filename so
    // SURVIVOR_MERGE downstream has a stable input path.
    """
    set -euo pipefail
    mkdir -p severus_out

    severus \\
        --target-bam ${bam} \\
        --out-dir severus_out \\
        --threads ${params.severus_threads} \\
        --vntr-bed ${params.severus_vntr_bed} \\
        --PON ${params.severus_pon} \\
        --min-support ${params.severus_min_support} \\
        --min-mapq ${params.severus_min_mapq}

    # Normalize: find whichever VCF exists, copy it to a canonical name.
    normalized="severus_out/${meta.id}.severus.vcf"
    for cand in \\
        severus_out/somatic_SVs/severus_somatic.vcf \\
        severus_out/all_SVs/severus_all.vcf \\
        severus_out/severus_somatic.vcf \\
        severus_out/severus_all.vcf
    do
        if [ -s "\$cand" ]; then
            cp "\$cand" "\$normalized"
            break
        fi
    done

    if [ ! -s "\$normalized" ]; then
        echo "ERROR: no Severus VCF produced under severus_out/" >&2
        ls -lR severus_out/ >&2 || true
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        severus: \$(severus --version 2>&1 | awk '{print \$NF}')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p severus_out
    touch severus_out/${meta.id}.severus.vcf
    echo '"${task.process}": stub' > versions.yml
    """
}
