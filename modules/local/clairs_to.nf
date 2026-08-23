process CLAIRS_TO {
    tag      "${meta.id}"
    label    'process_medium'

    input:
    tuple val(meta), path(bam), path(bai), path(panel_bed)

    output:
    tuple val(meta), path("clairs_to_out"),                                                emit: outdir
    tuple val(meta), path("clairs_to_out/snv_${meta.id}.vcf.gz"),
                     path("clairs_to_out/snv_${meta.id}.vcf.gz.tbi"),                     emit: snv_vcf
    tuple val(meta), path("clairs_to_out/indel_${meta.id}.vcf.gz"),
                     path("clairs_to_out/indel_${meta.id}.vcf.gz.tbi"),                   emit: indel_vcf
    path "versions.yml",                                                                    emit: versions

    script:
    // --disable_verdict because current ClairS-TO docker image is missing the
    // G1000 loci files needed for Verdict's CNA-aware classification. PoN
    // filtering still runs and tags germline candidates as NonSomatic.
    def verdict_flag = params.clairs_to_disable_verdict ? '--disable_verdict' : ''
    """
    set -euo pipefail
    mkdir -p clairs_to_out

    docker run --rm \\
        --user \$(id -u):\$(id -g) \\
        -v /goast:/goast \\
        -v \$(pwd):/work \\
        -w /work \\
        ${params.clairs_to_image} \\
        /opt/bin/run_clairs_to \\
            --tumor_bam_fn ${bam} \\
            --ref_fn ${params.hg38_fasta} \\
            --bed_fn ${panel_bed} \\
            --threads ${params.clairs_to_threads} \\
            --platform ${params.clairs_to_platform} \\
            --output_dir clairs_to_out \\
            --sample_name ${meta.id} \\
            ${verdict_flag}

    # ClairS-TO produces snv_<sample>.vcf.gz and indel_<sample>.vcf.gz with
    # their .tbi indices. Verify and report variant counts.
    for kind in snv indel; do
        vcf=clairs_to_out/\${kind}_${meta.id}.vcf.gz
        if [ -s "\$vcf" ]; then
            n=\$(bcftools view -H "\$vcf" 2>/dev/null | wc -l)
            n_pass=\$(bcftools view -H -f PASS "\$vcf" 2>/dev/null | wc -l)
            echo "[\${kind}] \$n total / \$n_pass PASS"
        else
            echo "[\${kind}] expected VCF missing: \$vcf" >&2
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clairs-to: ${params.clairs_to_image}
    END_VERSIONS
    """

    stub:
    """
    mkdir -p clairs_to_out
    touch clairs_to_out/snv_${meta.id}.vcf.gz clairs_to_out/snv_${meta.id}.vcf.gz.tbi
    touch clairs_to_out/indel_${meta.id}.vcf.gz clairs_to_out/indel_${meta.id}.vcf.gz.tbi
    echo '"${task.process}": stub' > versions.yml
    """
}
