process CLAIR3_PHASED {
    tag      "${meta.id}"
    label    'process_xhigh'
    label    'process_long'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("clair3_out"),                                                                emit: outdir
    tuple val(meta), path("clair3_out/merge_output.vcf.gz"),
                     path("clair3_out/merge_output.vcf.gz.tbi"),                                        emit: merge_output
    tuple val(meta), path("clair3_out/phased_merge_output.vcf.gz"),
                     path("clair3_out/phased_merge_output.vcf.gz.tbi"),     optional: true,            emit: phased
    path "versions.yml",                                                                                  emit: versions

    script:
    // Model location:
    //   - SUP v500 ships inside the docker at /opt/models/...
    //   - all other models live on the host under params.clair3_model_dir/<name>/
    //     containing pileup.pt and full_alignment.pt
    // We bind-mount /goast into the container so host paths are also valid
    // container paths.
    def bundled_models = ['r1041_e82_400bps_sup_v500']
    def is_bundled     = bundled_models.contains(params.clair3_model_name)
    def model_path     = is_bundled \
        ? "/opt/models/${params.clair3_model_name}" \
        : "${params.clair3_model_dir}/${params.clair3_model_name}"
    """
    set -euo pipefail
    mkdir -p clair3_out

    # Best-effort basecaller-vs-model consistency check (warn only)
    bc_line=\$(samtools view -H ${bam} 2>/dev/null | grep -i basecall_model | head -1 || true)
    if [ -n "\$bc_line" ]; then
        echo "BAM basecaller line: \$bc_line"
        bc_class="unknown"
        echo "\$bc_line" | grep -qi "_hac@" && bc_class="hac"
        echo "\$bc_line" | grep -qi "_sup@" && bc_class="sup"
        model_class="unknown"
        echo "${params.clair3_model_name}" | grep -q "_hac_" && model_class="hac"
        echo "${params.clair3_model_name}" | grep -q "_sup_" && model_class="sup"
        if [ "\$bc_class" != "unknown" ] && [ "\$model_class" != "unknown" ] && [ "\$bc_class" != "\$model_class" ]; then
            echo "WARNING: basecaller is \$bc_class but model is \$model_class. Variant accuracy may degrade." >&2
        fi
    fi

    docker run --rm \\
        --user \$(id -u):\$(id -g) \\
        -v /goast:/goast \\
        -v \$(pwd):/work \\
        -w /work \\
        ${params.clair3_image} \\
        /opt/bin/run_clair3.sh \\
            --bam_fn=${bam} \\
            --ref_fn=${params.hg38_fasta} \\
            --output=clair3_out \\
            --threads=${params.clair3_threads} \\
            --platform=${params.clair3_platform} \\
            --model_path=${model_path} \\
            --enable_phasing \\
            --longphase_for_phasing \\
            --sample_name=${meta.id}

    # Ensure expected outputs are present and tabix-indexed
    for v in clair3_out/merge_output.vcf.gz clair3_out/phased_merge_output.vcf.gz; do
        if [ -s "\$v" ] && [ ! -f "\$v.tbi" ]; then
            tabix -p vcf "\$v"
        fi
    done

    if [ ! -s clair3_out/merge_output.vcf.gz ]; then
        echo "ERROR: clair3 merge_output.vcf.gz missing" >&2
        ls -la clair3_out/ >&2 || true
        exit 1
    fi

    # Quick summary
    total=\$(bcftools view -H clair3_out/merge_output.vcf.gz | wc -l)
    pass=\$(bcftools view -f PASS -H clair3_out/merge_output.vcf.gz | wc -l)
    echo "Clair3 totals: \$total records / \$pass PASS"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: ${params.clair3_image}
        model: ${params.clair3_model_name}
    END_VERSIONS
    """

    stub:
    """
    mkdir -p clair3_out
    touch clair3_out/merge_output.vcf.gz clair3_out/merge_output.vcf.gz.tbi
    touch clair3_out/phased_merge_output.vcf.gz clair3_out/phased_merge_output.vcf.gz.tbi
    echo '"${task.process}": stub' > versions.yml
    """
}
