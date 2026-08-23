process ICHORCNA {
    tag      "${meta.id}"
    label    'process_medium'

    input:
    tuple val(meta), path(bam), path(bai), path(panel_bed)

    output:
    tuple val(meta), path("ichorcna_out"),                                            emit: outdir
    tuple val(meta), path("ichorcna_out/${meta.id}.params.txt"), optional: true, emit: params_file
    tuple val(meta), path("ichorcna_out/${meta.id}.seg.txt"),    optional: true, emit: seg
    tuple val(meta), path("ichorcna_out/${meta.id}.cna.seg"),    optional: true, emit: cna_seg
    path "versions.yml",                                                              emit: versions

    script:
    // Native conda env at params.ichorcna_env_prefix — DO NOT activate.
    // The docker image (molhemat/ichorcna:1.3) ships R 3.5.2 which can't load
    // the bundled optparse package. The native conda env has R 4.0.3 +
    // ichorCNA 0.3.4 working cleanly.
    //
    // Trick from production step 11: the panel BED is concatenated with the
    // bundled centromere file and passed as --centromere. This masks the
    // adaptive-sampling on-target regions out of ichorCNA's bin analysis,
    // since those regions are heavily over-enriched relative to the rest of
    // the genome and would otherwise distort the HMM segmentation.
    """
    set -euo pipefail
    mkdir -p ichorcna_out

    READCOUNTER=${params.ichorcna_env_prefix}/bin/readCounter
    RSCRIPT=${params.ichorcna_env_prefix}/bin/Rscript

    # Build the exclude regions file (centromeres + panel mask)
    {
        awk -v OFS='\\t' '!/^#/ && NF>=3 {print \$1, \$2, \$3}' ${params.ichorcna_centromere}
        awk -v OFS='\\t' '!/^#/ && NF>=3 {print \$1, \$2, \$3}' ${params.ichorcna_exclude}
        awk -v OFS='\\t' '!/^#/ && NF>=3 {print \$1, \$2, \$3}' ${panel_bed}
    } | sort -k1,1 -k2,2n > exclude.txt

    # Stage 1: readCounter binning at 1 Mb
    "\$READCOUNTER" \\
        --window ${params.ichorcna_bin_size} \\
        --quality ${params.ichorcna_qual} \\
        --chromosome ${params.ichorcna_chrs} \\
        ${bam} > ${meta.id}.wig

    n_bins=\$(grep -vc "^[fv]" ${meta.id}.wig || echo 0)
    echo "WIG bins: \$n_bins"
    if [ "\$n_bins" -lt 100 ]; then
        echo "ERROR: too few WIG bins (\$n_bins)" >&2
        exit 1
    fi

    # Stage 2: ichorCNA HMM CN calling
    "\$RSCRIPT" ${params.ichorcna_run_script} \\
        --id ${meta.id} \\
        --WIG ${meta.id}.wig \\
        --ploidy "${params.ichorcna_ploidy}" \\
        --normal "${params.ichorcna_normal}" \\
        --maxCN ${params.ichorcna_max_cn} \\
        --gcWig ${params.ichorcna_gc_wig} \\
        --mapWig ${params.ichorcna_map_wig} \\
        --normalPanel ${params.ichorcna_pon_rds} \\
        --centromere exclude.txt \\
        --includeHOMD False \\
        --chrs "c(1:22)" \\
        --chrTrain "${params.ichorcna_chrtrain}" \\
        --estimateNormal TRUE \\
        --estimatePloidy FALSE \\
        --estimateScPrevalence FALSE \\
        --scStates 'c()' \\
        --txnE 0.9999999 \\
        --txnStrength 10000000 \\
        --outDir ichorcna_out/

    # Move expected outputs into the publish dir at predictable paths
    if [ -f ichorcna_out/${meta.id}.params.txt ]; then
        tf=\$(awk -F'\\t' '/Tumor Fraction:/ {print \$2}' ichorcna_out/${meta.id}.params.txt | head -1)
        ploidy=\$(awk -F'\\t' '/Ploidy:/ {print \$2}' ichorcna_out/${meta.id}.params.txt | head -1)
        echo "Tumor fraction: \${tf:-?}  Ploidy: \${ploidy:-?}"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ichorCNA: \$("\$RSCRIPT" -e 'cat(as.character(packageVersion("ichorCNA")))' 2>/dev/null || echo unknown)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ichorcna_out
    touch ichorcna_out/${meta.id}.params.txt ichorcna_out/${meta.id}.seg.txt ichorcna_out/${meta.id}.cna.seg
    echo '"${task.process}": stub' > versions.yml
    """
}
