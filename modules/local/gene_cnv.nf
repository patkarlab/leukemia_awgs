process GENE_CNV {
    tag   "${meta.id} [${meta.panel}]"
    label 'process_medium'

    input:
    tuple val(meta), path(bam), path(bai), path(panel_bed), path(ref_bams), path(ref_bais)
    path gff

    output:
    tuple val(meta), path("${meta.id}.chrom_cnv.tsv"),          emit: chrom
    tuple val(meta), path("${meta.id}.gene_cnv.tsv"),           emit: genes
    tuple val(meta), path("${meta.id}.gene_cnv_segments.tsv"),  emit: segments
    tuple val(meta), path("${meta.id}.gene_cnv_plots"),         emit: plots, optional: true
    path "versions.yml",                                        emit: versions

    script:
    // Run on T2T: adaptive sampling enriched against that panel BED, so
    // on-target depth is defined in those coordinates.
    //
    // Other samples in the run are passed as references. They are used only
    // for intragenic segments and to tell a homozygous deletion from an
    // unenriched region, both of which need a cohort; chromosome and gene
    // calls are normalised within the sample so a copy-number change in a
    // reference member cannot propagate into them.
    def refs = ref_bams.findAll { it.name != 'NO_FILE' }
                       .collect { "--reference-bam ${it}" }.join(' ')
    def gff_arg = gff.name != 'NO_FILE' ? "--gff ${gff}" : ''
    """
    call_gene_cnv.py \\
        --bam        ${bam} \\
        --panel-bed  ${panel_bed} \\
        --sample     ${meta.id} \\
        --out-prefix ${meta.id} \\
        ${gff_arg} ${refs} \\
        --binsize                ${params.gene_cnv_binsize} \\
        --min-mapq               ${params.gene_cnv_min_mapq} \\
        --threshold              ${params.gene_cnv_threshold} \\
        --seg-z                  ${params.gene_cnv_seg_z} \\
        --min-bins               ${params.gene_cnv_min_bins} \\
        --exon-edge-tol          ${params.gene_cnv_exon_edge_tol} \\
        --min-chrom-regions      ${params.gene_cnv_min_chrom_regions} \\
        --plot-genes             ${params.gene_cnv_plot_genes}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        call_gene_cnv.py: \$(call_gene_cnv.py --version 2>&1 | awk '{print \$NF}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.chrom_cnv.tsv ${meta.id}.gene_cnv.tsv ${meta.id}.gene_cnv_segments.tsv
    mkdir -p ${meta.id}.gene_cnv_plots
    echo '"${task.process}": stub' > versions.yml
    """
}
