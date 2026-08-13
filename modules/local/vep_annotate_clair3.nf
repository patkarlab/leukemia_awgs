process VEP_ANNOTATE_CLAIR3 {
    tag      "${meta.id}"
    label    'process_medium'

    input:
    // From CLAIR3_PHASED.out.merge_output
    tuple val(meta), path(clair3_vcf), path(clair3_tbi)

    output:
    tuple val(meta), path("vep_out"),                                                emit: outdir
    tuple val(meta), path("vep_out/${meta.id}.all_annotated.tsv"),       optional: true, emit: all_tsv
    tuple val(meta), path("vep_out/${meta.id}.somatic_candidates.tsv"),  optional: true, emit: candidates_tsv
    path "versions.yml",                                                                  emit: versions

    script:
    // Mirrors step 12b's 4 stages:
    //   1. bcftools view -f PASS → pass.vcf.gz
    //   2. docker VEP v113 annotate (NO --hgvs / --hgvsg — known crash bug)
    //   3. bcftools +split-vep extract tidy TSV (19 columns, no HGVS)
    //   4. awk filter: pop_af < threshold AND non-synonymous → candidates TSV
    """
    set -euo pipefail
    mkdir -p vep_out

    pass_vcf=vep_out/${meta.id}.pass.vcf.gz
    ann_vcf=vep_out/${meta.id}.annotated.vcf.gz
    vep_html=vep_out/${meta.id}.vep_summary.html
    all_tsv=vep_out/${meta.id}.all_annotated.tsv
    cand_tsv=vep_out/${meta.id}.somatic_candidates.tsv

    # ----- Stage 1: PASS filter -----
    bcftools view -f PASS --threads ${params.vep_threads} \\
        -O z -o "\$pass_vcf" ${clair3_vcf}
    tabix -p vcf -f "\$pass_vcf"
    norm_vcf=vep_out/${meta.id}.norm.vcf.gz
    bcftools norm -m- -f ${params.hg38_fasta} "\$pass_vcf" -O z -o "\$norm_vcf"
    tabix -p vcf -f "\$norm_vcf"
    n_pass=\$(bcftools view -H "\$pass_vcf" | wc -l)
    echo "[stage 1] PASS variants: \$n_pass"

    if [ "\$n_pass" -eq 0 ]; then
        echo "No PASS variants — emitting empty TSVs and exiting."
        printf "chrom\\tpos\\tref\\talt\\tqual\\tvariant_type\\tgene\\ttranscript\\tbiotype\\tcanonical\\tconsequence\\timpact\\texon\\tdomains\\trs_id\\tpop_af_max\\tpop_af_max_source\\tclinvar_sig\\ttumor_af\\tREF_COUNT\\tALT_COUNT\\tDP\\thgvsc\\thgvsp\\n" > "\$all_tsv"
        cp "\$all_tsv" "\$cand_tsv"
        echo '"${task.process}": no-variants' > versions.yml
        exit 0
    fi

    # ----- Stage 2: VEP annotation via docker -----
    # --hgvs enabled; multi-allelics split by bcftools norm -m- above  [hgvs-annotation-applied-v2]
    docker run --rm \\
        --user \$(id -u):\$(id -g) \\
        -v /goast:/goast \\
        -v \$(pwd):/work \\
        -v ${params.vep_cache_dir}:/opt/vep/.vep \\
        -w /work \\
        ${params.vep_image} \\
        vep \\
            --input_file "\$norm_vcf" \\
            --output_file "\$ann_vcf" \\
            --vcf \\
            --compress_output bgzip \\
            --stats_file "\$vep_html" \\
            --cache \\
            --dir_cache /opt/vep/.vep \\
            --cache_version ${params.vep_cache_version} \\
            --assembly ${params.vep_assembly} \\
            --merged \\
            --offline \\
            --fasta ${params.hg38_fasta} \\
            --use_given_ref \\
            --hgvs \
            --pick \\
            --symbol --canonical --biotype --numbers --domains \\
            --check_existing \\
            --max_af --af \\
            --regulatory \\
            --fork ${params.vep_threads} \\
            --buffer_size 10000 \\
            --no_progress \\
            --force_overwrite

    tabix -p vcf -f "\$ann_vcf"
    echo "[stage 2] VEP done"

    # ----- Stage 3: tidy TSV via bcftools +split-vep -----
    export BCFTOOLS_PLUGINS=\${CONDA_PREFIX}/libexec/bcftools
    bcftools +split-vep \\
        "\$ann_vcf" \\
        -d \\
        -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t%SYMBOL\\t%Feature\\t%BIOTYPE\\t%CANONICAL\\t%Consequence\\t%IMPACT\\t%EXON\\t%DOMAINS\\t%Existing_variation\\t%MAX_AF\\t%MAX_AF_POPS\\t%CLIN_SIG\\t[%AF]\\t[%AD]\\t[%DP]\\t%HGVSc\\t%HGVSp\\n' \\
        -A tab \\
    | awk -v FS='\\t' -v OFS='\\t' '
        {
            if (length(\$3)==1 && length(\$4)==1) { vtype="SNV" } else { vtype="indel" }
            ad=\$19; dp=\$20; rc="-1"; ac="-1"
            if (ad != "" && ad != ".") { n=split(ad, a, ","); if (a[1] != "") rc=a[1]; if (n>=2 && a[2] != "") ac=a[2] }
            if (dp=="" || dp==".") dp="-1"
            print \$1, \$2, \$3, \$4, \$5, vtype, \$6, \$7, \$8, \$9, \$10, \$11, \$12, \$13, \$14, \$15, \$16, \$17, \$18, rc, ac, dp, \$21, \$22
        }
        ' > "\${all_tsv}.body"

    {
        printf "chrom\\tpos\\tref\\talt\\tqual\\tvariant_type\\tgene\\ttranscript\\tbiotype\\tcanonical\\tconsequence\\timpact\\texon\\tdomains\\trs_id\\tpop_af_max\\tpop_af_max_source\\tclinvar_sig\\ttumor_af\\tREF_COUNT\\tALT_COUNT\\tDP\\thgvsc\\thgvsp\\n"
        cat "\${all_tsv}.body"
    } > "\$all_tsv"
    rm -f "\${all_tsv}.body"

    n_all=\$(( \$(wc -l < "\$all_tsv") - 1 ))
    echo "[stage 3] Annotated TSV: \$n_all rows"

    # ----- Stage 4: somatic-candidate filter -----
    awk -v FS='\\t' -v OFS='\\t' \\
        -v af_max="${params.vep_gnomad_af_max}" \\
        -v exclude="${params.vep_exclude_consequences}" \\
        '
        NR==1 { print; next }
        {
            pop_af = \$16
            if (pop_af != "" && pop_af != ".") {
                if (pop_af + 0 > af_max) next
            }
            cons = \$11
            n = split(cons, parts, "&")
            keep = 0
            for (i=1; i<=n; i++) {
                if (parts[i] !~ "^(" exclude ")\$") { keep = 1; break }
            }
            if (!keep) next
            print
        }
        ' "\$all_tsv" > "\$cand_tsv"

    n_cand=\$(( \$(wc -l < "\$cand_tsv") - 1 ))
    echo "[stage 4] Somatic candidates: \$n_cand rows"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: ${params.vep_image}
        bcftools: \$(bcftools --version | head -1 | awk '{print \$NF}')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p vep_out
    touch vep_out/${meta.id}.all_annotated.tsv vep_out/${meta.id}.somatic_candidates.tsv
    echo '"${task.process}": stub' > versions.yml
    """
}
