#!/usr/bin/env python3
"""
qc_v6_sample.py

Per-sample on-target QC for the adaptive-WGS MM panel. For ONE T2T BAM and the
v6 chr-named panel BED, it produces:

  1. Per-region mean coverage for all panel regions (via mosdepth --by), plus
     the fraction of each region at >= a depth threshold (default 15x).
  2. On-target read-length statistics (median, N50, mean, count) computed over
     reads overlapping the panel, de-duplicated by read name so a read spanning
     two adjacent regions is counted once.
  3. On-target per-read quality statistics (median, mean) using the per-read
     mean of the BAM QUAL field. NOTE: this is the basecaller's per-read quality
     estimate (the ONT model at basecall time), not a mapping or consensus
     quality.
  4. Two PNG histograms: read length and per-read mean Q, for on-target reads.

Outputs (into --outdir):
  <sample>.region_coverage.tsv     one row per panel region
  <sample>.readlen_qscore.tsv      summary stats (one row, labelled metrics)
  <sample>.readlen_hist.png
  <sample>.qscore_hist.png

This script is variant-agnostic and finding-agnostic: it reports whatever the
data shows for whatever regions are in the BED. Region names come from the BED.

Requires (all in conda env awgs_sv): mosdepth, samtools on PATH; pysam, numpy,
matplotlib importable.

Usage:
  qc_v6_sample.py --bam <sample>.t2t.bam --bed panel_chr.bed \
      --sample <id> --outdir <dir> [--threshold 15] [--threads 4]
"""

import argparse
import gzip
import os
import statistics
import subprocess
import sys

import numpy as np
import pysam

import matplotlib
matplotlib.use("Agg")  # headless; no display on the server
import matplotlib.pyplot as plt


def eprint(*a):
    print(*a, file=sys.stderr)


def run_mosdepth(bam, bed, sample, outdir, threshold, threads):
    """Run mosdepth restricted to the panel; return path to regions.bed.gz.

    --by <bed>        : per-region summary (mean depth in col 5 of regions.bed.gz)
    --no-per-base     : we do not need base-resolution depth, much faster/smaller
    --thresholds T    : adds a thresholds.bed.gz with bases >= each cutoff
    --mapq 0          : keep MAPQ-0 reads (IGH-side reads are MAPQ-0 on T2T; the
                        project relies on these for translocation detection, so
                        excluding them would understate on-target coverage)
    """
    prefix = os.path.join(outdir, sample)
    cmd = [
        "mosdepth",
        "--by", bed,
        "--no-per-base",
        "--thresholds", str(threshold),
        "--mapq", "0",
        "--threads", str(threads),
        prefix,
        bam,
    ]
    eprint("[mosdepth]", " ".join(cmd))
    subprocess.run(cmd, check=True)
    regions = prefix + ".regions.bed.gz"
    thresholds = prefix + ".thresholds.bed.gz"
    if not os.path.exists(regions):
        raise SystemExit("mosdepth did not produce %s" % regions)
    return regions, thresholds


def parse_region_coverage(regions_gz, thresholds_gz, threshold):
    """Build per-region rows: name, chrom, start, end, span, mean_depth, pct_ge_thresh.

    mosdepth regions.bed.gz columns: chrom start end name mean
    thresholds.bed.gz columns:       chrom start end region <T>X  (bases >= T)
    """
    # mean depth keyed by (chrom,start,end)
    means = {}
    with gzip.open(regions_gz, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            chrom, start, end, name, mean = f[0], f[1], f[2], f[3], f[4]
            means[(chrom, start, end)] = (name, float(mean))

    # bases >= threshold keyed the same way (thresholds file present iff --thresholds)
    ge = {}
    if thresholds_gz and os.path.exists(thresholds_gz):
        with gzip.open(thresholds_gz, "rt") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            # last column is the count of bases >= threshold
            for line in fh:
                f = line.rstrip("\n").split("\t")
                chrom, start, end = f[0], f[1], f[2]
                bases_ge = int(f[-1])
                ge[(chrom, start, end)] = bases_ge

    rows = []
    for (chrom, start, end), (name, mean) in means.items():
        span = int(end) - int(start)
        bases_ge = ge.get((chrom, start, end), 0)
        pct = (100.0 * bases_ge / span) if span > 0 else 0.0
        rows.append({
            "region": name,
            "chrom": chrom,
            "start": int(start),
            "end": int(end),
            "span_bp": span,
            "mean_depth": round(mean, 2),
            "pct_ge_%dx" % threshold: round(pct, 1),
        })
    # keep BED order by (chrom index as encountered, start) -> just sort by start within file order
    rows.sort(key=lambda r: (r["chrom"], r["start"]))
    return rows


def ontarget_read_stats(bam, bed):
    """Collect read length and per-read mean Q for reads overlapping the panel.

    De-dup by read name: a read overlapping two adjacent regions is fetched
    once per region by pysam, so we track seen names. Uses each region from the
    BED to fetch, which is far faster than scanning the whole BAM.
    """
    # load regions
    regions = []
    with open(bed) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            regions.append((f[0], int(f[1]), int(f[2])))

    lengths = []
    qmeans = []
    seen = set()
    bf = pysam.AlignmentFile(bam, "rb")
    bam_refs = set(bf.references)
    for chrom, start, end in regions:
        if chrom not in bam_refs:
            eprint("[warn] %s not in BAM header; skipping region" % chrom)
            continue
        for r in bf.fetch(chrom, start, end):
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            name = r.query_name
            if name in seen:
                continue
            seen.add(name)
            # query_length is the read length as aligned-record sequence length;
            # prefer infer_read_length() to include hard-clipped bases (full molecule)
            rl = r.infer_read_length() or r.query_length
            if rl:
                lengths.append(int(rl))
            q = r.query_qualities
            if q is not None and len(q) > 0:
                qmeans.append(float(np.mean(q)))
    bf.close()
    return lengths, qmeans


def n50(values):
    if not values:
        return 0
    s = sorted(values, reverse=True)
    half = sum(s) / 2.0
    run = 0
    for v in s:
        run += v
        if run >= half:
            return v
    return s[-1]


def write_region_tsv(rows, path, threshold):
    cols = ["region", "chrom", "start", "end", "span_bp", "mean_depth",
            "pct_ge_%dx" % threshold]
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")


def write_summary_tsv(sample, lengths, qmeans, path, threshold, rows):
    n_reads = len(lengths)
    med_len = int(statistics.median(lengths)) if lengths else 0
    mean_len = int(statistics.mean(lengths)) if lengths else 0
    read_n50 = n50(lengths)
    med_q = round(statistics.median(qmeans), 2) if qmeans else 0.0
    mean_q = round(statistics.mean(qmeans), 2) if qmeans else 0.0
    # panel-level coverage rollups
    depths = [r["mean_depth"] for r in rows]
    panel_mean_depth = round(statistics.mean(depths), 2) if depths else 0.0
    panel_median_depth = round(statistics.median(depths), 2) if depths else 0.0
    n_regions = len(rows)
    n_below = sum(1 for r in rows if r["mean_depth"] < threshold)

    fields = [
        ("sample", sample),
        ("ontarget_reads", n_reads),
        ("median_read_length_bp", med_len),
        ("mean_read_length_bp", mean_len),
        ("read_length_N50_bp", read_n50),
        ("median_per_read_Q_basecaller", med_q),
        ("mean_per_read_Q_basecaller", mean_q),
        ("n_panel_regions", n_regions),
        ("panel_mean_region_depth", panel_mean_depth),
        ("panel_median_region_depth", panel_median_depth),
        ("n_regions_below_%dx" % threshold, n_below),
    ]
    with open(path, "w") as fh:
        fh.write("\t".join(k for k, _ in fields) + "\n")
        fh.write("\t".join(str(v) for _, v in fields) + "\n")
    return dict(fields)


def plot_hist(values, title, xlabel, path, sample, median_val=None, bins=60):
    plt.figure(figsize=(8, 4.5))
    if values:
        arr = np.array(values, dtype=float)
        plt.hist(arr, bins=bins, color="#3b6ea5", edgecolor="white", linewidth=0.3)
        if median_val is not None:
            plt.axvline(median_val, color="#c0392b", linestyle="--", linewidth=1.5,
                        label="median = %s" % (int(median_val) if median_val >= 10 else round(median_val, 2)))
            plt.legend(frameon=False)
    else:
        plt.text(0.5, 0.5, "no on-target reads", ha="center", va="center",
                 transform=plt.gca().transAxes)
    plt.title("%s  (%s, on-target)" % (title, sample))
    plt.xlabel(xlabel)
    plt.ylabel("read count")
    plt.tight_layout()
    plt.savefig(path, dpi=120)
    plt.close()


def main():
    ap = argparse.ArgumentParser(description="Per-sample on-target QC for the v6 MM panel.")
    ap.add_argument("--bam", required=True, help="T2T chr-named BAM (indexed)")
    ap.add_argument("--bed", required=True, help="v6 chr-named panel BED")
    ap.add_argument("--sample", required=True, help="sample id (output prefix)")
    ap.add_argument("--outdir", required=True, help="output directory")
    ap.add_argument("--threshold", type=int, default=15,
                    help="depth threshold for pct-covered and region flag (default 15)")
    ap.add_argument("--threads", type=int, default=4, help="mosdepth threads")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # 1. coverage
    regions_gz, thresholds_gz = run_mosdepth(
        args.bam, args.bed, args.sample, args.outdir, args.threshold, args.threads)
    rows = parse_region_coverage(regions_gz, thresholds_gz, args.threshold)
    region_tsv = os.path.join(args.outdir, "%s.region_coverage.tsv" % args.sample)
    write_region_tsv(rows, region_tsv, args.threshold)
    eprint("[ok] wrote %s (%d regions)" % (region_tsv, len(rows)))

    # 2 + 3. read length & quality (on-target)
    lengths, qmeans = ontarget_read_stats(args.bam, args.bed)
    summary_tsv = os.path.join(args.outdir, "%s.readlen_qscore.tsv" % args.sample)
    summary = write_summary_tsv(args.sample, lengths, qmeans, summary_tsv, args.threshold, rows)
    eprint("[ok] wrote %s" % summary_tsv)

    # 4. histograms
    med_len = statistics.median(lengths) if lengths else None
    med_q = statistics.median(qmeans) if qmeans else None
    len_png = os.path.join(args.outdir, "%s.readlen_hist.png" % args.sample)
    q_png = os.path.join(args.outdir, "%s.qscore_hist.png" % args.sample)
    plot_hist(lengths, "On-target read length", "read length (bp)", len_png, args.sample, med_len)
    plot_hist(qmeans, "On-target per-read mean Q (basecaller)", "per-read mean Q", q_png, args.sample, med_q)
    eprint("[ok] wrote %s , %s" % (len_png, q_png))

    # console summary
    eprint("\n=== %s on-target QC ===" % args.sample)
    for k, v in summary.items():
        eprint("  %-32s %s" % (k, v))


if __name__ == "__main__":
    main()
