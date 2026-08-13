#!/usr/bin/env python3
"""
plot_region_coverage.py

Per-region on-target coverage plot(s) from one or more *.region_coverage.tsv
files produced by qc_v6_sample.py.

Two modes, auto-selected by number of inputs:
  - ONE coverage TSV  -> horizontal bar chart of mean depth per region for that
    sample, regions ordered as in the file (genomic order), with the depth
    threshold drawn as a reference line.
  - MANY coverage TSVs -> heatmap of mean depth (regions x samples).

Optional --groups FILE: a 2-column TSV (region<TAB>group) used ONLY to color/
order bars by an external grouping you supply (e.g. panel version, pathway).
The script hardcodes no gene lists, findings, or panel assumptions; any grouping
is data you pass in. Regions with no group entry fall into "ungrouped".

Metric plotted is mean_depth by default; --metric pct_ge can plot the
percent-at-threshold column instead.

Usage:
  # single sample bar chart
  plot_region_coverage.py --out cov.png SAMPLE.region_coverage.tsv

  # single sample, colored by an external grouping
  plot_region_coverage.py --groups panel_groups.tsv --out cov.png SAMPLE.region_coverage.tsv

  # cohort heatmap
  plot_region_coverage.py --out cohort_cov.png S1.region_coverage.tsv S2.region_coverage.tsv ...
"""
import argparse
import os
import sys

import statistics

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_coverage(path, metric_col):
    """Return (sample_name, [(region, value), ...]) preserving file order."""
    rows = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            ri = header.index("region")
        except ValueError:
            raise SystemExit("no 'region' column in %s" % path)
        if metric_col not in header:
            raise SystemExit("no '%s' column in %s (have: %s)"
                             % (metric_col, path, ",".join(header)))
        mi = header.index(metric_col)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(ri, mi):
                continue
            rows.append((f[ri], float(f[mi])))
    # sample name from filename: strip .region_coverage.tsv
    base = os.path.basename(path)
    sample = base.replace(".region_coverage.tsv", "").replace("_region_coverage.tsv", "")
    return sample, rows


def load_groups(path):
    g = {}
    if not path:
        return g
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                g[f[0]] = f[1]
    return g


def bar_chart(sample, rows, groups, metric, threshold, out):
    regions = [r for r, _ in rows]
    vals = [v for _, v in rows]
    grp = [groups.get(r, "ungrouped") for r in regions]
    uniq = sorted(set(grp))
    # stable color per group
    cmap = plt.get_cmap("tab10")
    color_of = {g: cmap(i % 10) for i, g in enumerate(uniq)}
    colors = [color_of[g] for g in grp]

    n = len(regions)
    fig_h = max(4.0, 0.26 * n)
    plt.figure(figsize=(9, fig_h))
    y = np.arange(n)[::-1]  # top-to-bottom in file order
    plt.barh(y, vals, color=colors, edgecolor="white", linewidth=0.3)
    plt.yticks(y, regions, fontsize=8)
    # reference lines: sample median (over all plotted regions) and the
    # depth threshold. The median is computed over every region in the file
    # with no filtering, so on v5-era BAMs (where some panel regions were not
    # targeted) it reflects all bars as plotted; label states this plainly.
    handles_extra = []
    if vals:
        med = statistics.median(vals)
        ln_med = plt.axvline(med, color="#c0392b", linestyle="--", linewidth=1.5,
                             label="median (all regions) = %.2f x" % med)
        handles_extra.append(ln_med)
    if metric == "mean_depth" and threshold is not None:
        ln_thr = plt.axvline(threshold, color="#7f8c8d", linestyle=":", linewidth=1.4,
                             label="threshold = %g x" % threshold)
        handles_extra.append(ln_thr)
    xlabel = {"mean_depth": "mean depth (x)"}.get(metric, metric)
    plt.xlabel(xlabel)
    title = "On-target %s  (%s)" % (xlabel, sample)
    plt.title(title, fontsize=11)
    # legend: combine group swatches (if a real grouping was supplied) with the
    # reference lines, so neither overwrites the other.
    if groups and uniq != ["ungrouped"]:
        grp_handles = [plt.Rectangle((0, 0), 1, 1, color=color_of[g]) for g in uniq]
        all_handles = grp_handles + handles_extra
        all_labels = list(uniq) + [h.get_label() for h in handles_extra]
        plt.legend(all_handles, all_labels, frameon=False, fontsize=8,
                   loc="lower right")
    elif handles_extra:
        plt.legend(handles=handles_extra, frameon=False, fontsize=8,
                   loc="lower right")
    plt.tight_layout()
    plt.savefig(out, dpi=130)
    plt.close()
    sys.stderr.write("[ok] wrote %s (%d regions)\n" % (out, n))


def heatmap(samples_rows, metric, out):
    # union of regions in file order of the first sample, then any extras
    region_order = []
    seen = set()
    for _, rows in samples_rows:
        for r, _ in rows:
            if r not in seen:
                seen.add(r); region_order.append(r)
    samples = [s for s, _ in samples_rows]
    mat = np.full((len(region_order), len(samples)), np.nan)
    for j, (_, rows) in enumerate(samples_rows):
        d = dict(rows)
        for i, r in enumerate(region_order):
            if r in d:
                mat[i, j] = d[r]

    fig_h = max(4.0, 0.26 * len(region_order))
    fig_w = max(5.0, 1.1 * len(samples) + 3)
    plt.figure(figsize=(fig_w, fig_h))
    im = plt.imshow(mat, aspect="auto", cmap="viridis")
    plt.colorbar(im, label={"mean_depth": "mean depth (x)"}.get(metric, metric))
    plt.yticks(np.arange(len(region_order)), region_order, fontsize=8)
    plt.xticks(np.arange(len(samples)), samples, rotation=45, ha="right", fontsize=8)
    plt.title("On-target %s by region" % {"mean_depth": "mean depth"}.get(metric, metric),
              fontsize=11)
    plt.tight_layout()
    plt.savefig(out, dpi=130)
    plt.close()
    sys.stderr.write("[ok] wrote %s (%d regions x %d samples)\n"
                     % (out, len(region_order), len(samples)))


def main():
    ap = argparse.ArgumentParser(description="Plot per-region on-target coverage.")
    ap.add_argument("coverage_tsv", nargs="+", help="one or more *.region_coverage.tsv")
    ap.add_argument("--out", required=True, help="output PNG")
    ap.add_argument("--metric", default="mean_depth",
                    help="column to plot (default mean_depth)")
    ap.add_argument("--groups", default=None,
                    help="optional 2-col TSV: region<TAB>group (coloring only)")
    ap.add_argument("--threshold", type=float, default=15.0,
                    help="reference line for mean_depth bar chart (default 15)")
    args = ap.parse_args()

    groups = load_groups(args.groups)

    if len(args.coverage_tsv) == 1:
        sample, rows = load_coverage(args.coverage_tsv[0], args.metric)
        thr = args.threshold if args.metric == "mean_depth" else None
        bar_chart(sample, rows, groups, args.metric, thr, args.out)
    else:
        sr = [load_coverage(p, args.metric) for p in args.coverage_tsv]
        heatmap(sr, args.metric, args.out)


if __name__ == "__main__":
    main()
