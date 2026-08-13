#!/usr/bin/env python3
"""
sample_qc.py
---------------------------------------------------------------------------
Generate per-sample QC for one ONT adaptive-WGS BAM and emit it as an HTML
fragment suitable for igv-reports `create_report --header`, plus a metrics
JSON. Nothing is collected from pre-existing artifacts (there are none in the
results tree); QC is computed fresh from tools already in awgs_sv:

  - coverage         mosdepth (global mean + on-target mean over the panel BED
                     + per-region depth), parsed from mosdepth's summary/regions
  - read length      pysam pass (sampled), N50 / median / max + histogram
  - per-read mean Q   pysam pass (sampled), median Q + histogram

NOTE on "insert size": ONT reads are single-molecule; there is no paired-end
insert size. The ONT-appropriate analogs are the read-length distribution
(with N50) and the per-read mean quality distribution, which is what this
produces.

Plots are hand-rendered as inline SVG (no matplotlib / NanoPlot dependency),
so the fragment embeds directly in the report and stays self-contained.

Dependencies: Python stdlib + pysam (in awgs_sv). mosdepth + samtools on PATH.
No sample-specific findings are hardcoded; only the panel BED (reference
feature) and the BAM are read.
---------------------------------------------------------------------------
"""

import argparse
import gzip
import html
import json
import os
import subprocess
import sys

import pysam


def eprint(*a):
    print(*a, file=sys.stderr)


def fail(msg, code=2):
    eprint("ERROR: " + msg)
    sys.exit(code)


# ---------------------------------------------------------------------------
# Coverage via mosdepth
# ---------------------------------------------------------------------------

def run_mosdepth(bam, panel_bed, out_prefix, threads):
    """Run mosdepth (fast mode, panel as regions). Returns the prefix used."""
    md_prefix = out_prefix  # mosdepth writes <prefix>.mosdepth.summary.txt etc.
    cmd = ["mosdepth", "--no-per-base", "-t", str(threads)]
    if panel_bed:
        cmd += ["--by", panel_bed]
    cmd += [md_prefix, bam]
    eprint("[run] " + " ".join(cmd))
    subprocess.run(cmd, check=True)
    return md_prefix


def parse_mosdepth_summary(prefix):
    """From <prefix>.mosdepth.summary.txt return (global_mean, ontarget_mean).
    The 'total' row is genome-wide; 'total_region' appears when --by is used."""
    path = prefix + ".mosdepth.summary.txt"
    global_mean, ontarget_mean = None, None
    if not os.path.isfile(path):
        return None, None
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            mean_idx = header.index("mean")
        except ValueError:
            mean_idx = 3
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if not f:
                continue
            if f[0] == "total":
                global_mean = float(f[mean_idx])
            elif f[0] == "total_region":
                ontarget_mean = float(f[mean_idx])
    return global_mean, ontarget_mean


def parse_mosdepth_regions(prefix, top_n=40):
    """From <prefix>.regions.bed.gz return [(label, mean_depth), ...] for the
    per-panel-region depths (label taken from BED col 4 if present)."""
    path = prefix + ".regions.bed.gz"
    rows = []
    if not os.path.isfile(path):
        return rows
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            depth = f[-1]
            label = f[3] if len(f) >= 5 else f"{f[0]}:{f[1]}-{f[2]}"
            try:
                rows.append((label, float(depth)))
            except ValueError:
                continue
    # collapse duplicate region names (a gene may span multiple BED rows)
    agg = {}
    for label, d in rows:
        agg.setdefault(label, []).append(d)
    merged = [(k, sum(v) / len(v)) for k, v in agg.items()]
    merged.sort(key=lambda kv: kv[1], reverse=True)
    return merged[:top_n]


# ---------------------------------------------------------------------------
# Read length + per-read mean Q via pysam (sampled for speed)
# ---------------------------------------------------------------------------

def total_mapped(bam):
    try:
        total = 0
        for line in pysam.idxstats(bam).splitlines():
            parts = line.split("\t")
            if len(parts) >= 3:
                total += int(parts[2])
        return total
    except Exception:
        return 0


def collect_read_metrics(bam, max_reads):
    """Sample primary mapped reads evenly; return (lengths, mean_quals,
    n_total_mapped, n_sampled)."""
    n_total = total_mapped(bam)
    stride = max(1, n_total // max_reads) if n_total else 1
    lengths, quals = [], []
    seen = 0
    with pysam.AlignmentFile(bam, "rb") as af:
        for read in af.fetch(until_eof=True):
            if read.is_secondary or read.is_supplementary or read.is_unmapped:
                continue
            seen += 1
            if seen % stride != 0:
                continue
            ql = read.query_length or (len(read.query_sequence)
                                       if read.query_sequence else 0)
            if ql:
                lengths.append(ql)
                q = read.query_qualities
                if q is not None and len(q) > 0:
                    quals.append(sum(q) / len(q))
            if len(lengths) >= max_reads:
                break
    return lengths, quals, n_total, len(lengths)


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


def histogram(values, bin_width, max_edge=None):
    """Return (edges, counts). Last bin captures the overflow."""
    if not values:
        return [], []
    hi = max(values)
    if max_edge is None:
        max_edge = ((int(hi) // bin_width) + 1) * bin_width
    nbins = max(1, int(max_edge // bin_width))
    counts = [0] * (nbins + 1)  # +1 overflow
    for v in values:
        idx = int(v // bin_width)
        if idx >= nbins:
            idx = nbins
        counts[idx] += 1
    edges = [i * bin_width for i in range(nbins + 1)]
    return edges, counts


# ---------------------------------------------------------------------------
# Inline-SVG rendering (no external libs)
# ---------------------------------------------------------------------------

def esc(x):
    return html.escape("" if x is None else str(x))


def bar_svg(labels, values, title, unit="", width=560, height=190,
            color="#3a6ea5"):
    """Hand-rendered SVG bar chart."""
    if not values:
        return f'<p class="qc-empty">{esc(title)}: no data.</p>'
    pad_l, pad_b, pad_t, pad_r = 44, 34, 24, 8
    plot_w = width - pad_l - pad_r
    plot_h = height - pad_t - pad_b
    vmax = max(values) or 1
    n = len(values)
    bw = plot_w / n
    bars = []
    for i, v in enumerate(values):
        bh = (v / vmax) * plot_h
        x = pad_l + i * bw
        y = pad_t + (plot_h - bh)
        bars.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{max(0.5, bw-1):.1f}" '
                    f'height="{bh:.1f}" fill="{color}"><title>{esc(labels[i])}: '
                    f'{esc(v)}</title></rect>')
    # y axis ticks (0, max)
    yt = (f'<text x="{pad_l-6}" y="{pad_t+8}" text-anchor="end" '
          f'class="qc-tick">{esc(vmax)}</text>'
          f'<text x="{pad_l-6}" y="{pad_t+plot_h}" text-anchor="end" '
          f'class="qc-tick">0</text>')
    # x labels: first, middle, last (avoid clutter)
    xl = []
    for i in (0, n // 2, n - 1):
        if 0 <= i < n:
            x = pad_l + i * bw + bw / 2
            xl.append(f'<text x="{x:.1f}" y="{height-pad_b+14}" '
                      f'text-anchor="middle" class="qc-tick">{esc(labels[i])}</text>')
    axis = (f'<line x1="{pad_l}" y1="{pad_t}" x2="{pad_l}" y2="{pad_t+plot_h}" '
            f'class="qc-axis"/><line x1="{pad_l}" y1="{pad_t+plot_h}" '
            f'x2="{pad_l+plot_w}" y2="{pad_t+plot_h}" class="qc-axis"/>')
    return (f'<svg viewBox="0 0 {width} {height}" class="qc-svg" '
            f'role="img" aria-label="{esc(title)}">'
            f'<text x="{pad_l}" y="14" class="qc-title">{esc(title)}'
            f'{(" ("+esc(unit)+")") if unit else ""}</text>'
            f'{axis}{"".join(bars)}{yt}{"".join(xl)}</svg>')


def kb_label(bp):
    return f"{bp/1000:.0f}k"


def build_fragment(sample, metrics, rl_svg, q_svg, cov_svg):
    """Assemble the QC HTML fragment (a div) for create_report --header."""
    m = metrics
    def cell(label, value):
        return f'<div class="qc-metric"><span>{esc(label)}</span><b>{esc(value)}</b></div>'
    cov_global = "NA" if m["global_mean"] is None else f'{m["global_mean"]:.2f}x'
    cov_on = "NA" if m["ontarget_mean"] is None else f'{m["ontarget_mean"]:.2f}x'
    enr = "NA"
    if m["global_mean"] and m["ontarget_mean"] and m["global_mean"] > 0:
        enr = f'{m["ontarget_mean"]/m["global_mean"]:.1f}x'
    cards = "".join([
        cell("Reads (mapped)", f'{m["n_total_mapped"]:,}'),
        cell("Sampled", f'{m["n_sampled"]:,}'),
        cell("Read N50", f'{m["read_n50"]:,} bp'),
        cell("Median length", f'{m["read_median"]:,} bp'),
        cell("Max length", f'{m["read_max"]:,} bp'),
        cell("Median read Q", f'{m["q_median"]:.1f}'),
        cell("Global coverage", cov_global),
        cell("On-target coverage", cov_on),
        cell("Enrichment", enr),
    ])
    style = """
<style>
.qc-wrap{font:13px/1.45 -apple-system,Segoe UI,Roboto,Arial,sans-serif;
  border:1px solid #e3e3e3;border-radius:8px;padding:14px 16px;margin:0 0 14px;}
.qc-wrap h3{margin:0 0 8px;font-size:15px;}
.qc-metrics{display:flex;flex-wrap:wrap;gap:8px;margin-bottom:10px;}
.qc-metric{background:#f2f4f5;border-radius:6px;padding:4px 10px;font-size:12px;}
.qc-metric span{color:#666;margin-right:6px;}
.qc-plots{display:flex;flex-wrap:wrap;gap:14px;}
.qc-svg{width:560px;max-width:100%;height:auto;background:#fafafa;
  border:1px solid #eee;border-radius:6px;}
.qc-title{font-size:12px;font-weight:600;fill:#244;}
.qc-tick{font-size:10px;fill:#666;}
.qc-axis{stroke:#bbb;stroke-width:1;}
.qc-empty{color:#888;font-style:italic;}
.qc-note{color:#888;font-size:11px;margin-top:8px;}
</style>"""
    return (f'{style}<div class="qc-wrap"><h3>QC &mdash; {esc(sample)}</h3>'
            f'<div class="qc-metrics">{cards}</div>'
            f'<div class="qc-plots">{rl_svg}{q_svg}{cov_svg}</div>'
            f'<p class="qc-note">ONT single-molecule reads: read-length (N50) and '
            f'per-read mean Q stand in for paired-end insert size. Coverage from '
            f'mosdepth; distributions from a sample of {m["n_sampled"]:,} reads.</p>'
            f'</div>')


def main():
    ap = argparse.ArgumentParser(description="Per-sample ONT QC -> HTML fragment + JSON.")
    ap.add_argument("--bam", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--panel-bed", default="", help="Panel BED for on-target coverage")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--max-reads", type=int, default=200000,
                    help="Reads sampled for length/Q histograms (default 200000)")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--rl-bin", type=int, default=2000, help="Read-length bin bp")
    ap.add_argument("--skip-mosdepth", action="store_true",
                    help="Reuse existing <out-prefix>.mosdepth.* instead of rerunning")
    args = ap.parse_args()

    if not os.path.isfile(args.bam):
        fail(f"BAM not found: {args.bam}")
    if not (os.path.isfile(args.bam + ".bai") or os.path.isfile(args.bam + ".csi")
            or os.path.isfile(args.bam[:-4] + ".bai")):
        eprint("[warn] no BAM index found; idxstats sampling may be slow")

    # 1) coverage
    if not args.skip_mosdepth:
        try:
            run_mosdepth(args.bam, args.panel_bed or None, args.out_prefix, args.threads)
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            eprint(f"[warn] mosdepth failed ({e}); coverage will be reported as NA")
    global_mean, ontarget_mean = parse_mosdepth_summary(args.out_prefix)
    per_region = parse_mosdepth_regions(args.out_prefix) if args.panel_bed else []

    # 2) read length + Q
    lengths, quals, n_total, n_sampled = collect_read_metrics(args.bam, args.max_reads)
    read_median = sorted(lengths)[len(lengths) // 2] if lengths else 0
    q_median = sorted(quals)[len(quals) // 2] if quals else 0.0

    metrics = {
        "sample": args.sample,
        "n_total_mapped": n_total,
        "n_sampled": n_sampled,
        "read_n50": n50(lengths),
        "read_median": read_median,
        "read_max": max(lengths) if lengths else 0,
        "q_median": q_median,
        "global_mean": global_mean,
        "ontarget_mean": ontarget_mean,
        "per_region_depth": per_region,
    }

    # 3) plots
    rl_edges, rl_counts = histogram(lengths, args.rl_bin)
    rl_labels = [kb_label(e) for e in rl_edges]
    rl_svg = bar_svg(rl_labels, rl_counts, "Read length", unit="bp", color="#3a6ea5")

    q_edges, q_counts = histogram(quals, 2, max_edge=40)
    q_labels = [str(e) for e in q_edges]
    q_svg = bar_svg(q_labels, q_counts, "Per-read mean Q", color="#6a8d3a")

    if per_region:
        cov_labels = [r[0] for r in per_region]
        cov_vals = [round(r[1], 1) for r in per_region]
        cov_svg = bar_svg(cov_labels, cov_vals, "On-target depth by region",
                          unit="x", color="#9a5a3a")
    else:
        cov_svg = bar_svg(
            ["global", "on-target"],
            [round(global_mean or 0, 1), round(ontarget_mean or 0, 1)],
            "Coverage", unit="x", color="#9a5a3a")

    # 4) write outputs
    frag = build_fragment(args.sample, metrics, rl_svg, q_svg, cov_svg)
    html_path = args.out_prefix + ".qc.html"
    json_path = args.out_prefix + ".qc.json"
    with open(html_path, "w") as fh:
        fh.write(frag)
    with open(json_path, "w") as fh:
        json.dump(metrics, fh, indent=2)
    eprint(f"[done] QC fragment: {html_path}")
    eprint(f"[done] QC metrics : {json_path}")
    print(html_path)


if __name__ == "__main__":
    main()
