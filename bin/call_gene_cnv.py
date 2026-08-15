#!/usr/bin/env python3
"""
call_gene_cnv.py
================

Gene-level and intragenic copy number from panel-enriched long-read data.

WHY NOT ichorCNA
----------------
ichorCNA works from off-target reads in 1 Mb bins and answers a chromosome-arm
question well. It cannot see a 50 kb intragenic deletion, and on a flat genome
its ploidy/normal-fraction fit becomes degenerate. This answers the other
question: what is the copy number of each panel gene, and is any part of a gene
deleted. The two are complementary and both are reported.

METHOD
------
Depth is binned within each panel region, then normalised twice.

1. Against the sample's own autosomal on-target median, which cancels
   library yield and flow-cell differences.
2. Against the per-bin median across a reference set of samples, when one is
   supplied. This is the step that matters: enrichment efficiency, GC content
   and mappability vary hugely between bins but are reproducible across
   samples on the same panel, so dividing by the cohort median removes bias
   that no within-sample correction can reach. It is the same principle as a
   CNVkit pooled reference.

Without a reference set the caller still runs, using within-sample
normalisation only. Gene-level calls remain usable; intragenic calls are
weaker, because a bin that is intrinsically poorly enriched looks the same as a
deleted one.

PLOIDY
------
chrX bins are scaled by the observed X:autosome ratio rather than assumed
diploid, so hemizygous regions in a male sample are not reported as losses.
The ratio is measured, not declared.

Dependencies: pysam, numpy, matplotlib (all in the awgs_sv env).
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

try:
    import pysam
except ImportError:
    sys.stderr.write("ERROR: pysam is required. Activate the awgs_sv env.\n")
    sys.exit(1)

__version__ = "0.1.0"

SUFFIX = {"LOCUS", "CLUSTER", "ENHANCER", "REGION", "INTERVAL"}


# -----------------------------------------------------------------------------
# Inputs
# -----------------------------------------------------------------------------
def load_panel(path: Path) -> List[Tuple[str, int, int, str]]:
    out = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith(("#", "track")):
            continue
        f = line.split("\t")
        if len(f) >= 4:
            out.append((f[0], int(f[1]), int(f[2]), f[3]))
    return out


def bin_depth(bam, chrom: str, start: int, end: int, binsize: int,
              min_mapq: int) -> Tuple[np.ndarray, np.ndarray]:
    """Mean depth per bin. Returns (bin_starts, depths).

    Counted from alignment spans rather than pileup: at this bin size the
    difference is immaterial and the span walk is an order of magnitude
    faster over a 40 Mb panel.
    """
    n = max(1, (end - start) // binsize)
    edges = start + np.arange(n + 1) * binsize
    cov = np.zeros(n, dtype=np.float64)
    try:
        it = bam.fetch(chrom, start, end)
    except ValueError:
        return edges[:-1], cov
    for a in it:
        if a.is_unmapped or a.is_secondary or a.is_supplementary:
            continue
        if a.mapping_quality < min_mapq:
            continue
        s = max(a.reference_start, start)
        e = min(a.reference_end or s, end)
        if e <= s:
            continue
        lo = (s - start) // binsize
        hi = min(n - 1, (e - 1 - start) // binsize)
        for b in range(lo, hi + 1):
            bs = start + b * binsize
            be = bs + binsize
            cov[b] += (min(e, be) - max(s, bs))
    return edges[:-1], cov / binsize


def sample_bins(bam_path: Path, panel, binsize: int, min_mapq: int) -> Dict:
    """All per-bin depths for one sample, keyed by (chrom, bin_start)."""
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    bins: Dict[Tuple[str, int], float] = {}
    region_of: Dict[Tuple[str, int], str] = {}
    for chrom, start, end, name in panel:
        if chrom not in bam.references:
            continue
        starts, depths = bin_depth(bam, chrom, start, end, binsize, min_mapq)
        for s, d in zip(starts, depths):
            bins[(chrom, int(s))] = float(d)
            region_of[(chrom, int(s))] = name
    bam.close()
    return {"bins": bins, "region": region_of}


# -----------------------------------------------------------------------------
# Normalisation
# -----------------------------------------------------------------------------
def autosomal_median(bins: Dict[Tuple[str, int], float]) -> float:
    vals = [v for (c, _), v in bins.items()
            if c not in ("chrX", "chrY") and v > 0]
    return float(np.median(vals)) if vals else 0.0


def x_ratio(bins: Dict[Tuple[str, int], float], auto_med: float) -> float:
    """Observed chrX depth relative to autosomes, so hemizygous regions in a
    male sample are not reported as losses."""
    xs = [v for (c, _), v in bins.items() if c == "chrX" and v > 0]
    if not xs or auto_med <= 0:
        return 1.0
    return float(np.median(xs)) / auto_med


def build_reference(ref_samples: List[Dict], max_disagree: float
                    ) -> Tuple[Dict[Tuple[str, int], float], Dict[str, int]]:
    """Per-bin baseline from the reference set, with disagreeing bins dropped.

    A reference is only a baseline where its members agree. If one reference
    sample carries a copy-number change at a locus, the summary of that locus
    is pulled toward it and every other sample reads as the opposite change:
    a single trisomic member makes diploid samples look deleted by
    log2(2/2.5) = -0.32, which is large enough to call.

    The median offers no protection at small n, since the median of two values
    is their mean. Instead a bin is used only when its reference values agree
    within max_disagree of each other, measured as max/min. Discordant bins
    are dropped rather than averaged, so the locus goes uncalled instead of
    miscalled.
    """
    stacked: Dict[Tuple[str, int], List[float]] = defaultdict(list)
    for s in ref_samples:
        med = autosomal_median(s["bins"])
        if med <= 0:
            continue
        for k, v in s["bins"].items():
            stacked[k].append(v / med)

    ref, stats = {}, {"total": 0, "kept": 0, "discordant": 0, "too_few": 0}
    for k, vals in stacked.items():
        stats["total"] += 1
        vals = [v for v in vals if v > 0]
        if len(vals) < 2:
            stats["too_few"] += 1
            continue
        if max_disagree > 0 and max(vals) / min(vals) > max_disagree:
            stats["discordant"] += 1
            continue
        ref[k] = float(np.median(vals))
        stats["kept"] += 1
    return ref, stats


# -----------------------------------------------------------------------------
# Calling
# -----------------------------------------------------------------------------
def mad(a: np.ndarray) -> float:
    """Median absolute deviation, scaled to a standard-deviation equivalent."""
    a = a[~np.isnan(a)]
    if a.size == 0:
        return 0.0
    return 1.4826 * float(np.median(np.abs(a - np.median(a))))


def segment(log2: np.ndarray, positions: np.ndarray, threshold: float,
            min_bins: int) -> List[Tuple[int, int, float, int]]:
    """Contiguous runs of bins on the same side of the threshold.

    A plain run-length scan rather than circular binary segmentation. The
    events this needs to find are focal and sharp-edged, IKZF1 exon 4-7 or a
    CDKN2A homozygous loss, and the panel gives at most a few hundred bins per
    gene, so the extra machinery buys nothing.
    """
    # Deviation is measured from the gene's own baseline, not from zero. A
    # gene sitting on a trisomic chromosome has every bin at log2 0.5; the
    # question for an intragenic call is whether part of it departs from the
    # rest of it.
    baseline = float(np.nanmedian(log2)) if len(log2) else 0.0
    log2 = log2 - baseline

    out = []
    n = len(log2)
    i = 0
    while i < n:
        if np.isnan(log2[i]) or abs(log2[i]) < threshold:
            i += 1
            continue
        sign = np.sign(log2[i])
        j = i
        while j < n and not np.isnan(log2[j]) and np.sign(log2[j]) == sign \
                and abs(log2[j]) >= threshold:
            j += 1
        if j - i >= min_bins:
            seg = log2[i:j]
            out.append((int(positions[i]), int(positions[j - 1]),
                        float(np.median(seg)) + baseline, j - i))
        i = j
    return out


def cn_from_log2(l2: float, purity: float) -> float:
    """Absolute copy number from a log2 ratio at a given tumour fraction.

    ratio = (purity * CN + 2 * (1 - purity)) / 2, solved for CN. At purity 1
    this is the textbook 2 * 2^log2.
    """
    ratio = 2.0 ** l2
    if purity <= 0:
        return 2.0 * ratio
    return (2.0 * ratio - 2.0 * (1.0 - purity)) / purity


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True, type=Path)
    ap.add_argument("--panel-bed", required=True, type=Path)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out-prefix", required=True, type=Path)
    ap.add_argument("--reference-bam", action="append", default=[], type=Path,
                    help="BAM of another sample on the same panel, used to "
                         "build a per-bin reference. Repeatable. Two or more "
                         "makes intragenic calls far more reliable; without "
                         "any, only within-sample normalisation is applied.")
    ap.add_argument("--binsize", type=int, default=1000,
                    help="Bin width within panel regions [1000]. A 50 kb "
                         "intragenic deletion spans 50 bins at this setting.")
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--max-reference-disagreement", type=float, default=0.0,
                    help="Drop a bin when its reference values differ by more "
                         "than this ratio, max over min. 0 disables [0]. Off "
                         "by default: adaptive sampling depth is not "
                         "reproducible bin to bin the way capture is, so at "
                         "small n this rejects ordinary variance and the "
                         "surviving bins are a biased subset. It also cannot "
                         "catch a whole-chromosome change, whose bins agree "
                         "with each other at the wrong level.")
    ap.add_argument("--homdel-log2", type=float, default=3.0,
                    help="log2 magnitude recorded for a homozygous loss [3.0]. "
                         "The true value is unbounded; this is a floor for "
                         "reporting.")
    ap.add_argument("--min-covered-fraction", type=float, default=0.25,
                    help="A bin counts as covered at this fraction of the "
                         "sample's autosomal median depth [0.25].")
    ap.add_argument("--min-covered-bins", type=float, default=0.5,
                    help="A region needs this fraction of its bins covered "
                         "before a call is attempted [0.5]. Below it the "
                         "region is reported as uninformative.")
    ap.add_argument("--min-chrom-regions", type=int, default=3,
                    help="Panel regions a chromosome needs before its own "
                         "median is used as the baseline for its genes [3]. "
                         "Below this the panel baseline is used.")
    ap.add_argument("--min-bins", type=int, default=30,
                    help="Shortest reportable segment, in bins [30]. At the "
                         "default bin size that is 30 kb, below which real "
                         "focal events are rare and depth noise is not.")
    ap.add_argument("--seg-sd", type=float, default=4.0,
                    help="Segment threshold in units of the sample's own "
                         "bin-level MAD [4.0]. Derived per sample rather than "
                         "fixed, since bin scatter depends on depth.")
    ap.add_argument("--threshold", type=float, default=0.3,
                    help="Absolute log2 ratio for a bin to join a segment "
                         "[0.4]. A clonal heterozygous loss is -1.0; 0.4 "
                         "catches subclonal events without chasing noise.")
    ap.add_argument("--purity", type=float, default=1.0,
                    help="Tumour fraction for absolute CN [1.0]. Blast count "
                         "or the on-target estimate is better than ichorCNA's "
                         "on a near-flat genome.")
    ap.add_argument("--plot-genes", default="",
                    help="Comma-separated genes to plot, or 'called' for "
                         "every gene with a segment, or 'all'.")
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args = ap.parse_args()

    panel = load_panel(args.panel_bed)
    if not panel:
        sys.stderr.write(f"ERROR: no regions in {args.panel_bed}\n")
        return 1

    sys.stderr.write(f"binning {len(panel)} regions at {args.binsize} bp\n")
    smp = sample_bins(args.bam, panel, args.binsize, args.min_mapq)
    auto_med = autosomal_median(smp["bins"])
    if auto_med <= 0:
        sys.stderr.write("ERROR: zero autosomal on-target depth.\n")
        return 1
    xr = x_ratio(smp["bins"], auto_med)
    sys.stderr.write(f"autosomal median {auto_med:.1f}x, X:A {xr:.2f} "
                     f"({'one X' if xr < 0.75 else 'two X'})\n")

    ref, ref_xr = {}, 1.0
    if args.reference_bam:
        refs = []
        for rp in args.reference_bam:
            sys.stderr.write(f"  reference: {rp.name}\n")
            refs.append(sample_bins(rp, panel, args.binsize, args.min_mapq))
        ref, rstats = build_reference(refs, args.max_reference_disagreement)
        ref_xr = float(np.median([x_ratio(r["bins"], autosomal_median(r["bins"]))
                                  for r in refs])) or 1.0
        sys.stderr.write(
            f"reference built from {len(refs)} sample(s), X:A {ref_xr:.2f}: "
            f"{rstats['kept']:,} bins kept, {rstats['discordant']:,} dropped "
            f"for disagreement, {rstats['too_few']:,} for too few values\n")
        if len(refs) < 4:
            sys.stderr.write(
                f"  NOTE: {len(refs)} reference sample(s). A copy-number change "
                f"in one member biases every other sample at that locus in the "
                f"opposite direction. Bins where the references disagree are "
                f"dropped, but with few samples a shared change cannot be "
                f"distinguished from a shared baseline. Treat calls at loci "
                f"altered in a reference member with caution.\n")
    else:
        sys.stderr.write(
            "no reference samples: within-sample normalisation only. "
            "Gene-level calls stand; intragenic calls are weaker, because a "
            "poorly enriched bin and a deleted bin look alike.\n")

    # Per-bin log2, grouped by region.
    by_region: Dict[str, List[Tuple[str, int, float]]] = defaultdict(list)
    for key, depth in smp["bins"].items():
        chrom, pos = key
        # The X correction is the ratio of this sample's X dosage to the
        # reference's. Comparing a two-X sample against a reference containing
        # a one-X sample otherwise reports every X-linked gene as a gain of
        # exactly the dosage difference.
        # Two values per bin, kept apart deliberately.
        #
        # l2_self depends on no other sample, so a copy-number change in a
        # reference member cannot reach it. Copy-number statements, meaning
        # the chromosome and gene calls, are made from this.
        #
        # l2_ref additionally divides by the cohort median for the bin, which
        # cancels enrichment, GC and mappability bias that no within-sample
        # correction can reach, at the cost of importing whatever the
        # reference members carry. Used only for intragenic segments, where
        # that bias is the limiting factor and the gene's own baseline absorbs
        # any offset the reference introduces.
        exp_self = auto_med * (xr if chrom == "chrX" else 1.0)
        if exp_self <= 0:
            continue
        l2_self = np.log2(max(depth, 0.01) / exp_self)

        l2_ref = l2_self
        if ref:
            r = ref.get(key)
            if r is not None and r > 0:
                x_adj = (xr / ref_xr) if chrom == "chrX" and ref_xr > 0 else 1.0
                exp_ref = auto_med * r * x_adj
                if exp_ref > 0:
                    l2_ref = np.log2(max(depth, 0.01) / exp_ref)
        by_region[smp["region"][key]].append((chrom, pos, l2_self, l2_ref))

    # Chromosome baselines. A gene is focally altered when it departs from its
    # own chromosome, not from the panel: on a trisomic chromosome every gene
    # is gained, which is one finding and not six. Requires enough regions for
    # the median to be a baseline rather than an average of the thing being
    # measured; below that the panel baseline is used, since a chromosome
    # carrying two panel genes cannot describe itself.
    per_chrom: Dict[str, List[float]] = defaultdict(list)
    for name, vals in by_region.items():
        if not vals:
            continue
        c = vals[0][0]
        d = np.array([smp["bins"][(c, int(v[1]))] for v in vals])
        if float(np.mean(d >= auto_med * args.min_covered_fraction)) < args.min_covered_bins:
            continue
        per_chrom[c].append(float(np.median([v[2] for v in vals])))
    chrom_base = {c: float(np.median(v)) for c, v in per_chrom.items()
                  if len(v) >= args.min_chrom_regions}

    chrom_rows = []
    for c in sorted(per_chrom, key=lambda x: (len(x), x)):
        med = float(np.median(per_chrom[c]))
        chrom_rows.append({
            "sample": args.sample, "chrom": c, "n_regions": len(per_chrom[c]),
            "median_log2": f"{med:.3f}",
            "copy_number": f"{cn_from_log2(med, args.purity):.2f}",
            "call": ("LOSS" if med <= -args.threshold else
                     "GAIN" if med >= args.threshold else "NEUTRAL"),
            "used_as_baseline": "yes" if c in chrom_base else "no",
        })

    gene_rows, seg_rows = [], []
    uninformative = []
    for name, vals in sorted(by_region.items()):
        vals.sort(key=lambda v: (v[0], v[1]))
        chrom = vals[0][0]
        pos = np.array([v[1] for v in vals])
        l2 = np.array([v[2] for v in vals])          # within-sample
        l2r = np.array([v[3] for v in vals])         # reference-adjusted
        # A region is only callable where it was actually enriched. Where
        # most of a window sits at off-target depth the median is a statement
        # about the panel, not about copy number: a target whose BED interval
        # is wider than what the instrument enriched reads as a homozygous
        # deletion in every sample. Reported as uninformative rather than
        # dropped, so the gap is visible instead of absent.
        depths = np.array([smp["bins"][(chrom, int(pp))] for pp in pos])
        frac_cov = float(np.mean(depths >= auto_med * args.min_covered_fraction))
        if frac_cov < args.min_covered_bins:
            # A homozygous deletion is also empty, and CDKN2A/B is the case
            # that matters most. The two are separable only against a
            # reference: an unenriched region is empty in every sample, a
            # deleted one is empty in this sample alone. With no reference the
            # region is reported as uninformative and the ambiguity is stated,
            # rather than guessed either way.
            ref_cov = None
            if ref:
                rv = [ref.get((chrom, int(pp))) for pp in pos]
                rv = [x for x in rv if x is not None]
                if rv:
                    ref_cov = float(np.mean(
                        np.array(rv) >= args.min_covered_fraction))
            if ref_cov is not None and ref_cov >= args.min_covered_bins:
                # Covered in the reference, empty here: a real loss.
                med = -args.homdel_log2 - chrom_base.get(chrom, 0.0)
                gene_rows.append({
                    "sample": args.sample, "region": name, "chrom": chrom,
                    "start": int(pos.min()),
                    "end": int(pos.max() + args.binsize),
                    "n_bins": len(l2),
                    "median_log2": f"{med:.3f}",
                    "log2_vs_panel": f"{med:.3f}",
                    "chrom_baseline": f"{chrom_base.get(chrom, 0.0):.3f}",
                    "ratio": "0.000", "copy_number": "0.00",
                    "call": "HOMOZYGOUS_LOSS",
                    "n_reference_samples": str(len(args.reference_bam)),
                })
                continue
            uninformative.append({
                "sample": args.sample, "region": name, "chrom": chrom,
                "n_bins": len(l2),
                "frac_bins_covered": f"{frac_cov:.3f}",
                "median_depth": f"{float(np.median(depths)):.1f}",
                "autosomal_median_depth": f"{auto_med:.1f}",
                "reason": ("empty here and in the reference, so unenriched"
                           if ref else
                           "too little of the region reached on-target depth. "
                           "Without a reference an unenriched region and a "
                           "homozygous deletion are indistinguishable"),
            })
            continue

        raw_med = float(np.median(l2))
        base = chrom_base.get(chrom, 0.0)
        med = raw_med - base
        gene_rows.append({
            "sample": args.sample, "region": name, "chrom": chrom,
            "start": int(pos.min()), "end": int(pos.max() + args.binsize),
            "n_bins": len(l2),
            "median_log2": f"{med:.3f}",
            "log2_vs_panel": f"{raw_med:.3f}",
            "chrom_baseline": f"{base:.3f}",
            "ratio": f"{2 ** med:.3f}",
            "copy_number": f"{cn_from_log2(med, args.purity):.2f}",
            "call": ("LOSS" if med <= -args.threshold else
                     "GAIN" if med >= args.threshold else "NEUTRAL"),
            "n_reference_samples": str(len(args.reference_bam)),
        })
        # Segment threshold is the larger of the global floor and this
        # gene's own scatter, so a noisy region has to clear a higher bar.
        seg_thr = max(args.threshold, args.seg_sd * mad(l2r))
        for s, e, sl2, nb in segment(l2r, pos, seg_thr, args.min_bins):
            if nb == len(l2):
                continue        # whole-region shift, already in the gene row
            seg_rows.append({
                "sample": args.sample, "region": name, "chrom": chrom,
                "start": s, "end": e + args.binsize, "n_bins": nb,
                "span_bp": e + args.binsize - s,
                "median_log2": f"{sl2:.3f}",
                "copy_number": f"{cn_from_log2(sl2, args.purity):.2f}",
                "call": "LOSS" if sl2 < 0 else "GAIN",
                "note": "intragenic; confirm against the alignment",
            })

    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    gpath = Path(f"{args.out_prefix}.gene_cnv.tsv")
    spath = Path(f"{args.out_prefix}.gene_cnv_segments.tsv")
    for path, rows, cols in (
        (gpath, gene_rows, ["sample", "region", "chrom", "start", "end",
                            "n_bins", "median_log2", "log2_vs_panel",
                            "chrom_baseline", "ratio", "copy_number",
                            "call", "n_reference_samples"]),
        (spath, seg_rows, ["sample", "region", "chrom", "start", "end", "n_bins",
                           "span_bp", "median_log2", "copy_number", "call", "note"]),
    ):
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                               lineterminator="\n")
            w.writeheader(); w.writerows(rows)

    if uninformative:
        upath = Path(f"{args.out_prefix}.gene_cnv_uninformative.tsv")
        with open(upath, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(uninformative[0].keys()),
                               delimiter="\t", lineterminator="\n")
            w.writeheader(); w.writerows(uninformative)
        sys.stderr.write(
            f"  {len(uninformative)} region(s) uninformative -> {upath}\n")
        for u in uninformative:
            sys.stderr.write(f"    {u['region']:<22} "
                             f"{float(u['frac_bins_covered'])*100:5.1f}% of bins covered, "
                             f"median depth {u['median_depth']}x\n")

    cpath = Path(f"{args.out_prefix}.chrom_cnv.tsv")
    with open(cpath, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["sample", "chrom", "n_regions",
                                           "median_log2", "copy_number", "call",
                                           "used_as_baseline"],
                           delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(chrom_rows)
    for r in chrom_rows:
        if r["call"] != "NEUTRAL":
            sys.stderr.write(f"    {r['chrom']:<22} {r['call']:<8} "
                             f"log2 {r['median_log2']:>7}  CN {r['copy_number']} "
                             f"({r['n_regions']} regions)\n")

    n_called = sum(1 for r in gene_rows if r["call"] != "NEUTRAL")
    sys.stderr.write(
        f"{args.sample}: {len(gene_rows)} regions, {n_called} non-neutral, "
        f"{len(seg_rows)} intragenic segment(s)\n  {gpath}\n  {spath}\n")
    for r in gene_rows:
        if r["call"] != "NEUTRAL":
            sys.stderr.write(f"    {r['region']:<22} {r['call']:<8} "
                             f"log2 {r['median_log2']:>7}  CN {r['copy_number']}\n")
    for r in seg_rows:
        sys.stderr.write(f"    {r['region']:<22} {r['call']:<8} "
                         f"{r['span_bp']:>8,} bp  log2 {r['median_log2']}\n")

    # Plots.
    want = {g.strip().upper() for g in args.plot_genes.split(",") if g.strip()}
    if want:
        try:
            import matplotlib
            matplotlib.use("Agg")
            from matplotlib import pyplot as plt
        except ImportError:
            sys.stderr.write("matplotlib unavailable; skipping plots\n")
            return 0
        called = {r["region"].upper() for r in gene_rows if r["call"] != "NEUTRAL"}
        called |= {r["region"].upper() for r in seg_rows}
        pdir = Path(f"{args.out_prefix}.gene_cnv_plots")
        pdir.mkdir(parents=True, exist_ok=True)
        n = 0
        for name, vals in sorted(by_region.items()):
            u = name.upper()
            if not ("ALL" in want or (u in want)
                    or ("CALLED" in want and u in called)):
                continue
            vals.sort(key=lambda v: (v[0], v[1]))
            pos = np.array([v[1] for v in vals]) / 1e6
            l2 = np.array([v[2] for v in vals])
            plt.figure(figsize=(9, 3.2), dpi=140)
            plt.axhline(0, color="#999", lw=0.8)
            for y, c in ((-1.0, "#c00"), (0.585, "#06c")):
                plt.axhline(y, color=c, lw=0.6, ls=":")
            plt.scatter(pos, l2, s=6,
                        c=["#c00" if v <= -args.threshold else
                           "#06c" if v >= args.threshold else "#888" for v in l2])
            plt.ylim(-3.2, 2.2)
            plt.ylabel("log2 ratio")
            plt.xlabel(f"{vals[0][0]} (Mb)")
            plt.title(f"{args.sample}  {name}", fontsize=10)
            plt.tight_layout()
            plt.savefig(pdir / f"{args.sample}.{name.replace('/', '_')}.png")
            plt.close()
            n += 1
        sys.stderr.write(f"  {n} plot(s) -> {pdir}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
