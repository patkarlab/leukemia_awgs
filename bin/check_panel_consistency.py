#!/usr/bin/env python3
"""
check_panel_consistency.py
==========================

Cross-checks a panel's four descriptions against each other and reports the
mismatches that would otherwise fail silently.

WHY
---
A panel is described in four places, each consumed by a different part of the
pipeline:

  <PANEL>_panel_targets.tsv    the design table; what was intended
  <PANEL>_panel_t2t_chr.bed    fusion calling, on-target QC, intragenic
                               duplication
  <PANEL>_panel_t2t_NC.bed     the MinKNOW adaptive-sampling configuration
  <PANEL>_panel_hg38.bed       ClairS-TO calling region, ichorCNA mask, the
                               SNV panel filter, tandem duplication

Nothing in the pipeline compares them. A gene present in the T2T BED but not
the hg38 BED still gets its fusions called and never appears in the SNV
report, and no log line says so. That asymmetry is what makes this worth an
explicit check rather than a code review.

Entirely data-driven: it reads whatever files it is given and holds no gene
list of its own.

Usage:
  check_panel_consistency.py --panel AML
  check_panel_consistency.py --panel ALL --dictionary --strict
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

# Label components that are not gene symbols. Compound labels arise where
# regions merged during panel construction ('TAL1/STIL'), so labels are
# tokenised rather than compared as strings.
SUFFIX = {"LOCUS", "CLUSTER", "ENHANCER", "REGION", "INTERVAL"}


def tokens(label: str) -> set:
    up = label.upper().replace("+", "/").replace("_", "/")
    return {t for t in up.split("/") if t and t not in SUFFIX}


def read_bed(path: Path):
    """Return (gene tokens, region count, total bases, bases per chromosome).

    Per-chromosome totals matter because token comparison alone cannot see a
    partial coverage difference: a region can carry the right gene label on
    both references while spanning a much smaller interval on one of them.
    """
    tok, n, bases = set(), 0, 0
    per_chrom: dict = {}
    if not path.exists():
        return tok, n, bases, per_chrom
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith(("#", "track")):
            continue
        f = line.split("\t")
        if len(f) < 4:
            continue
        n += 1
        width = int(f[2]) - int(f[1])
        bases += width
        per_chrom[f[0]] = per_chrom.get(f[0], 0) + width
        tok |= tokens(f[3])
    return tok, n, bases, per_chrom


def compare_chrom_coverage(a: dict, b: dict, label_a: str, label_b: str,
                           tolerance: float):
    """Chromosomes whose covered bases differ by more than `tolerance`.

    Coordinates differ between references, so exact agreement is not expected;
    a large relative difference means the two BEDs are describing different
    designs rather than the same design in different coordinates.
    """
    out = []
    for chrom in sorted(set(a) | set(b), key=lambda c: (len(c), c)):
        va, vb = a.get(chrom, 0), b.get(chrom, 0)
        ref = max(va, vb)
        if ref == 0:
            continue
        if abs(va - vb) / ref > tolerance:
            out.append(f"{chrom}: {label_a} {va:,} bp vs {label_b} {vb:,} bp")
    return out


def read_targets(path: Path):
    tok, n = set(), 0
    if not path.exists():
        return tok, n
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if (row.get("KEEP") or "").strip().upper() != "Y":
                continue
            n += 1
            tok |= tokens(row.get("GENE_OR_REGION") or "")
    return tok, n


def read_dictionary(path: Path, panel: str):
    """Gene tokens the dictionary names for this panel.

    Rows whose off-panel partner carries an expected cytoband stay resolvable
    without that partner on the panel, so only the on-panel side is required.
    """
    tok = set()
    if not path.exists():
        return tok
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            disease = (row.get("disease") or "BOTH").strip().upper()
            if disease not in ("BOTH", panel.upper()):
                continue
            tok |= tokens(row.get("partner_a") or "")
            if not (row.get("partner_b_band") or "").strip():
                tok |= tokens(row.get("partner_b") or "")
    return tok


def report(title: str, items) -> None:
    items = sorted(items)
    if items:
        print(f"\n{title} ({len(items)}):")
        print("  " + ", ".join(items))


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel", required=True,
                    help="Panel name, matching the assets/<PANEL>_panel_* prefix.")
    ap.add_argument("--assets", type=Path, default=Path("assets"))
    ap.add_argument("--dictionary", action="store_true",
                    help="Also check that every gene the fusion dictionary "
                         "names for this panel is actually on it.")
    ap.add_argument("--coverage-tolerance", type=float, default=0.30,
                    help="Relative per-chromosome base difference between the "
                         "T2T and hg38 BEDs above which a chromosome is "
                         "flagged [0.30]. Coordinates differ between "
                         "references, so exact agreement is not expected.")
    ap.add_argument("--strict", action="store_true",
                    help="Exit non-zero on any mismatch. Without it, only a "
                         "missing or empty BED is an error.")
    args = ap.parse_args()

    a, p = args.assets, args.panel
    paths = {
        "targets": a / f"{p}_panel_targets.tsv",
        "t2t_chr": a / f"{p}_panel_t2t_chr.bed",
        "t2t_NC":  a / f"{p}_panel_t2t_NC.bed",
        "hg38":    a / f"{p}_panel_hg38.bed",
    }

    tgt_tok, tgt_n = read_targets(paths["targets"])
    chr_tok, chr_n, chr_b, chr_pc = read_bed(paths["t2t_chr"])
    nc_tok,  nc_n,  nc_b, _nc_pc = read_bed(paths["t2t_NC"])
    h_tok,   h_n,   h_b, h_pc = read_bed(paths["hg38"])

    print(f"{p} panel")
    print(f"  design table {paths['targets'].name:<30} {tgt_n:>4} targets")
    print(f"  T2T chr      {paths['t2t_chr'].name:<30} {chr_n:>4} regions  "
          f"{chr_b / 1e6:6.2f} Mb")
    print(f"  T2T NC_      {paths['t2t_NC'].name:<30} {nc_n:>4} regions  "
          f"{nc_b / 1e6:6.2f} Mb")
    print(f"  hg38         {paths['hg38'].name:<30} {h_n:>4} regions  "
          f"{h_b / 1e6:6.2f} Mb")

    errors, warnings = [], []
    counts = {"t2t_chr": chr_n, "t2t_NC": nc_n, "hg38": h_n}

    for key in ("t2t_chr", "t2t_NC", "hg38"):
        if not paths[key].exists():
            errors.append(f"{paths[key].name} does not exist")
        elif counts[key] == 0:
            errors.append(
                f"{paths[key].name} has no intervals. If this is the shipped "
                f"placeholder, replace it: ClairS-TO would call nothing and "
                f"the SNV report would be empty.")

    if chr_n and nc_n:
        if chr_tok != nc_tok:
            warnings.append("T2T chr- and NC_-named BEDs disagree on genes")
            report("In chr but not NC_", chr_tok - nc_tok)
            report("In NC_ but not chr", nc_tok - chr_tok)
        if chr_b != nc_b:
            warnings.append(
                f"T2T chr and NC_ BEDs cover different base counts "
                f"({chr_b:,} vs {nc_b:,}); they should be the same intervals "
                f"under two naming conventions")

    if chr_n and h_n:
        only_t2t = chr_tok - h_tok
        only_hg38 = h_tok - chr_tok
        if only_t2t:
            warnings.append("genes on T2T but not hg38")
            report("In T2T but not hg38 - fusions called, no SNV reporting",
                   only_t2t)
        if only_hg38:
            warnings.append("genes on hg38 but not T2T")
            report("In hg38 but not T2T - SNVs reported, no fusion detection",
                   only_hg38)

        drift = compare_chrom_coverage(h_pc, chr_pc, "hg38", "T2T",
                                       args.coverage_tolerance)
        if drift:
            warnings.append("per-chromosome coverage differs between "
                            "references beyond tolerance")
            print(f"\nCoverage differs by more than "
                  f"{args.coverage_tolerance:.0%} on ({len(drift)}):")
            for d in drift:
                print(f"  {d}")

    if tgt_n:
        if chr_n and (missing := tgt_tok - chr_tok):
            warnings.append("design targets absent from the T2T BED")
            report("In the design table but not the T2T BED", missing)
        if h_n and (missing := tgt_tok - h_tok):
            warnings.append("design targets absent from the hg38 BED")
            report("In the design table but not the hg38 BED", missing)

    if args.dictionary:
        dict_tok = read_dictionary(a / "al_fusion_dictionary.tsv", p)
        if dict_tok and chr_n and (orphan := dict_tok - chr_tok):
            warnings.append("dictionary genes absent from the panel")
            report("Named by the fusion dictionary but not on this panel "
                   "(these pairs can never be called here)", orphan)

    print()
    for e in errors:
        print(f"ERROR: {e}")
    for w in warnings:
        print(f"WARNING: {w}")
    if not errors and not warnings:
        print("All descriptions of this panel agree.")

    if errors:
        return 1
    return 1 if (warnings and args.strict) else 0


if __name__ == "__main__":
    sys.exit(main())
