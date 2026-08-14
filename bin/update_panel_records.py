#!/usr/bin/env python3
"""
update_panel_records.py
=======================

Regenerates the two derived records that drift every time a panel changes:

  <PANEL>_panel_targets.tsv   the BASES column
  panel_summary.tsv           per-group target and base counts, both panels

WHY
---
BASES and panel_summary.tsv are computed, not authored, but nothing computed
them. They were maintained by hand and went stale on every panel revision:
after RANBP2 was removed and CDKN2A/CDKN2B added, panel_summary still reported
the counts from two revisions earlier, and the two new rows carried a
placeholder of 0. Neither is read by the pipeline, which is exactly why the
drift went unnoticed -- and also why it matters, since these are the files
someone reads later to reconstruct why a panel looks the way it does.

This imports build_panel.py rather than reimplementing it, so FLANK, CENTERED
and NAMED_REGIONS cannot diverge between the builder and the records. If the
coverage rules change, both follow automatically.

BASES is the per-target span BEFORE merging, matching what build_panel.py
accumulates into its per-group counters. It is therefore not the sum of the
emitted BED regions: where two targets overlap they merge into one region,
and the merged span is smaller than the two individual spans added together.

Usage:
  update_panel_records.py \\
      --assets  assets \\
      --refgene /path/refGene.txt.gz \\
      --sizes   /path/hg38.chrom.sizes

  update_panel_records.py --assets assets ... --dry-run
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import sys
from collections import OrderedDict, defaultdict
from pathlib import Path

__version__ = "0.1.0"


def load_builder(path: Path):
    """Import build_panel.py as a module so its coverage rules are reused.

    A plain `import build_panel` would depend on sys.path and on the file
    being importable by name; loading it by location keeps this working from
    any working directory.
    """
    spec = importlib.util.spec_from_file_location("build_panel", path)
    if spec is None or spec.loader is None:
        sys.stderr.write(f"ERROR: cannot load {path}\n")
        sys.exit(1)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    for attr in ("FLANK", "CENTERED", "NAMED_REGIONS", "load_refgene",
                 "load_sizes", "clip"):
        if not hasattr(mod, attr):
            sys.stderr.write(
                f"ERROR: {path} has no '{attr}'. This script expects the "
                f"coverage rules to live in build_panel.py.\n")
            sys.exit(1)
    return mod


def target_span(name: str, bp, rg, sizes):
    """Pre-merge span for one target, by exactly build_panel.py's rules.

    Returns (bases, note). note is empty on success, otherwise why it is 0.
    """
    if name in bp.NAMED_REGIONS:
        chrom, s, e = bp.NAMED_REGIONS[name]
        s, e = bp.clip(s, e, chrom, sizes)
        return e - s, ""
    if name not in rg:
        return 0, "not in refGene"
    chrom, tx_s, tx_e = rg[name]
    if name in bp.CENTERED:
        hw = bp.CENTERED[name]
        mid = (tx_s + tx_e) // 2
        s, e = bp.clip(mid - hw, mid + hw, chrom, sizes)
    else:
        s, e = bp.clip(tx_s - bp.FLANK, tx_e + bp.FLANK, chrom, sizes)
    if e <= s:
        # clip() bounds the end to the contig length but not the start, so a
        # coordinate past the contig end yields an inverted interval. Real
        # data should never hit this; a mismatched refGene and chrom.sizes
        # pair will.
        return 0, f"interval collapsed ({chrom}:{s}-{e}); refGene and "
        f"chrom.sizes may be from different assemblies"
    return e - s, ""


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--assets", type=Path, default=Path("assets"))
    ap.add_argument("--builder", type=Path, default=None,
                    help="Path to build_panel.py [<assets>/../bin/build_panel.py]")
    ap.add_argument("--refgene", required=True)
    ap.add_argument("--sizes", required=True)
    ap.add_argument("--panels", default="AML,ALL")
    ap.add_argument("--dry-run", action="store_true",
                    help="Report what would change, write nothing.")
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args = ap.parse_args()

    builder = args.builder or (args.assets.parent / "bin" / "build_panel.py")
    bp = load_builder(builder)
    rg = bp.load_refgene(args.refgene)
    sizes = bp.load_sizes(args.sizes)
    panels = [p.strip() for p in args.panels.split(",") if p.strip()]

    summary_rows = []
    for panel in panels:
        path = args.assets / f"{panel}_panel_targets.tsv"
        if not path.exists():
            sys.stderr.write(f"skipping {panel}: {path} not found\n")
            continue

        rows = list(csv.DictReader(open(path), delimiter="\t"))
        cols = list(rows[0].keys())
        changed, missing = [], []
        groups = OrderedDict()

        for r in rows:
            if (r.get("KEEP") or "").strip().upper() != "Y":
                continue
            name = r["GENE_OR_REGION"].strip()
            bases, note = target_span(name, bp, rg, sizes)
            if note:
                missing.append(f"{name} ({note})")
            old = (r.get("BASES") or "").strip()
            if str(bases) != old:
                changed.append((name, old or "-", bases))
                r["BASES"] = str(bases)
            g = r["GROUP"].strip()
            if g not in groups:
                groups[g] = [0, 0]
            groups[g][0] += 1
            groups[g][1] += bases

        print(f"\n{panel}: {sum(v[0] for v in groups.values())} targets, "
              f"{sum(v[1] for v in groups.values()):,} bases (pre-merge)")
        if changed:
            print(f"  BASES corrected on {len(changed)} row(s):")
            for name, old, new in changed:
                print(f"    {name:<24} {old:>10} -> {new:>10}")
        else:
            print("  BASES already correct on every row")
        if missing:
            print(f"  UNRESOLVED, BASES left at 0: {', '.join(missing)}")
            print("  These are dropped by build_panel.py too. Check for a "
                  "symbol rename and add it to assets/gene_aliases.tsv.")

        for g, (n, b) in groups.items():
            summary_rows.append({"PANEL": panel, "GROUP": g,
                                 "TARGETS": str(n), "BASES": str(b)})

        if not args.dry_run:
            with open(path, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                                   lineterminator="\n")
                w.writeheader(); w.writerows(rows)

    if summary_rows and not args.dry_run:
        out = args.assets / "panel_summary.tsv"
        with open(out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["PANEL", "GROUP", "TARGETS", "BASES"],
                               delimiter="\t", lineterminator="\n")
            w.writeheader(); w.writerows(summary_rows)
        print(f"\nwrote {out} ({len(summary_rows)} group rows)")
    elif args.dry_run:
        print("\ndry run: nothing written")
    return 0


if __name__ == "__main__":
    sys.exit(main())
