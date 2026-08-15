#!/usr/bin/env python3
"""
make_focal_dup_sites.py
=======================

Reshape focal-duplication caller output into the chrom/pos sites table that
igv_snapshots.py --mode somatic consumes, so FLT3-ITD, UBTF-TD and KMT2A-PTD
each get an alignment page.

A standalone script rather than a heredoc inside the Nextflow module: escapes
in an embedded script are interpolated by Groovy before the shell sees them, so
a \\t written for Python arrives as a literal tab and the source no longer
parses. Keeping it in a file means the escapes are the file's own.

Usage:
  make_focal_dup_sites.py OUT.tsv IN1.tsv [IN2.tsv ...]

Inputs may be either duplication table: the tandem caller names its position
column ref_pos and the intragenic one acceptor_pos, so the position is resolved
by column name rather than by index. Paths ending in NO_FILE are the pipeline's
placeholder for "this sample had none" and are skipped.
"""

import csv
import sys

COLS = ["chrom", "pos", "gene", "label", "dup_length_bp", "confidence", "hotspot"]


def main() -> int:
    if len(sys.argv) < 2:
        sys.stderr.write(__doc__)
        return 2
    out = sys.argv[1]
    srcs = [s for s in sys.argv[2:] if not s.endswith("NO_FILE")]

    sites = []
    for src in srcs:
        try:
            with open(src, newline="") as fh:
                rows = list(csv.DictReader(fh, delimiter="\t"))
        except OSError as exc:
            sys.stderr.write(f"  skipping {src}: {exc}\n")
            continue
        for r in rows:
            pos = r.get("ref_pos") or r.get("acceptor_pos")
            if not r.get("chrom") or not pos:
                continue
            sites.append({
                "chrom": r["chrom"],
                "pos": pos,
                "gene": r.get("gene", ""),
                "label": r.get("label", ""),
                "dup_length_bp": (r.get("dup_length_bp")
                                  or r.get("duplication_span_bp", "")),
                "confidence": r.get("confidence", ""),
                "hotspot": r.get("hotspot", ""),
            })

    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(sites)

    sys.stderr.write(
        f"focal duplication sites: {len(sites)} from {len(srcs)} table(s)\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
