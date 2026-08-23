#!/usr/bin/env python3
"""
build_gene_model.py

Emit one BED feature per gene, bare gene bodies, no flank and no merging.

WHY THIS EXISTS
---------------
The panel BED serves capture: intervals are flanked and overlapping ones are
merged, which is correct for enrichment and wrong for naming a breakpoint.
Merging produces compound labels ("PAX5/ZCCHC7") and, worse, puts both ends of
an intragenic fusion in one feature so gene_a == gene_b and no dictionary pair
can match. PAX5::ZCCHC7 and STIL::TAL1 are unmatchable for this reason.

This model is for annotation only. It never gates panel membership; that stays
with the panel BED. It only answers "which gene is this coordinate in".

Genes absent from the GFF (P2RY8 has no chrX record in T2T RefSeq) are taken
from --extra-regions, whose fourth column must be a bare gene symbol.

Usage:
  build_gene_model.py \\
      --targets  assets/ALL_panel_targets.tsv \\
      --gff      <T2T RefSeq GFF> \\
      --rename   assets/t2t_rename_map.tsv \\
      --extra-regions assets/extra_regions_ALL_t2t.bed \\
      --output   assets/ALL_gene_model_t2t.bed
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

__version__ = "0.1.0"

NAME_RE = re.compile(r"(?:^|;)Name=([^;]+)")


def load_targets(path: Path) -> set:
    """Gene symbols from the target table. Named intervals are skipped: they
    are not single genes and have no GFF span."""
    wanted = set()
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            keep, group, symbol = parts[0].strip(), parts[3].strip(), parts[4].strip()
            if keep.upper() == "KEEP":          # header row
                continue
            if keep.upper() != "Y":
                continue
            if group == "named_intervals" or not symbol:
                continue
            wanted.add(symbol.upper())
    return wanted


def load_rename(path):
    """NC_ accession -> chr name. Optional; identity if absent."""
    mapping = {}
    if path is None:
        return mapping
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[0] and parts[1]:
                mapping[parts[0]] = parts[1]
    return mapping


def parse_gff(path: Path, wanted: set, rename: dict):
    """Return {symbol: (chrom, start, end)} for wanted genes.

    A symbol on more than one contig keeps the longer record and reports the
    collision; PAR genes are annotated on both chrX and chrY.
    """
    found = {}
    collisions = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            m = NAME_RE.search(f[8])
            if not m:
                continue
            sym = m.group(1).strip()
            if sym.upper() not in wanted:
                continue
            chrom = rename.get(f[0], f[0])
            start, end = int(f[3]) - 1, int(f[4])
            prev = found.get(sym)
            if prev is None:
                found[sym] = (chrom, start, end)
            else:
                collisions.append(f"{sym}: {prev[0]} and {chrom}")
                if (end - start) > (prev[2] - prev[1]):
                    found[sym] = (chrom, start, end)
    return found, collisions


def load_extra(path, wanted):
    """Bare gene symbols supplied by coordinate. Fourth column must be a single
    symbol; compound labels are rejected, since a compound entry here would
    reintroduce exactly the defect this file exists to remove."""
    out = {}
    if path is None:
        return out
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith(("#", "track")):
                continue
            f = line.split("\t")
            if len(f) < 4:
                continue
            name = f[3].strip()
            if "/" in name or name.upper() not in wanted:
                continue
            out[name] = (f[0], int(f[1]), int(f[2]))
    return out


def chrom_key(chrom):
    name = chrom[3:] if chrom.lower().startswith("chr") else chrom
    if name.isdigit():
        return (0, int(name), "")
    return (1, {"X": 0, "Y": 1, "M": 2, "MT": 2}.get(name.upper(), 99), name.upper())


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--targets", required=True, type=Path)
    ap.add_argument("--gff", required=True, type=Path)
    ap.add_argument("--rename", type=Path, default=None,
                    help="TSV mapping GFF contig accessions to chr names.")
    ap.add_argument("--extra-regions", type=Path, default=None,
                    help="BED of genes the GFF cannot place. Bare symbols only.")
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--version", action="version",
                    version=f"%(prog)s {__version__}")
    args = ap.parse_args()

    wanted = load_targets(args.targets)
    rename = load_rename(args.rename)
    genes, collisions = parse_gff(args.gff, wanted, rename)
    extra = load_extra(args.extra_regions, wanted)

    # Extras win: they exist because the GFF placed the gene wrongly or not
    # at all.
    overridden = sorted(set(extra) & set(genes))
    genes.update(extra)

    missing = sorted(w for w in wanted
                     if w not in {g.upper() for g in genes})

    rows = sorted(((c, s, e, sym) for sym, (c, s, e) in genes.items()),
                  key=lambda r: (chrom_key(r[0]), r[1], r[2]))
    with open(args.output, "w") as out:
        for c, s, e, sym in rows:
            out.write(f"{c}\t{s}\t{e}\t{sym}\n")

    w = sys.stderr.write
    w(f"targets (genes)  : {len(wanted)}\n")
    w(f"from GFF         : {len(genes) - len(extra)}\n")
    w(f"from extras      : {len(extra)}\n")
    w(f"written          : {len(rows)} -> {args.output}\n")
    if overridden:
        w(f"extras overrode GFF ({len(overridden)}): {', '.join(overridden)}\n")
    if collisions:
        w(f"\nmulti-contig symbols ({len(collisions)}), longer record kept:\n  "
          + "\n  ".join(collisions) + "\n")
    if missing:
        w(f"\nNOT RESOLVED ({len(missing)}): absent from the GFF and from "
          f"--extra-regions.\n  " + ", ".join(missing) + "\n"
          "  A breakpoint in one of these falls back to the panel label.\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
