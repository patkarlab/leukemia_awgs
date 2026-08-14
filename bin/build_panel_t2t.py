#!/usr/bin/env python3
"""
build_panel_t2t.py
==================

Build the T2T-CHM13v2.0 panel BEDs from a target table, an existing T2T BED,
and the NCBI RefSeq T2T GFF. Emits both namings from one build:

    <PANEL>_panel_t2t_chr.bed   chr-named, for analysis on gandalf
    <PANEL>_panel_t2t_NC.bed    NC_-named, for MinKNOW on the P2i

WHY THIS EXISTS
---------------
build_panel.py makes the hg38 BED from hg38 RefSeq. The T2T BEDs had no
generator, yet they are the ones the P2i enriches against, so they decide
whether reads exist at all. That asymmetry is how the T2T BEDs drifted out of
step with the design twice: sixteen genes stayed on the T2T panel after being
dropped from the target table, and PAR1_CRLF2_P2RY8_Y was added to the design
and never reached T2T.

This does NOT lift over from hg38. It reads gene spans directly from the T2T
RefSeq annotation, the same approach mm-awgs-nextflow uses in
bin/build_v6_panel.py. A chain has to thread coordinates through assembly
differences, and segmental duplications are exactly where that goes wrong.

CARRY FORWARD, DERIVE ONLY WHAT IS NEW
--------------------------------------
Regions already present in --prev-bed are copied byte-for-byte, so anything
the lab has already sequenced against keeps its exact coordinates. Only
targets absent from --prev-bed are derived from the GFF. This is deliberate:
a rebuild must not silently move an interval that a run was enriched against.

It also solves a problem a pure GFF build cannot. The named intervals in
build_panel.py (IGH_locus, TRB_locus, HOXA_cluster, PAR1_CRLF2_P2RY8 and the
rest) are not single RefSeq genes, so they have no GFF span to look up. They
already carry validated T2T coordinates in the existing BED, and carrying them
forward preserves those.

Targets in the table, absent from --prev-bed, and not resolvable from the GFF
are reported as unresolved. Named intervals newly added to the design fall
here: they need coordinates supplied via --extra-regions, because nothing can
derive them.

FORCED RE-DERIVATION
--------------------
--rederive SYMBOL[,SYMBOL...] drops a label from the carry-forward set and
rebuilds it from the GFF. Use when a carried region is wrong, e.g. a span
inflated by a segmental duplication.

SIZE GATE
---------
--hg38-bed enables a span comparison against the hg38 build. Assembly
differences are real, so exact agreement is not expected, but a region several
times its hg38 span is almost always an annotation artefact rather than a
coverage requirement. Reported for every region, carried or derived.
--max-ratio makes it fatal.

Usage
-----
  build_panel_t2t.py \\
      --panel     ALL \\
      --targets   assets/ALL_panel_targets.tsv \\
      --prev-bed  assets/ALL_panel_t2t_chr.bed \\
      --gff       /goast/hemat_data/references/T2T/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff \\
      --fai       /goast/nikhil_awgs_testing/t2t/refs/chm13v2.0.ucsc.fa.fai \\
      --hg38-bed  assets/ALL_panel_hg38.bed \\
      --outdir    /tmp/t2t_build

Writes nothing outside --outdir.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

__version__ = "0.1.0"

# NCBI RefSeq accessions for T2T-CHM13v2.0, autosomes and sex chromosomes.
# The GFF is NC_-named; analysis tooling wants chr-naming. Same coordinates.
NC_TO_CHR = {
    "NC_060925.1": "chr1",  "NC_060926.1": "chr2",  "NC_060927.1": "chr3",
    "NC_060928.1": "chr4",  "NC_060929.1": "chr5",  "NC_060930.1": "chr6",
    "NC_060931.1": "chr7",  "NC_060932.1": "chr8",  "NC_060933.1": "chr9",
    "NC_060934.1": "chr10", "NC_060935.1": "chr11", "NC_060936.1": "chr12",
    "NC_060937.1": "chr13", "NC_060938.1": "chr14", "NC_060939.1": "chr15",
    "NC_060940.1": "chr16", "NC_060941.1": "chr17", "NC_060942.1": "chr18",
    "NC_060943.1": "chr19", "NC_060944.1": "chr20", "NC_060945.1": "chr21",
    "NC_060946.1": "chr22", "NC_060947.1": "chrX",  "NC_060948.1": "chrY",
}
CHR_TO_NC = {v: k for k, v in NC_TO_CHR.items()}

# Label components that are not gene symbols. Compound labels arise where
# regions merged during a previous build ("TAL1/STIL").
SUFFIX = {"LOCUS", "CLUSTER", "ENHANCER", "REGION", "INTERVAL"}


@dataclass
class Region:
    chrom: str          # chr-named
    start: int          # 0-based
    end: int            # exclusive
    name: str
    origin: str         # "carried" | "derived" | "extra"


def tokens(label: str) -> set:
    up = label.upper().replace("+", "/").replace("_", "/")
    return {t for t in up.split("/") if t and t not in SUFFIX}


def open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def chrom_key(c: str) -> Tuple[int, object]:
    c2 = c.replace("chr", "")
    return (0, int(c2)) if c2.isdigit() else (1, c2)


# -----------------------------------------------------------------------------
# Inputs
# -----------------------------------------------------------------------------
def load_centered(builder: Optional[Path]):
    """Import CENTERED and its panel-aware resolver from build_panel.py.

    The coverage rules live in one place. Importing them rather than
    duplicating them is what stops the two references disagreeing: a gene given
    a centred window on hg38 must get the same window on T2T, or the panel says
    two different things depending on which BED you read.

    Returns a callable (name, panel) -> half-width or None. If build_panel.py
    cannot be found, returns a resolver that always says "not centred", so this
    degrades to the previous flat-flank behaviour rather than failing.
    """
    import importlib.util
    if not builder or not builder.exists():
        sys.stderr.write(
            f"WARNING: {builder} not found; no centred windows will be applied. "
            f"Any gene with a CENTERED entry will be built at the flat "
            f"--flank instead, and will disagree with the hg38 BED.\n")
        return lambda name, panel: None
    spec = importlib.util.spec_from_file_location("build_panel", builder)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    if hasattr(mod, "centered_halfwidth"):
        return mod.centered_halfwidth
    centered = getattr(mod, "CENTERED", {})
    return lambda name, panel: centered.get(name)


def load_targets(path: Path) -> List[Tuple[str, str]]:
    """(gene_or_region, group) for KEEP == Y rows, first nomination wins."""
    out, seen = [], set()
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if (row.get("KEEP") or "").strip().upper() != "Y":
                continue
            name = (row.get("GENE_OR_REGION") or "").strip()
            if name and name not in seen:
                seen.add(name)
                out.append((name, (row.get("GROUP") or "").strip()))
    return out


def load_bed(path: Path, origin: str) -> List[Region]:
    out = []
    if not path or not path.exists():
        return out
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith(("#", "track")):
            continue
        f = line.split("\t")
        if len(f) < 4:
            continue
        chrom = NC_TO_CHR.get(f[0], f[0])
        out.append(Region(chrom, int(f[1]), int(f[2]), f[3], origin))
    return out


def load_fai(path: Path) -> Dict[str, int]:
    sizes = {}
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        f = line.split("\t")
        sizes[NC_TO_CHR.get(f[0], f[0])] = int(f[1])
    return sizes


def parse_gff_attrs(field: str) -> Dict[str, str]:
    out = {}
    for part in field.rstrip(";").split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def parse_gff(path: Path, wanted: set) -> Dict[str, Region]:
    """Whole-gene spans from an NCBI RefSeq GFF3, keyed by uppercase symbol.

    A symbol on more than one contig prefers the NC_ primary placement, which
    keeps an alt or unplaced scaffold from winning.
    """
    wanted_u = {s.upper() for s in wanted}
    span: Dict[str, dict] = defaultdict(
        lambda: {"chrom": None, "start": None, "end": None, "symbol": None})
    with open_text(path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "gene":
                continue
            attrs = parse_gff_attrs(cols[8])
            sym = attrs.get("Name") or attrs.get("gene") or ""
            if sym.upper() not in wanted_u:
                continue
            chrom, start0, end1 = cols[0], int(cols[3]) - 1, int(cols[4])
            rec = span[sym.upper()]
            if rec["chrom"] is None:
                rec["chrom"], rec["symbol"] = chrom, sym
            elif rec["chrom"] != chrom:
                if chrom.startswith("NC_") and not rec["chrom"].startswith("NC_"):
                    rec["chrom"], rec["start"], rec["end"] = chrom, None, None
                else:
                    continue
            if rec["start"] is None or start0 < rec["start"]:
                rec["start"] = start0
            if rec["end"] is None or end1 > rec["end"]:
                rec["end"] = end1
    return {
        u: Region(NC_TO_CHR.get(r["chrom"], r["chrom"]),
                  r["start"], r["end"], r["symbol"], "derived")
        for u, r in span.items()
        if None not in (r["chrom"], r["start"], r["end"])
    }


def load_aliases(path: Optional[Path]) -> Dict[str, str]:
    out = {}
    if not path or not path.exists():
        return out
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            a = (row.get("alias") or "").strip().upper()
            s = (row.get("symbol") or "").strip().upper()
            if a and s:
                out[a] = s
                out[s] = a          # both directions: the GFF may use either
    return out


# -----------------------------------------------------------------------------
# Assembly
# -----------------------------------------------------------------------------
def clip(start: int, end: int, chrom: str, sizes: Dict[str, int]) -> Tuple[int, int]:
    return max(0, start), min(end, sizes.get(chrom, end))


def merge_overlaps(regions: List[Region]) -> List[Region]:
    """Merge overlapping intervals, joining names with '/'.

    Matches build_panel.py so the two references produce comparable labels.
    """
    by_chrom: Dict[str, List[Region]] = defaultdict(list)
    for r in regions:
        by_chrom[r.chrom].append(r)
    out: List[Region] = []
    for chrom in sorted(by_chrom, key=chrom_key):
        iv = sorted(by_chrom[chrom], key=lambda r: (r.start, r.end))
        cur = iv[0]
        names, origins = [cur.name], {cur.origin}
        cs, ce = cur.start, cur.end
        for r in iv[1:]:
            if r.start <= ce:
                ce = max(ce, r.end)
                if r.name not in names:
                    names.append(r.name)
                origins.add(r.origin)
            else:
                out.append(Region(chrom, cs, ce, "/".join(names),
                                  "mixed" if len(origins) > 1 else origins.pop()))
                cs, ce, names, origins = r.start, r.end, [r.name], {r.origin}
        out.append(Region(chrom, cs, ce, "/".join(names),
                          "mixed" if len(origins) > 1 else origins.pop()))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel", required=True, help="Panel name, used for output filenames.")
    ap.add_argument("--targets", required=True, type=Path)
    ap.add_argument("--prev-bed", required=True, type=Path,
                    help="Existing T2T BED. Regions here are carried forward "
                         "byte-for-byte, never re-derived.")
    ap.add_argument("--gff", required=True, type=Path,
                    help="NCBI RefSeq T2T-CHM13v2.0 GFF (NC_-named).")
    ap.add_argument("--fai", required=True, type=Path,
                    help="T2T FASTA index, for clipping to contig ends.")
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--flank", type=int, default=50000,
                    help="Flank either side of a derived gene body [50000]. "
                         "Matches FLANK in build_panel.py.")
    ap.add_argument("--aliases", type=Path, default=None,
                    help="assets/gene_aliases.tsv, for symbols the GFF names "
                         "differently from the target table.")
    ap.add_argument("--extra-regions", type=Path, default=None,
                    help="BED of coordinates for targets that are neither in "
                         "--prev-bed nor resolvable from the GFF, i.e. newly "
                         "added named intervals.")
    ap.add_argument("--builder", type=Path, default=None,
                    help="Path to build_panel.py, whose CENTERED table supplies "
                         "centred-window half-widths [alongside this script].")
    ap.add_argument("--rederive", default="",
                    help="Comma-separated labels to rebuild from the GFF even "
                         "if present in --prev-bed.")
    ap.add_argument("--hg38-bed", type=Path, default=None,
                    help="hg38 BED to compare spans against. Reporting only "
                         "unless --max-ratio is set.")
    ap.add_argument("--max-ratio", type=float, default=None,
                    help="Fail if any region exceeds its hg38 span by more "
                         "than this factor. Suggested 2.0.")
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args = ap.parse_args()

    targets = load_targets(args.targets)
    prev = load_bed(args.prev_bed, "carried")
    sizes = load_fai(args.fai)
    aliases = load_aliases(args.aliases)
    extra = load_bed(args.extra_regions, "extra") if args.extra_regions else []
    rederive = {s.strip().upper() for s in args.rederive.split(",") if s.strip()}

    # Which targets does the previous BED already cover?
    #
    # Matching is by EXACT label first, then by token. Exact-first matters
    # because labels can share every token and still be different regions:
    # PAR1_CRLF2_P2RY8 on chrX and PAR1_CRLF2_P2RY8_Y on chrY differ only by
    # the "Y". A token index maps PAR1, CRLF2 and P2RY8 to whichever of the two
    # is read second, so the other is claimed by nothing and silently dropped
    # from the rebuild. Token matching is still needed for compound labels a
    # previous merge produced ("TAL1/STIL"), so it is kept as a fallback and
    # only consulted when no exact label matches.
    ambiguous: List[str] = []
    prev_exact: Dict[str, Region] = {}
    prev_tokens: Dict[str, List[Region]] = defaultdict(list)
    for r in prev:
        prev_exact.setdefault(r.name.upper(), r)
        for tk in tokens(r.name):
            prev_tokens[tk].append(r)

    def find_prev(name: str) -> Optional[Region]:
        hit = prev_exact.get(name.upper())
        if hit is not None:
            return hit
        # Token fallback, for compound labels a previous merge produced.
        cands = []
        for tk in tokens(name):
            for r in prev_tokens.get(tk, []):
                if r not in cands:
                    cands.append(r)
        if len(cands) == 1:
            return cands[0]
        if len(cands) > 1:
            # Prefer a region whose label is exactly this target plus merge
            # partners, i.e. one of its "/"-separated components is the target
            # verbatim. That resolves PAR1_CRLF2_P2RY8 to
            # "PAR1_CRLF2_P2RY8/CRLF2" rather than to the chrY region, which
            # merely shares tokens.
            want = name.upper()
            component = [r for r in cands
                         if want in [c.strip().upper()
                                     for c in r.name.split("/")]]
            if len(component) == 1:
                return component[0]
            # Otherwise prefer the candidate whose full token set is the
            # closest superset; an exact token-set equality wins outright.
            exact_set = [r for r in cands if tokens(r.name) == tokens(name)]
            if len(exact_set) == 1:
                return exact_set[0]
            ambiguous.append(
                f"{name} -> {', '.join(sorted(c.name for c in cands))}")
        return None

    carried, to_derive, unresolved = [], [], []
    seen_prev: set = set()
    for name, _group in targets:
        if tokens(name) & rederive:
            to_derive.append(name)
            continue
        hit = find_prev(name)
        if hit is not None:
            if id(hit) not in seen_prev:
                seen_prev.add(id(hit))
                carried.append(hit)
        else:
            to_derive.append(name)

    # Anything in the previous BED whose target row is gone has been dropped
    # from the design and must not be carried.
    target_tokens = set()
    for name, _ in targets:
        target_tokens |= tokens(name)
    dropped = [r for r in prev
               if not (tokens(r.name) & target_tokens) and id(r) not in seen_prev]

    # Derive the rest from the GFF.
    wanted = set()
    for n in to_derive:
        wanted |= tokens(n)
        for t in tokens(n):
            if t in aliases:
                wanted.add(aliases[t])
    coords = parse_gff(args.gff, wanted) if wanted else {}
    centered = load_centered(args.builder or (Path(__file__).parent / "build_panel.py"))

    extra_exact = {r.name.upper(): r for r in extra}
    extra_tokens: Dict[str, List[Region]] = defaultdict(list)
    for r in extra:
        for tk in tokens(r.name):
            extra_tokens[tk].append(r)

    def find_extra(name: str) -> Optional[Region]:
        hit = extra_exact.get(name.upper())
        if hit is not None:
            return hit
        cands = []
        for tk in tokens(name):
            for r in extra_tokens.get(tk, []):
                if r not in cands:
                    cands.append(r)
        return cands[0] if len(cands) == 1 else None

    derived: List[Region] = []
    centred_note: List[str] = []
    for name in to_derive:
        tk = tokens(name)
        hit = find_extra(name)
        if hit:
            derived.append(Region(hit.chrom, hit.start, hit.end, name, "extra"))
            continue
        gc = None
        for t in tk:
            gc = coords.get(t) or coords.get(aliases.get(t, ""))
            if gc:
                break
        if not gc:
            unresolved.append(name)
            continue
        hw = centered(name, args.panel)
        if hw is not None:
            # Centred on the gene-body midpoint, matching build_panel.py, so
            # the same locus gets the same width on both references.
            mid = (gc.start + gc.end) // 2
            s, e = clip(mid - hw, mid + hw, gc.chrom, sizes)
            centred_note.append(f"{name} (+/- {hw:,})")
        else:
            s, e = clip(gc.start - args.flank, gc.end + args.flank,
                        gc.chrom, sizes)
        derived.append(Region(gc.chrom, s, e, name, "derived"))

    final = merge_overlaps(carried + derived)
    final.sort(key=lambda r: (chrom_key(r.chrom), r.start))

    args.outdir.mkdir(parents=True, exist_ok=True)
    chr_path = args.outdir / f"{args.panel}_panel_t2t_chr.bed"
    nc_path = args.outdir / f"{args.panel}_panel_t2t_NC.bed"
    chr_path.write_text("".join(
        f"{r.chrom}\t{r.start}\t{r.end}\t{r.name}\n" for r in final))
    nc_path.write_text("".join(
        f"{CHR_TO_NC.get(r.chrom, r.chrom)}\t{r.start}\t{r.end}\t{r.name}\n"
        for r in final))

    total = sum(r.end - r.start for r in final)
    genome = sum(v for k, v in sizes.items() if k in CHR_TO_NC)
    w = sys.stderr.write
    w(f"panel            : {args.panel}\n")
    w(f"targets in list  : {len(targets)}\n")
    w(f"carried forward  : {len(carried)}\n")
    w(f"derived from GFF : {len([r for r in derived if r.origin == 'derived'])}\n")
    w(f"from extra BED   : {len([r for r in derived if r.origin == 'extra'])}\n")
    if centred_note:
        w(f"centred windows  : {', '.join(centred_note)}\n")
    w(f"regions emitted  : {len(final)}\n")
    w(f"total bases      : {total:,} ({100 * total / genome:.3f}% of T2T primary)\n")
    w(f"written          : {chr_path}\n                   {nc_path}\n")

    if ambiguous:
        w(f"\nAMBIGUOUS ({len(ambiguous)}): target matched more than one region "
          f"by token and none exactly, so it was not carried forward.\n")
        for a in ambiguous:
            w(f"  {a}\n")
        w("  Give the target a label matching its region exactly, or supply "
          "the interval via --extra-regions.\n")
    if dropped:
        w(f"\ndropped from the previous BED, no longer in the target table "
          f"({len(dropped)}):\n  " + ", ".join(sorted(r.name for r in dropped)) + "\n")
    if unresolved:
        w(f"\nUNRESOLVED ({len(unresolved)}): not in --prev-bed, not a gene in "
          f"the GFF, not in --extra-regions.\n  " + ", ".join(unresolved) + "\n"
          f"  Named intervals must be supplied via --extra-regions; nothing "
          f"can derive them.\n")

    # Span comparison against hg38.
    fatal = False
    if args.hg38_bed and args.hg38_bed.exists():
        hg = {}
        for r in load_bed(args.hg38_bed, "hg38"):
            for t in tokens(r.name):
                hg[t] = r.end - r.start
        rows = []
        for r in final:
            ref = max((hg[t] for t in tokens(r.name) if t in hg), default=None)
            if not ref:
                continue
            span = r.end - r.start
            ratio = span / ref
            if ratio >= 1.5 or ratio <= 0.67:
                rows.append((ratio, r.name, span, ref, r.origin))
        if rows:
            w(f"\nspan differs from hg38 by more than 1.5x ({len(rows)}):\n")
            for ratio, name, span, ref, origin in sorted(rows, reverse=True):
                flag = ""
                if args.max_ratio and ratio > args.max_ratio:
                    flag, fatal = "  EXCEEDS --max-ratio", True
                w(f"  {name:<30} T2T {span:>9,}  hg38 {ref:>9,}  "
                  f"x{ratio:.2f}  [{origin}]{flag}\n")
            w("  Assembly differences are real; a large ratio on a carried "
              "region usually means an annotation artefact worth re-deriving "
              "with --rederive.\n")

    if fatal:
        w("\nFAILED: at least one region exceeds --max-ratio. Nothing was "
          "deleted; inspect the output before using it.\n")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
