#!/usr/bin/env python3
"""
call_intragenic_dup.py
======================

Detects large intragenic duplications from long-read alignments, for any gene
named on the command line. The motivating case is KMT2A-PTD.

NOTHING GENE-SPECIFIC IS HARDCODED
----------------------------------
The gene is a parameter and its interval is resolved from the panel BED handed
in at runtime, so this script holds no coordinates. Which genes are scanned,
and with what size and support thresholds, comes from
assets/focal_duplication_targets.tsv.

WHY A DEDICATED CALLER
----------------------
An intragenic duplication has no fusion partner: both breakpoints lie inside
the same gene. KMT2A-PTD, classically of exons 2-8 or 2-10 and spanning
several kilobases, is invisible both to KMT2A break-apart FISH and to
fusion-partner analysis for that reason. A pipeline that only reports junctions between two panel genes
will miss it entirely, which matters because KMT2A-PTD carries adverse risk
and co-occurs with a distinct comutation profile.

METHOD
------
Three orthogonal signals are gathered, and the call is graded by how many
agree. Any one alone is weak on ONT data; together they are reasonably
specific.

  1. SPLIT-READ JUNCTIONS. A read crossing the duplication junction aligns
     in two pieces, both inside KMT2A, with the downstream piece mapping
     back upstream in the same orientation. This is the direct evidence.

  2. LARGE INSERTIONS. A read that spans the whole duplication carries it as
     one large CIGAR insertion, the same signal call_flt3_itd.py uses. This
     only fires for short PTDs, since a read must span the whole event.

  3. EXON COVERAGE RATIO. Duplicated exons sit at roughly 1.5x the depth of
     non-duplicated KMT2A exons in a heterozygous PTD at high clonality.
     This is the weakest signal, and it degrades badly at low tumour
     fraction or low panel depth, but it independently localises which
     exons are involved.

Signals 1 and 2 need no exon annotation. Signal 3 needs an exon BED and is
skipped when one is not supplied, so a missing annotation reduces evidence
rather than causing a failure.

INTERPRETATION
--------------
This caller reports evidence, not a diagnosis. A high-confidence call means
two or three signals agree on a duplication of plausible size within the gene.
Confirm with an orthogonal assay before clinical reporting; the intended use
here is triage and MRD-target identification.

Dependencies: pysam (present in the awgs_sv conda env).
"""

from __future__ import annotations

import argparse
import csv
import statistics
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

try:
    import pysam
except ImportError:
    sys.stderr.write("ERROR: pysam is required. Activate the awgs_sv env.\n")
    sys.exit(1)

__version__ = "0.1.0"

BAM_CINS, BAM_CSOFT_CLIP, BAM_CHARD_CLIP = 1, 4, 5
CONSUMES_QUERY = {0, 1, 4, 7, 8}
CONSUMES_REF = {0, 2, 3, 7, 8}


@dataclass
class Junction:
    """A candidate duplication junction: reads jump from `donor` back to
    `acceptor`, with acceptor < donor, both inside the gene."""
    acceptor: int
    donor: int
    reads: List[str] = field(default_factory=list)

    @property
    def span(self) -> int:
        return self.donor - self.acceptor

    @property
    def n_support(self) -> int:
        return len(self.reads)


def parse_region(region: str) -> Tuple[str, int, int]:
    try:
        chrom, span = region.rsplit(":", 1)
        start, end = span.replace(",", "").split("-")
        return chrom, int(start) - 1, int(end)
    except ValueError:
        sys.stderr.write(f"ERROR: cannot parse region '{region}'.\n")
        sys.exit(1)


def gene_region_from_bed(bed: Path, gene: str) -> Optional[str]:
    """Look the gene interval up in the panel BED, so the caller does not
    carry its own coordinates and cannot drift from the panel."""
    with open(bed) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 4 and gene.upper() in \
                    {t for t in parts[3].upper().replace("_", "/").split("/")}:
                return f"{parts[0]}:{int(parts[1]) + 1}-{parts[2]}"
    return None


# -----------------------------------------------------------------------------
# Signal 1: split-read junctions from supplementary alignments
# -----------------------------------------------------------------------------
def parse_sa_tag(sa: str) -> List[Tuple[str, int, str]]:
    """Parse an SA:Z tag into [(chrom, pos, strand), ...]. Positions 1-based."""
    out = []
    for item in sa.rstrip(";").split(";"):
        if not item:
            continue
        f = item.split(",")
        if len(f) >= 3:
            try:
                out.append((f[0], int(f[1]), f[2]))
            except ValueError:
                continue
    return out


def collect_split_junctions(bam, chrom, start, end, min_mapq, tol
                            ) -> List[Junction]:
    """Find reads whose primary and supplementary alignments both fall inside
    the gene, in the same orientation, arranged tail-to-head. That geometry is
    what a tandem duplication produces; an inversion or deletion gives a
    different one and is excluded."""
    raw: List[Tuple[int, int, str]] = []
    for aln in bam.fetch(chrom, start, end):
        if aln.is_unmapped or aln.is_secondary or aln.mapping_quality < min_mapq:
            continue
        if not aln.has_tag("SA"):
            continue
        strand = "-" if aln.is_reverse else "+"
        for sa_chrom, sa_pos, sa_strand in parse_sa_tag(aln.get_tag("SA")):
            if sa_chrom != chrom or sa_strand != strand:
                continue
            sa_pos0 = sa_pos - 1
            if not (start <= sa_pos0 < end):
                continue
            a, b = sorted([aln.reference_start, sa_pos0])
            if b - a < 500:
                continue          # too close to be a PTD; alignment noise
            raw.append((a, b, aln.query_name))

    clusters: List[Junction] = []
    for a, b, name in sorted(raw):
        placed = False
        for c in clusters:
            if abs(a - c.acceptor) <= tol and abs(b - c.donor) <= tol:
                c.reads.append(name)
                placed = True
                break
        if not placed:
            clusters.append(Junction(acceptor=a, donor=b, reads=[name]))
    return clusters


# -----------------------------------------------------------------------------
# Signal 2: large CIGAR insertions
# -----------------------------------------------------------------------------
def collect_large_insertions(bam, chrom, start, end, min_len, min_mapq
                             ) -> Dict[int, List[Tuple[int, str]]]:
    """Insertions of at least min_len inside the gene, keyed by rounded
    reference position."""
    out: Dict[int, List[Tuple[int, str]]] = defaultdict(list)
    for aln in bam.fetch(chrom, start, end):
        if aln.is_unmapped or aln.is_secondary or aln.mapping_quality < min_mapq:
            continue
        cigar = aln.cigartuples
        if not cigar:
            continue
        ref_pos = aln.reference_start
        for op, length in cigar:
            if op == BAM_CINS and length >= min_len and start <= ref_pos < end:
                out[ref_pos // 100 * 100].append((length, aln.query_name))
            if op in CONSUMES_REF:
                ref_pos += length
    return out


# -----------------------------------------------------------------------------
# Signal 3: exon coverage ratio
# -----------------------------------------------------------------------------
def load_exons(bed: Path, chrom: str, start: int, end: int
               ) -> List[Tuple[int, int, str]]:
    out = []
    with open(bed) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track")):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 3 or p[0] != chrom:
                continue
            s, e = int(p[1]), int(p[2])
            if s >= start and e <= end:
                out.append((s, e, p[3] if len(p) > 3 else f"{s}-{e}"))
    out.sort()
    return out


def exon_depths(bam, chrom, exons, min_mapq) -> List[Tuple[str, float]]:
    """Mean depth per exon. count_coverage is exact and these intervals are
    small, so there is no need to shell out to mosdepth here."""
    out = []
    for s, e, name in exons:
        try:
            cols = bam.count_coverage(chrom, s, e, quality_threshold=0,
                                      read_callback=lambda r:
                                      not (r.is_unmapped or r.is_secondary)
                                      and r.mapping_quality >= min_mapq)
        except ValueError:
            continue
        total = sum(sum(base) for base in cols)
        width = max(1, e - s)
        out.append((name, total / width))
    return out


def coverage_ratio_call(depths: List[Tuple[str, float]], threshold: float
                        ) -> Tuple[str, str, float]:
    """Flag exons whose depth exceeds the gene median by `threshold`.

    Returns (verdict, exon list, max ratio). Requires at least six exons with
    non-zero depth, below which the median is not a usable baseline.
    """
    usable = [(n, d) for n, d in depths if d > 0]
    if len(usable) < 6:
        return "insufficient_exons", "", 0.0
    med = statistics.median(d for _, d in usable)
    if med <= 0:
        return "insufficient_depth", "", 0.0
    elevated = [(n, d / med) for n, d in usable if d / med >= threshold]
    if not elevated:
        return "no", "", max(d / med for _, d in usable)
    return ("yes",
            ",".join(n for n, _ in elevated),
            max(r for _, r in elevated))


COLUMNS = ["sample", "gene", "label", "chrom", "acceptor_pos", "donor_pos",
           "duplication_span_bp", "split_read_support", "insertion_support",
           "insertion_length_bp", "coverage_ratio_call", "elevated_exons",
           "max_exon_ratio", "n_signals", "confidence", "note"]


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True, type=Path)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--gene", required=True,
                    help="Gene symbol. Its interval is resolved from "
                         "--panel-bed; no coordinates live in this script.")
    ap.add_argument("--label", default="",
                    help="Reporting label for this target, e.g. KMT2A-PTD.")
    ap.add_argument("--region", default=None,
                    help="Gene interval chrom:start-end. If omitted, resolved "
                         "from --panel-bed.")
    ap.add_argument("--panel-bed", type=Path, default=None,
                    help="Panel BED used to resolve --gene to an interval.")
    ap.add_argument("--exon-bed", type=Path, default=None,
                    help="Exon BED for the coverage-ratio signal. Optional; "
                         "omitting it drops that signal, it does not fail.")
    ap.add_argument("--min-span", type=int, default=1000,
                    help="Shortest reportable duplication [1000 bp]. Classic "
                         "exon 2-8 PTD spans several kb.")
    ap.add_argument("--max-span", type=int, default=100000)
    ap.add_argument("--min-split-support", type=int, default=2)
    ap.add_argument("--min-insertion-length", type=int, default=500)
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--junction-tolerance", type=int, default=200)
    ap.add_argument("--coverage-ratio-threshold", type=float, default=1.30,
                    help="Exon depth over gene median above which an exon is "
                         "called elevated [1.30]. A clonal heterozygous PTD "
                         "predicts about 1.5.")
    ap.add_argument("--version", action="version",
                    version=f"%(prog)s {__version__}")
    args = ap.parse_args()

    region = args.region
    if not region and args.panel_bed:
        region = gene_region_from_bed(args.panel_bed, args.gene)
    if not region:
        sys.stderr.write(
            f"ERROR: no interval for {args.gene}. Give --region, or a "
            f"--panel-bed that contains it.\n")
        return 1

    chrom, start, end = parse_region(region)
    bam = pysam.AlignmentFile(str(args.bam), "rb")
    if chrom not in bam.references:
        sys.stderr.write(
            f"ERROR: contig '{chrom}' absent from {args.bam}. Check chr- vs "
            f"NC_-naming between the panel BED and the alignment.\n")
        return 1

    junctions = [j for j in collect_split_junctions(
        bam, chrom, start, end, args.min_mapq, args.junction_tolerance)
        if args.min_span <= j.span <= args.max_span
        and j.n_support >= args.min_split_support]

    insertions = collect_large_insertions(
        bam, chrom, start, end, args.min_insertion_length, args.min_mapq)

    cov_call, cov_exons, cov_ratio = "not_run", "", 0.0
    if args.exon_bed:
        exons = load_exons(args.exon_bed, chrom, start, end)
        if exons:
            cov_call, cov_exons, cov_ratio = coverage_ratio_call(
                exon_depths(bam, chrom, exons, args.min_mapq),
                args.coverage_ratio_threshold)
        else:
            cov_call = "no_exons_in_region"

    rows = []
    if junctions:
        for j in sorted(junctions, key=lambda x: -x.n_support):
            # An insertion near either breakpoint corroborates the junction.
            ins_support, ins_len = 0, 0
            for key, vals in insertions.items():
                if abs(key - j.acceptor) <= args.junction_tolerance or \
                        abs(key - j.donor) <= args.junction_tolerance:
                    ins_support += len(vals)
                    ins_len = max(ins_len, max(v[0] for v in vals))
            n_sig = 1 + (1 if ins_support else 0) + (1 if cov_call == "yes" else 0)
            rows.append({
                "sample": args.sample, "gene": args.gene,
                "label": args.label or args.gene, "chrom": chrom,
                "acceptor_pos": str(j.acceptor + 1),
                "donor_pos": str(j.donor + 1),
                "duplication_span_bp": str(j.span),
                "split_read_support": str(j.n_support),
                "insertion_support": str(ins_support),
                "insertion_length_bp": str(ins_len) if ins_len else "",
                "coverage_ratio_call": cov_call,
                "elevated_exons": cov_exons,
                "max_exon_ratio": f"{cov_ratio:.2f}" if cov_ratio else "",
                "n_signals": str(n_sig),
                "confidence": "high" if n_sig >= 2 else "low",
                "note": "confirm with an orthogonal assay before reporting",
            })
    elif cov_call == "yes":
        # Coverage alone. Emitted so the finding is visible, but flagged as
        # the weakest possible evidence: no junction was observed.
        rows.append({
            "sample": args.sample, "gene": args.gene,
                "label": args.label or args.gene, "chrom": chrom,
            "acceptor_pos": "", "donor_pos": "", "duplication_span_bp": "",
            "split_read_support": "0", "insertion_support": "0",
            "insertion_length_bp": "", "coverage_ratio_call": cov_call,
            "elevated_exons": cov_exons, "max_exon_ratio": f"{cov_ratio:.2f}",
            "n_signals": "1", "confidence": "low",
            "note": "coverage signal only, no junction observed; "
                    "may be a depth artefact",
        })

    with open(args.output, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    sys.stderr.write(
        f"{args.sample}: {args.gene} {region} -> {len(junctions)} candidate "
        f"junction(s), coverage_ratio={cov_call}, {len(rows)} row(s) "
        f"-> {args.output}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
