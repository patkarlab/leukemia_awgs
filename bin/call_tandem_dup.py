#!/usr/bin/env python3
"""
call_tandem_dup.py
==================

Detects short internal tandem duplications from long-read alignments, for any
gene named on the command line.

WHY A DEDICATED CALLER
----------------------
ClairS-TO and Clair3 are tuned for SNVs and short indels. Internal tandem
duplications run from about 15 bp to several hundred, often sit in GC-rich
sequence, and are frequently missed or split into several small indels by
general-purpose callers. For FLT3-ITD in particular, both duplication length
and allelic ratio are ELN 2022 risk inputs, so under-calling changes the risk
group.

Long ONT reads suit this well: a read spanning the whole duplication carries
it as a single CIGAR insertion, so the event is observed rather than inferred
from coverage.

NOTHING GENE-SPECIFIC IS HARDCODED
----------------------------------
The gene is a parameter. Its interval is resolved from the panel BED handed in
at runtime, so the caller cannot drift out of step with the panel and holds no
coordinates of its own. Which genes get called, and with what size and support
thresholds, comes from assets/focal_duplication_targets.tsv.

Hotspot annotation works the same way: an optional BED of named windows is
intersected with each call and its col-4 label reported. A missing or wrong
hotspot BED costs an annotation column, never a call.

METHOD
------
1. Take every alignment overlapping the gene interval.
2. Collect CIGAR insertions of at least --min-len. Insertions are clustered by
   reference position and length, because ONT reads place the same junction a
   few bases apart and vary reported insert length by a few percent.
3. Verify each cluster is a genuine tandem duplication: the inserted sequence
   must align to reference sequence immediately adjacent to the insertion
   point at or above --min-tandem-identity. This is what separates a tandem
   duplication from a random insertion or an alignment artefact.
4. Report allelic ratio as supporting reads over reads spanning the insertion
   point. This is a read-level ratio; the clinical mutant:wild-type ratio from
   fragment-analysis PCR is a different measurement, and the two are not
   interchangeable without local validation.

Dependencies: pysam (present in the awgs_sv conda env).
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

try:
    import pysam
except ImportError:
    sys.stderr.write("ERROR: pysam is required. Activate the awgs_sv env.\n")
    sys.exit(1)

__version__ = "0.1.0"

# CIGAR operation codes as used by pysam.
BAM_CMATCH, BAM_CINS, BAM_CDEL, BAM_CREF_SKIP = 0, 1, 2, 3
BAM_CSOFT_CLIP, BAM_CHARD_CLIP, BAM_CPAD = 4, 5, 6
BAM_CEQUAL, BAM_CDIFF = 7, 8

CONSUMES_QUERY = {BAM_CMATCH, BAM_CINS, BAM_CSOFT_CLIP, BAM_CEQUAL, BAM_CDIFF}
CONSUMES_REF = {BAM_CMATCH, BAM_CDEL, BAM_CREF_SKIP, BAM_CEQUAL, BAM_CDIFF}


@dataclass
class RawInsertion:
    ref_pos: int          # reference position immediately before the insert
    length: int
    seq: str
    read_name: str


@dataclass
class ItdCluster:
    ref_pos: int
    length: int
    reads: List[str] = field(default_factory=list)
    seqs: List[str] = field(default_factory=list)

    @property
    def n_support(self) -> int:
        return len(self.reads)


_SUFFIX = {"LOCUS", "CLUSTER", "ENHANCER", "REGION", "INTERVAL"}


def gene_region_from_bed(bed: Path, gene: str) -> Optional[str]:
    """Resolve a gene symbol to chrom:start-end using the panel BED col-4 label.

    Labels may be compound where regions merged during panel construction
    ('TAL1/STIL'), so the label is tokenised rather than compared as a string.
    Keeping the lookup here means no genomic coordinate is written into code.
    """
    want = gene.strip().upper()
    for line in bed.read_text().splitlines():
        if not line.strip() or line.startswith(("#", "track")):
            continue
        f = line.split("\t")
        if len(f) < 4:
            continue
        toks = {x for x in f[3].upper().replace("+", "/").replace("_", "/").split("/")
                if x and x not in _SUFFIX}
        if want in toks:
            return f"{f[0]}:{int(f[1]) + 1}-{f[2]}"
    return None


def load_hotspots(bed: Optional[Path]) -> List[Tuple[str, int, int, str]]:
    """Named annotation windows: (chrom, start, end, label). Optional."""
    if not bed or not bed.exists():
        return []
    out = []
    for line in bed.read_text().splitlines():
        if not line.strip() or line.startswith(("#", "track")):
            continue
        f = line.split("\t")
        if len(f) >= 3:
            out.append((f[0], int(f[1]), int(f[2]),
                        f[3] if len(f) > 3 else f"{f[0]}:{f[1]}-{f[2]}"))
    return out


def hotspot_label(hotspots, chrom: str, pos: int) -> str:
    for hc, hs, he, label in hotspots:
        if hc == chrom and hs <= pos < he:
            return label
    return ""


def parse_region(region: str) -> Tuple[str, int, int]:
    """Parse 'chrN:start-end' into (chrom, start, end), 0-based start."""
    try:
        chrom, span = region.rsplit(":", 1)
        start, end = span.replace(",", "").split("-")
        return chrom, int(start) - 1, int(end)
    except ValueError:
        sys.stderr.write(f"ERROR: cannot parse region '{region}'. "
                         f"Expected chrom:start-end.\n")
        sys.exit(1)


def collect_insertions(bam: "pysam.AlignmentFile", chrom: str, start: int,
                       end: int, min_len: int, min_mapq: int
                       ) -> Tuple[List[RawInsertion], int]:
    """Walk each alignment's CIGAR and record insertions inside the window."""
    out: List[RawInsertion] = []
    n_reads = 0
    for aln in bam.fetch(chrom, start, end):
        if aln.is_unmapped or aln.is_secondary:
            continue
        if aln.mapping_quality < min_mapq:
            continue
        n_reads += 1
        cigar = aln.cigartuples
        if not cigar:
            continue
        ref_pos = aln.reference_start
        query_pos = 0
        seq = aln.query_sequence or ""
        for op, length in cigar:
            if op == BAM_CINS:
                if length >= min_len and start <= ref_pos < end:
                    out.append(RawInsertion(
                        ref_pos=ref_pos, length=length,
                        seq=seq[query_pos:query_pos + length] if seq else "",
                        read_name=aln.query_name))
            if op in CONSUMES_QUERY:
                query_pos += length
            if op in CONSUMES_REF:
                ref_pos += length
    return out, n_reads


def cluster_insertions(raws: List[RawInsertion], pos_tol: int,
                       len_tol_frac: float) -> List[ItdCluster]:
    """Group insertions that plausibly represent the same event.

    ONT alignments place an identical junction a few bases apart and vary the
    reported insert length by a few percent, so exact grouping would split one
    ITD into many singletons.
    """
    clusters: List[ItdCluster] = []
    for r in sorted(raws, key=lambda x: (x.ref_pos, x.length)):
        placed = False
        for c in clusters:
            len_tol = max(2, int(c.length * len_tol_frac))
            if abs(r.ref_pos - c.ref_pos) <= pos_tol and \
                    abs(r.length - c.length) <= len_tol:
                # Running mean keeps the cluster centred as reads accumulate.
                n = c.n_support
                c.ref_pos = int((c.ref_pos * n + r.ref_pos) / (n + 1))
                c.length = int((c.length * n + r.length) / (n + 1))
                c.reads.append(r.read_name)
                if r.seq:
                    c.seqs.append(r.seq)
                placed = True
                break
        if not placed:
            clusters.append(ItdCluster(ref_pos=r.ref_pos, length=r.length,
                                       reads=[r.read_name],
                                       seqs=[r.seq] if r.seq else []))
    return clusters


def consensus_seq(seqs: List[str]) -> str:
    """Longest common length-normalised representative of the cluster.

    A true consensus would need a multiple alignment; the median-length read
    is adequate for the tandem check and keeps this dependency-free.
    """
    if not seqs:
        return ""
    return sorted(seqs, key=len)[len(seqs) // 2]


def identity(a: str, b: str) -> float:
    """Ungapped identity over the overlapping prefix. ONT insert sequences
    carry indel noise, so this is a permissive check by design; it exists to
    exclude non-tandem insertions, not to polish the sequence."""
    if not a or not b:
        return 0.0
    n = min(len(a), len(b))
    if n == 0:
        return 0.0
    matches = sum(1 for i in range(n) if a[i].upper() == b[i].upper())
    return matches / n


def verify_tandem(cluster: ItdCluster, fasta: "pysam.FastaFile", chrom: str,
                  min_identity: float) -> Tuple[bool, float, str]:
    """Check the inserted sequence duplicates adjacent reference sequence.

    An ITD duplicates a stretch immediately upstream or downstream of the
    insertion point, so the insert should align to one or the other. Both
    orientations are tested and the better score reported.
    """
    ins = consensus_seq(cluster.seqs)
    if not ins:
        return False, 0.0, "no_sequence"
    n = len(ins)
    try:
        upstream = fasta.fetch(chrom, max(0, cluster.ref_pos - n),
                               cluster.ref_pos)
        downstream = fasta.fetch(chrom, cluster.ref_pos,
                                 cluster.ref_pos + n)
    except (ValueError, KeyError):
        return False, 0.0, "reference_fetch_failed"

    id_up = identity(ins, upstream)
    id_down = identity(ins, downstream)
    if id_up >= id_down:
        best, which = id_up, "upstream"
    else:
        best, which = id_down, "downstream"
    return best >= min_identity, best, which


def spanning_depth(bam: "pysam.AlignmentFile", chrom: str, pos: int,
                   min_mapq: int) -> int:
    """Reads whose alignment spans the insertion point, ITD-carrying or not.

    This is the denominator for the allelic ratio. Reads that merely start or
    end near the junction are excluded because they cannot report on it.
    """
    n = 0
    for aln in bam.fetch(chrom, max(0, pos - 1), pos + 1):
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
            continue
        if aln.mapping_quality < min_mapq:
            continue
        if aln.reference_start <= pos < aln.reference_end:
            n += 1
    return n


COLUMNS = ["sample", "gene", "label", "chrom", "ref_pos", "dup_length_bp",
           "in_frame", "support_reads", "spanning_reads", "allelic_ratio",
           "tandem_verified", "tandem_identity", "duplicated_side",
           "hotspot", "callers", "note"]


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True, type=Path)
    ap.add_argument("--fasta", required=True, type=Path,
                    help="Reference FASTA matching the BAM, indexed (.fai).")
    ap.add_argument("--gene", required=True,
                    help="Gene symbol. Its interval is resolved from "
                         "--panel-bed; no coordinates live in this script.")
    ap.add_argument("--panel-bed", type=Path, default=None,
                    help="Panel BED whose col-4 labels resolve --gene.")
    ap.add_argument("--region", default=None,
                    help="Explicit chrom:start-end override. Normally omitted "
                         "so the panel BED stays the single source of truth.")
    ap.add_argument("--hotspot-bed", type=Path, default=None,
                    help="Optional BED of named annotation windows. A call "
                         "landing in one reports that window's col-4 label. "
                         "Annotation only; never filters, so a missing or "
                         "stale hotspot BED cannot cause a false negative.")
    ap.add_argument("--label", default="",
                    help="Reporting label for this target, e.g. FLT3-ITD.")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--min-len", type=int, default=15,
                    help="Shortest reportable duplication [15].")
    ap.add_argument("--max-len", type=int, default=500)
    ap.add_argument("--min-support", type=int, default=3,
                    help="Minimum supporting reads to emit a call [3].")
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--pos-tolerance", type=int, default=20,
                    help="Clustering tolerance on insertion position [20 bp].")
    ap.add_argument("--len-tolerance-frac", type=float, default=0.15,
                    help="Clustering tolerance on insert length, as a "
                         "fraction of cluster length [0.15].")
    ap.add_argument("--min-tandem-identity", type=float, default=0.70,
                    help="Identity to adjacent reference required to call a "
                         "cluster a tandem duplication [0.70]. ONT insert "
                         "sequence carries indel noise, so this is permissive.")
    ap.add_argument("--emit-unverified", action="store_true",
                    help="Also emit clusters that failed the tandem check, "
                         "flagged tandem_verified=no. Audit aid.")
    ap.add_argument("--version", action="version",
                    version=f"%(prog)s {__version__}")
    args = ap.parse_args()

    region = args.region
    if not region and args.panel_bed:
        region = gene_region_from_bed(args.panel_bed, args.gene)
    if not region:
        sys.stderr.write(
            f"ERROR: no interval for {args.gene}. Supply a --panel-bed that "
            f"carries it, or an explicit --region.\n")
        return 1
    chrom, start, end = parse_region(region)
    hotspots = load_hotspots(args.hotspot_bed)

    bam = pysam.AlignmentFile(str(args.bam), "rb")
    if chrom not in bam.references:
        sys.stderr.write(
            f"ERROR: contig '{chrom}' absent from {args.bam}. "
            f"Check that --region naming matches the BAM "
            f"(chr-named vs NC_-named).\n")
        return 1
    fasta = pysam.FastaFile(str(args.fasta))

    raws, n_window_reads = collect_insertions(
        bam, chrom, start, end, args.min_len, args.min_mapq)
    clusters = cluster_insertions(raws, args.pos_tolerance,
                                  args.len_tolerance_frac)

    rows = []
    for c in sorted(clusters, key=lambda x: -x.n_support):
        if c.n_support < args.min_support:
            continue
        if not (args.min_len <= c.length <= args.max_len):
            continue
        verified, ident, side = verify_tandem(
            c, fasta, chrom, args.min_tandem_identity)
        if not verified and not args.emit_unverified:
            continue
        depth = spanning_depth(bam, chrom, c.ref_pos, args.min_mapq)
        ratio = (c.n_support / depth) if depth else 0.0
        hs = hotspot_label(hotspots, chrom, c.ref_pos)
        rows.append({
            "sample": args.sample,
            "gene": args.gene,
            "label": args.label or args.gene,
            "chrom": chrom,
            "ref_pos": str(c.ref_pos + 1),          # 1-based for reporting
            "dup_length_bp": str(c.length),
            "in_frame": "yes" if c.length % 3 == 0 else "no",
            "support_reads": str(c.n_support),
            "spanning_reads": str(depth),
            "allelic_ratio": f"{ratio:.3f}",
            "tandem_verified": "yes" if verified else "no",
            "tandem_identity": f"{ident:.2f}",
            "duplicated_side": side,
            "hotspot": hs,
            "callers": "call_tandem_dup",
            "note": "" if verified else "failed tandem check; not a called duplication",
        })

    with open(args.output, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    n_called = sum(1 for r in rows if r["tandem_verified"] == "yes")
    sys.stderr.write(
        f"{args.sample}: {args.gene} {region}, {n_window_reads} reads, "
        f"{len(raws)} raw insertions >= {args.min_len} bp, "
        f"{len(clusters)} clusters, {n_called} verified TD -> {args.output}\n")
    if n_window_reads == 0:
        sys.stderr.write(
            f"WARNING: no reads at {region}. Either the region is "
            f"off-panel for this run or the contig naming is wrong.\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
