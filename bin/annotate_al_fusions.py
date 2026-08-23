#!/usr/bin/env python3
"""
annotate_al_fusions.py
======================

Annotates a SURVIVOR-merged SV VCF against an acute-leukaemia panel.

Ported from annotate_mm_translocations.py (mm-awgs-nextflow) with three
changes that acute leukaemia requires and myeloma did not:

  1. PROMISCUOUS ANCHORS. KMT2A has >100 described partners, NUP98 >30, and
     any RARA junction is APL until proven otherwise. A dictionary of named
     pairs cannot cover these. Anchors listed in the anchor table are
     reported whenever they appear on either side of a junction, whatever
     the partner is, including when the partner is off-panel.

  2. SINGLE-SIDED DICTIONARY MATCHING. Some defining partners are not on the
     panel by design: PBX1 (TCF3::PBX1) and DUX4 (IGH::DUX4) are the two
     that matter most. When one side is on-panel and the other resolves only
     to a cytoband, a dictionary row carrying a partner_b_band that equals
     that cytoband still names the event, flagged as a partial match.

  3. DISEASE SCOPE. The dictionary is shared between the AML and ALL panels;
     rows are filtered to the panel in use so an ALL run does not report AML
     entities and vice versa. Rows marked BOTH always apply.

The script holds no biological priors of its own. Everything reportable
comes from --dictionary and --anchors.

Usage:
  annotate_al_fusions.py \
      --vcf merged.vcf.gz \
      --panel-bed AML_panel_t2t_chr.bed \
      --cytoband-bed chm13v2.0_cytobands_allchrs.bed \
      --dictionary al_fusion_dictionary.tsv \
      --anchors al_fusion_anchors.tsv \
      --panel AML \
      --sample SAMPLE_ID \
      --output SAMPLE_ID.al_annotated.tsv
"""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

__version__ = "0.1.0"


# -----------------------------------------------------------------------------
# Panel and cytoband tables
# -----------------------------------------------------------------------------
@dataclass
class PanelRegion:
    chrom: str
    start: int
    end: int
    name: str

    def contains(self, chrom: str, pos: int) -> bool:
        return chrom == self.chrom and self.start <= pos < self.end


def load_panel(bed_path: Path) -> List[PanelRegion]:
    out: List[PanelRegion] = []
    with open(bed_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith(("#", "track")):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            name = parts[3] if len(parts) >= 4 else f"{chrom}:{start}-{end}"
            out.append(PanelRegion(chrom, start, end, name))
    if not out:
        sys.stderr.write(f"ERROR: no panel regions parsed from {bed_path}\n")
        sys.exit(1)
    return out


@dataclass
class CytobandTable:
    bands: Dict[str, List[Tuple[int, int, str]]]

    def band_for(self, chrom: Optional[str], pos: Optional[int]) -> Optional[str]:
        if chrom is None or pos is None:
            return None
        for start, end, name in self.bands.get(chrom, []):
            if start <= pos < end:
                return name
        return None


def load_cytobands(bed_path: Path) -> CytobandTable:
    bands: Dict[str, List[Tuple[int, int, str]]] = {}
    with open(bed_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith(("#", "track")):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            bands.setdefault(parts[0], []).append(
                (int(parts[1]), int(parts[2]), parts[3]))
    for chrom in bands:
        bands[chrom].sort(key=lambda t: t[0])
    if not bands:
        sys.stderr.write(f"ERROR: no cytobands parsed from {bed_path}\n")
        sys.exit(1)
    return CytobandTable(bands)


# -----------------------------------------------------------------------------
# Token normalisation
# -----------------------------------------------------------------------------
# Words that appear in region or partner labels but are not gene symbols.
# Dropped so that IGH_locus matches IGH, and HOXA_cluster matches HOXA.
_SUFFIX_TOKENS = {"LOCUS", "CLUSTER", "ENHANCER", "REGION", "INTERVAL"}


def norm_tokens(name: str) -> frozenset:
    """Reduce a partner or region label to a set of gene tokens.

    Splits on '/', '_' and '+', uppercases, and drops region-suffix words.
    'IGH_locus'      -> {IGH}
    'TAL1/STIL'      -> {TAL1, STIL}
    'PAR1_CRLF2_P2RY8' -> {PAR1, CRLF2, P2RY8}
    Two labels refer to the same locus when their token sets intersect.

    NKX2-1 and NKX2-5 keep their hyphen deliberately: splitting on '-' would
    collapse both to {NKX2, 1} / {NKX2, 5} and make them cross-match.
    """
    raw = str(name).strip().upper().replace("+", "/").replace("_", "/")
    return frozenset(t for t in raw.split("/") if t and t not in _SUFFIX_TOKENS)


def strip_chr(chrom: str) -> str:
    return chrom[3:] if chrom.startswith("chr") else chrom


# -----------------------------------------------------------------------------
# Dictionary and anchors
# -----------------------------------------------------------------------------
@dataclass
class DictEntry:
    tok_a: frozenset
    tok_b: frozenset
    band_b: str
    row: dict


def load_dictionary(path: Path, panel: str) -> List[DictEntry]:
    """Load named fusion pairs, keeping rows whose disease matches the panel.

    A missing dictionary is fatal here (unlike the MM version): running an
    acute-leukaemia panel with no priors would silently report every junction
    as unnamed, which is worse than failing at launch.
    """
    if not path.exists():
        sys.stderr.write(f"ERROR: dictionary not found: {path}\n")
        sys.exit(1)
    out: List[DictEntry] = []
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            disease = (row.get("disease") or "BOTH").strip().upper()
            if disease not in ("BOTH", panel.upper()):
                continue
            a = (row.get("partner_a") or "").strip()
            b = (row.get("partner_b") or "").strip()
            if not a or not b:
                continue
            out.append(DictEntry(norm_tokens(a), norm_tokens(b),
                                 (row.get("partner_b_band") or "").strip(), row))
    return out


def load_anchors(path: Path, panel: str) -> Dict[str, dict]:
    """Map gene token -> anchor row, for anchors in scope for this panel."""
    out: Dict[str, dict] = {}
    if not path.exists():
        sys.stderr.write(f"WARNING: anchor table not found: {path}\n")
        return out
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            disease = (row.get("disease") or "BOTH").strip().upper()
            if disease not in ("BOTH", panel.upper()):
                continue
            for tok in norm_tokens(row.get("anchor") or ""):
                out[tok] = row
    return out


def _band_levels(band: Optional[str]) -> List[str]:
    """Progressively coarser forms of a cytoband label, most precise first.

    '1q23.3' -> ['1q23.3', '1q23', '1q']

    Exact sub-band equality is too brittle to rely on: T2T-CHM13 band
    boundaries differ from GRCh38's, and a breakpoint in an adjacent sub-band
    would otherwise silently lose a defining entity such as TCF3::PBX1. The
    levels are compared in order so the match records how precise it was.
    """
    if not band:
        return []
    out = [band]
    if "." in band:
        out.append(band.split(".", 1)[0])
    for i, ch in enumerate(out[-1]):
        if ch in "pq":
            arm = out[-1][:i + 1]
            if arm != out[-1]:
                out.append(arm)
            break
    return out


def _band_match(observed: Optional[str], expected: str) -> Optional[str]:
    """Return the match quality if observed and expected agree at any level."""
    if not observed or not expected:
        return None
    obs_levels = _band_levels(observed)
    exp_levels = _band_levels(expected)
    if not obs_levels or not exp_levels:
        return None
    # Walk coarse -> fine; remember the finest level that agrees and stop
    # at the first shared level that disagrees. 4q25 vs 4q35 agree on the
    # arm (level 0) but differ at the major band (level 1): that is a
    # contradiction, not a partial match.
    depth = -1
    for i in range(min(len(obs_levels), len(exp_levels))):
        if obs_levels[i] == exp_levels[i]:
            depth = i
        else:
            if depth < 0:
                return None
            break
    if depth >= 1:
        return "partial_band"   # major band (or finer) agrees
    if depth == 0:
        return "partial_arm"    # arm only: a lead, not a call
    return None


def dictionary_lookup(dictionary: List[DictEntry], label_a: str, label_b: str,
                      band_a: Optional[str], band_b: Optional[str]
                      ) -> Tuple[Optional[dict], str]:
    """Return (row, match_quality) for an unordered pair.

    'full'          both sides matched on gene identity
    'partial_band'  one side matched a gene, the other matched the
                    dictionary's expected cytoband for its off-panel partner
    'partial_arm'   as above but agreeing only at chromosome-arm level;
                    treat as a lead, not a call
    Returns (None, '') when nothing matches.
    """
    ta, tb = norm_tokens(label_a), norm_tokens(label_b)

    for e in dictionary:
        if (ta & e.tok_a and tb & e.tok_b) or (ta & e.tok_b and tb & e.tok_a):
            return e.row, "full"

    # One side on-panel, other side resolvable only to a band. The
    # dictionary's partner_b_band is what makes these nameable at all.
    best = None
    for e in dictionary:
        if not e.band_b:
            continue
        for gene_tok, other_band in ((ta, band_b), (tb, band_a)):
            if not (gene_tok & e.tok_a):
                continue
            q = _band_match(other_band, e.band_b)
            if q == "partial_band":
                return e.row, q
            if q and best is None:
                best = (e.row, q)
    return best if best else (None, "")


def anchor_hits(anchors: Dict[str, dict], label_a: str, label_b: str) -> List[dict]:
    """Anchor rows triggered by either side of a junction.

    A self-pair returns nothing. An anchor is a claim about partners: the
    gene rearranges with many of them, so any partner is worth surfacing.
    A gene joined to itself has no partner, and the premise does not apply.
    Without this guard, every intragenic indel at an anchor inherits
    reportability from a rule that was never about it: on one validation
    sample 513 of 589 reportable rows were self-pairs, 335 of them V(D)J
    and somatic hypermutation products inside IGH, IGK, IGL and the TCR
    loci, which recombine physiologically and are not lesions.

    Dictionary matching is untouched. No row in the dictionary pairs a gene
    with itself or one locus with another, so no self-pair can carry a tier;
    if one ever could, the caller's `hit or hits` keeps that path open.

    Ordinary-gene self-pairs are real intragenic events and stay in the
    table as reportable=no. They are not fusions and should not be graded
    by a partner-pair dictionary; representing them properly is a separate
    change.
    """
    if label_a and label_b and label_a == label_b:
        return []
    hits, seen = [], set()
    for label in (label_a, label_b):
        for tok in norm_tokens(label):
            row = anchors.get(tok)
            if row and row["anchor"] not in seen:
                seen.add(row["anchor"])
                hits.append(row)
    return hits


# -----------------------------------------------------------------------------
# VCF parsing (unchanged from the MM pipeline)
# -----------------------------------------------------------------------------
@dataclass
class SvRecord:
    chrom: str
    pos: int
    sv_id: str
    sv_type: str
    mate_chrom: Optional[str]
    mate_pos: Optional[int]
    filt: str
    info: Dict[str, str]
    callers: List[str] = field(default_factory=list)
    support: str = ""


def parse_info(info_field: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for part in info_field.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
        else:
            out[part] = "True"
    return out


def parse_bnd_alt(alt: str) -> Tuple[Optional[str], Optional[int]]:
    """Parse a BND ALT such as N]chr11:69500000] or [chr9:133600000[N."""
    for bracket in ("]", "["):
        if bracket in alt:
            try:
                inner = alt.split(bracket)[1].split(bracket)[0]
                chrom, pos = inner.rsplit(":", 1)
                return chrom, int(pos)
            except (IndexError, ValueError):
                continue
    return None, None


def infer_callers(info: Dict[str, str]) -> List[str]:
    """Decode SURVIVOR's SUPP_VEC. Bit order is the order of input VCFs,
    which in this pipeline is [Sniffles, CuteSV, Severus]."""
    order = ["Sniffles", "CuteSV", "Severus"]
    return [n for bit, n in zip(info.get("SUPP_VEC", ""), order) if bit == "1"]


def support_reads_from(info: Dict[str, str], fmt: str,
                       sample_cols: List[str]) -> str:
    """Variant-supporting read count. The merged VCF carries one sample column
    per input caller; only the caller that found the junction holds a real
    FORMAT/DV, so take the maximum across columns, then fall back to
    caller-specific INFO tags."""
    if fmt:
        keys = fmt.split(":")
        best = None
        for col in sample_cols:
            if not col or col in (".", "./."):
                continue
            f = dict(zip(keys, col.split(":")))
            dv = f.get("DV", "")
            if dv.isdigit():
                best = int(dv) if best is None else max(best, int(dv))
        if best is not None:
            return str(best)
    for tag in ("RE", "SUPPORT", "SR"):
        v = str(info.get(tag, ""))
        if v.isdigit():
            return v
    first = str(info.get("SUPP_READS", "")).split(":")[0]
    return first if first.isdigit() else ""


def open_vcf(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def parse_vcf(vcf_path: Path) -> List[SvRecord]:
    out: List[SvRecord] = []
    with open_vcf(vcf_path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            chrom, pos, sv_id, _ref, alt, _qual, filt, info_field = cols[:8]
            info = parse_info(info_field)
            sv_type = info.get("SVTYPE", "")
            mate_chrom, mate_pos = None, None

            if sv_type in ("BND", "TRA"):
                mate_chrom, mate_pos = parse_bnd_alt(alt)
                if mate_chrom is None:
                    chr2, end = info.get("CHR2"), info.get("END")
                    if chr2 and end:
                        try:
                            mate_chrom, mate_pos = chr2, int(end)
                        except ValueError:
                            pass
            elif sv_type in ("DEL", "DUP", "INV", "INS"):
                end = info.get("END")
                if end:
                    try:
                        mate_chrom, mate_pos = chrom, int(end)
                    except ValueError:
                        pass

            fmt = cols[8] if len(cols) > 8 else ""
            sample_cols = cols[9:] if len(cols) > 9 else []
            out.append(SvRecord(
                chrom=chrom, pos=int(pos), sv_id=sv_id, sv_type=sv_type,
                mate_chrom=mate_chrom, mate_pos=mate_pos, filt=filt, info=info,
                callers=infer_callers(info),
                support=support_reads_from(info, fmt, sample_cols)))
    return out


# -----------------------------------------------------------------------------
# Annotation
# -----------------------------------------------------------------------------
def region_for(chrom, pos, panel) -> Optional[PanelRegion]:
    if chrom is None or pos is None:
        return None
    for r in panel:
        if r.contains(chrom, pos):
            return r
    return None


def gene_for(chrom, pos, model):
    """Tightest gene-model feature containing this coordinate, or None.

    Tightest wins because gene models legitimately overlap (a gene inside
    another's intron); the smaller feature is the more specific answer."""
    if not model or chrom is None or pos is None:
        return None
    best = None
    for reg in model:
        if reg.chrom == chrom and reg.start <= pos < reg.end:
            if best is None or (reg.end - reg.start) < (best.end - best.start):
                best = reg
    return best.name if best else None


def dist_to_gene(chrom, pos, name, model):
    """Bases from this coordinate to the named gene's body; 0 if inside.

    None when the label is not a gene in the model - a named interval such as
    IGH_locus, or a cytoband - because there is no body to measure from.
    """
    if not model or chrom is None or pos is None or not name:
        return None
    best = None
    for reg in model:
        if reg.chrom == chrom and reg.name == name:
            d = max(reg.start - pos, pos - reg.end, 0)
            if best is None or d < best:
                best = d
    return best


def characterize_side(chrom, pos, region, cytobands, gene_model=None):
    """Return (label, source, band) for one breakpoint side.

    Resolution order: gene model, panel interval label, cytoband, coordinate.
    The gene model comes first because a panel interval may be a merged
    compound ("PAX5/ZCCHC7") or a named locus, neither of which names the gene
    a breakpoint actually fell in. Panel membership is decided separately, by
    region_for against the panel BED; this only supplies the label."""
    band = cytobands.band_for(chrom, pos)
    band_label = f"{strip_chr(chrom)}{band}" if band and chrom else None
    gene = gene_for(chrom, pos, gene_model)
    if gene:
        return gene, "gene_model", band_label
    if region is not None:
        return region.name, "panel", band_label
    if band_label:
        return band_label, "cytoband", band_label
    if chrom is not None and pos is not None:
        return f"{chrom}:{pos / 1e6:.1f}Mb", "coordinate", None
    return "OFF_PANEL", "coordinate", None


def annotate(records, panel, dictionary, anchors, sample, cytobands, panel_name,
             gene_model=None):
    out = []
    for r in records:
        side_a = region_for(r.chrom, r.pos, panel)
        side_b = region_for(r.mate_chrom, r.mate_pos, panel)
        if side_a is None and side_b is None:
            continue

        gene_a, src_a, band_a = characterize_side(
            r.chrom, r.pos, side_a, cytobands, gene_model)
        gene_b, src_b, band_b = characterize_side(
            r.mate_chrom, r.mate_pos, side_b, cytobands, gene_model)

        dist_a = dist_to_gene(r.chrom, r.pos, gene_a, gene_model)
        dist_b = dist_to_gene(r.mate_chrom, r.mate_pos, gene_b, gene_model)

        hit, quality = dictionary_lookup(dictionary, gene_a, gene_b, band_a, band_b)

        # Span guard. A cryptic-deletion pair whose two partners share one
        # panel interval matches on tokens alone, so any small intra-interval
        # event is stamped with a defining entity (observed: ~75 false CRLF2
        # calls across PAR1; a 36 bp del labelled del(1)(p33p33)). When the
        # dictionary row declares min_span_bp and the event is
        # intrachromosomal with both positions known, an undersized span
        # demotes the match: the record stays, the entity claim goes.
        span_note = ""
        if hit and quality == "full":
            _ms = (hit.get("min_span_bp") or "").strip()
            if (_ms and r.mate_chrom and r.mate_chrom == r.chrom
                    and r.mate_pos is not None):
                _span = abs(r.mate_pos - r.pos)
                if _span < int(_ms):
                    span_note = (f"span {_span}bp < min {_ms}bp for "
                                 f"{hit.get('name','pair')}; entity removed")
                    hit, quality = None, "below_span"

        # An anchor names a gene, so the breakend should be in that gene,
        # not merely inside the panel interval surrounding it. 22.8 Mb of the
        # 41.6 Mb panel is flank rather than gene body, and on one sample 326
        # of 431 gene-named breakends sat more than 20 kb from the gene they
        # were named after. The consequence was concrete: an ART1 germline CNV
        # reported as a NUP98 deletion in all three validation samples,
        # because the NUP98 panel interval begins 50 kb before the gene and
        # covers ART1 entirely.
        #
        # A named interval (IGH_locus, HOXA_cluster, PAR1_CRLF2_P2RY8) has no
        # gene body to measure from, so dist is None and it is exempt: those
        # labels describe the whole interval accurately, and IGH must still be
        # able to fire an anchor for IGH::CRLF2.
        anchor_a = gene_a if dist_a in (0, None) else None
        anchor_b = gene_b if dist_b in (0, None) else None
        hits = anchor_hits(anchors, anchor_a, anchor_b)

        # An event is reportable when the dictionary names it, or when it
        # touches a promiscuous anchor. Everything else is emitted too but
        # carries reportable=no, so the on-panel callset stays auditable.
        reportable = "yes" if (hit or hits) else "no"

        out.append({
            "sample":          sample,
            "panel":           panel_name,
            "sv_id":           r.sv_id,
            "sv_type":         r.sv_type,
            "filter":          r.filt,
            "chrom_a":         r.chrom,
            "pos_a":           str(r.pos),
            "gene_a":          gene_a,
            "chrom_b":         r.mate_chrom or "",
            "pos_b":           str(r.mate_pos) if r.mate_pos is not None else "",
            "gene_b":          gene_b,
            "gene_a_dist":     "" if dist_a is None else str(dist_a),
            "gene_b_dist":      "" if dist_b is None else str(dist_b),
            "gene_a_source":   src_a,
            "gene_b_source":   src_b,
            "band_a":          band_a or "",
            "band_b":          band_b or "",
            "known_pair":      hit.get("name", "") if hit else "",
            "entity":          hit.get("entity", "") if hit else "",
            "tier":            ((hit.get("tier", "") if quality == "full"
                                 else (hit.get("tier", "") + "?" if hit.get("tier") else ""))
                                if hit else ""),
            "known_freq":      hit.get("frequency", "") if hit else "",
            "match_quality":   quality,
            "anchor":          ",".join(h["anchor"] for h in hits),
            "anchor_class":    ",".join(h["anchor_class"] for h in hits),
            "reportable":      reportable,
            "dict_notes":      (hit.get("notes", "") if hit else span_note),
            "callers":         ",".join(r.callers) or "unknown",
            "n_callers":       str(len(r.callers)),
            "supp_vec":        r.info.get("SUPP_VEC", ""),
            "support_reads":   r.support,
        })
    return out


COLUMNS = [
    "sample", "panel", "sv_id", "sv_type", "filter",
    "chrom_a", "pos_a", "gene_a", "chrom_b", "pos_b", "gene_b",
    "gene_a_source", "gene_b_source", "gene_a_dist", "gene_b_dist",
    "band_a", "band_b",
    "known_pair", "entity", "tier", "known_freq", "match_quality",
    "anchor", "anchor_class", "reportable", "dict_notes",
    "callers", "n_callers", "supp_vec", "support_reads",
]


def load_excluded_junctions(path):
    """Junctions to drop, keyed by chromosome pair, both orientations stored.

    Each entry was observed at base-identical coordinates in unrelated
    patients. Somatic breakpoints do not recur to the nucleotide across
    individuals: repair at a real junction is imprecise, so two patients
    sharing a rearrangement share the intron, not the base.

    Recurrence alone is not the criterion. Acute leukaemia's defining events
    are recurrent by definition, and a KMT2A or MECOM rearrangement in several
    patients is a finding. Coordinate identity is a different claim.
    """
    out = {}
    if not path or not Path(str(path)).is_file():
        return out
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            ca, pa, cb, pb, tol = f[0], int(f[1]), f[2], int(f[3]), int(f[4])
            note = f[6] if len(f) > 6 else ""
            out.setdefault((ca, cb), []).append((pa, pb, tol, note))
            out.setdefault((cb, ca), []).append((pb, pa, tol, note))
    return out


def excluded_reason(excl, row):
    """Why this junction is excluded, or None.

    A dictionary-named pair or a defining-tier entity is never excluded: the
    list is a noise filter, and an entity the dictionary recognises has cleared
    a higher bar than coordinate recurrence can overturn.
    """
    if not excl:
        return None
    if str(row.get("known_pair") or "").strip():
        return None
    if str(row.get("tier") or "").strip().lower() == "defining":
        return None
    for pa, pb, tol, note in excl.get((row.get("chrom_a"), row.get("chrom_b")), []):
        try:
            if abs(int(row["pos_a"]) - pa) <= tol and abs(int(row["pos_b"]) - pb) <= tol:
                return note or "base-identical junction seen in unrelated patients"
        except (TypeError, ValueError, KeyError):
            continue
    return None


def load_gene_model(bed_path):
    """One feature per gene, unmerged and unflanked. Optional: absent means
    naming falls back to the panel label, which is the pre-gene-model
    behaviour. Unlike load_panel this does not exit on an empty file, since
    the model is not required for the run to be valid."""
    if bed_path is None:
        return []
    out = []
    with open(bed_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith(("#", "track")):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            out.append(PanelRegion(parts[0], int(parts[1]), int(parts[2]),
                                   parts[3]))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", required=True, type=Path)
    ap.add_argument("--panel-bed", required=True, type=Path)
    ap.add_argument("--gene-model", type=Path, default=None,
                    help="One feature per gene, unmerged and unflanked. Used "
                         "for naming a breakpoint's gene; panel membership "
                         "still comes from --panel-bed.")
    ap.add_argument("--cytoband-bed", required=True, type=Path)
    ap.add_argument("--dictionary", required=True, type=Path)
    ap.add_argument("--anchors", required=True, type=Path)
    ap.add_argument("--excluded-junctions", default=None,
                    help="TSV of recurrent artefact junctions to drop. Matched "
                         "on both breakpoints within --excluded-tol, in either "
                         "orientation. A dictionary-named pair is never "
                         "dropped.")
    ap.add_argument("--excluded-tol", type=int, default=50,
                    help="Bases either breakpoint may differ from a listed "
                         "artefact and still match [50].")
    ap.add_argument("--panel", required=True, choices=["AML", "ALL"])
    ap.add_argument("--sample", required=True)
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--version", action="version",
                    version=f"%(prog)s {__version__}")
    args = ap.parse_args()

    panel = load_panel(args.panel_bed)
    gene_model = load_gene_model(args.gene_model)
    if gene_model:
        sys.stderr.write(f"gene model: {len(gene_model)} genes from "
                         f"{args.gene_model}\n")
    dictionary = load_dictionary(args.dictionary, args.panel)
    anchors = load_anchors(args.anchors, args.panel)
    cytobands = load_cytobands(args.cytoband_bed)
    records = parse_vcf(args.vcf)
    excl = load_excluded_junctions(args.excluded_junctions)
    if excl:
        sys.stderr.write(
            f"excluded junctions: {len(set(k[0] for k in excl))} chromosome "
            f"pair(s) from {args.excluded_junctions}\n")

    rows = annotate(records, panel, dictionary, anchors, args.sample,
                    cytobands, args.panel, gene_model)

    with open(args.output, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        dropped = []
        keep = []
        for r in rows:
            why = excluded_reason(excl, r)
            (dropped if why else keep).append((r, why))
        w.writeheader()
        w.writerows([r for r, _ in keep])

    if dropped:
        sys.stderr.write(f"dropped {len(dropped)} junction(s) on the exclusion list:\n")
        seen = set()
        for r, why in dropped:
            k = (r.get("chrom_a"), r.get("pos_a"), r.get("chrom_b"), r.get("pos_b"))
            if k in seen:
                continue
            seen.add(k)
            sys.stderr.write(
                f"  {r.get('gene_a','')} x {r.get('gene_b','')}  "
                f"{k[0]}:{k[1]} x {k[2]}:{k[3]}  — {why}\n")

    # Count the table that was written, not the one before exclusion. These
    # were computed over `rows` while `keep` is what reaches the file, so the
    # summary described a table that did not exist: dropping four junctions
    # left the reportable figure unchanged on all three validation samples.
    written = [r for r, _ in keep]
    n_report = sum(1 for r in written if r["reportable"] == "yes")
    n_named = sum(1 for r in written if r["known_pair"])
    sys.stderr.write(
        f"{args.sample} [{args.panel}]: {len(written)} on-panel SV records "
        f"({len(dropped)} excluded), "
        f"{n_report} reportable, {n_named} named by dictionary "
        f"({len(dictionary)} pairs, {len(anchors)} anchor tokens in scope) "
        f"-> {args.output}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
