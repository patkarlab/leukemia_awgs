#!/usr/bin/env python3
"""
augment_sv_support.py

The SURVIVOR-merged VCF drops per-caller read support (it keeps only SUPP/
SUPP_VEC and coordinates), so support_reads cannot be recovered from the merged
file or from mm_annotated.tsv. This script layers the real per-caller support
back on, read directly from each caller's own VCF:

  Sniffles : INFO/SUPPORT
  CuteSV   : INFO/RE
  Severus  : FORMAT/DV   (fallback INFO/SUPP_READS first colon-field)

For every row in an mm_annotated.tsv it matches the junction (orientation-
agnostic, within --tol bp on BOTH breakpoints to absorb the few-bp shifts
SURVIVOR introduces) against each caller's calls and writes:

  support_sniffles, support_cutesv, support_severus   per-caller read counts
  support_reads (overwritten)                          = max of those present

Matching is by coordinate only; no specific variants or findings are hardcoded.
Same-chromosome SV types (DEL/DUP/INV/INS) are matched on the primary breakpoint
plus END; translocations (TRA/BND) on both breakends.

Standard-library only.

Usage:
  python3 augment_sv_support.py \\
      --annotated  <sample>.mm_annotated.tsv \\
      --sniffles   <sample>.sniffles.t2t.vcf.gz \\
      --cutesv     <sample>.cutesv.t2t.vcf.gz \\
      --severus    <sample>.severus.vcf \\
      --output     <sample>.mm_annotated.tsv      (in place is fine)
      [--tol 25]
"""

import argparse
import csv
import gzip
import os
import sys
from typing import Dict, List, Optional, Tuple


def open_any(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def chrom_key(chrom):
    name = chrom[3:] if chrom.lower().startswith("chr") else chrom
    if name.isdigit():
        return (0, int(name), "")
    return (1, {"X": 0, "Y": 1, "M": 2, "MT": 2}.get(name.upper(), 99), name.upper())


def parse_info(info_field):
    d = {}
    for part in info_field.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            d[k] = v
        elif part:
            d[part] = "True"
    return d


def parse_bnd_alt(alt):
    for opener, closer in (("]", "]"), ("[", "[")):
        if opener in alt:
            try:
                inner = alt.split(opener)[1].split(closer)[0]
                c, p = inner.rsplit(":", 1)
                return c, int(p)
            except (IndexError, ValueError):
                continue
    return None, None


def canon(ca, pa, cb, pb):
    """Canonical orientation-agnostic ordering of the two ends."""
    ka, kb = (chrom_key(ca), pa), (chrom_key(cb), pb)
    return (ca, pa, cb, pb) if ka <= kb else (cb, pb, ca, pa)


def caller_support(value_kind, info, fmt, sample):
    """Extract this caller's own support count from one of its records."""
    if value_kind == "RE":
        v = info.get("RE", "")
        return int(v) if str(v).isdigit() else None
    if value_kind == "SUPPORT":
        v = info.get("SUPPORT", "")
        return int(v) if str(v).isdigit() else None
    if value_kind == "DV":
        if fmt and sample:
            f = dict(zip(fmt.split(":"), sample.split(":")))
            dv = f.get("DV", "")
            if dv.isdigit():
                return int(dv)
        sr = info.get("SUPP_READS", "")
        first = sr.split(":")[0] if sr else ""
        return int(first) if first.isdigit() else None
    return None


def load_caller(path, value_kind):
    """Return a list of (canon_a_chrom, a_pos, b_chrom, b_pos, support) for
    every SV record in a caller VCF."""
    out = []
    if not path or not os.path.isfile(path):
        return out
    with open_any(path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            chrom, pos, _id, _ref, alt, _q, _filt, info_field = cols[:8]
            pos = int(pos)
            info = parse_info(info_field)
            svtype = info.get("SVTYPE", "")
            fmt = cols[8] if len(cols) > 8 else ""
            sample = cols[9] if len(cols) > 9 else ""
            sup = caller_support(value_kind, info, fmt, sample)
            if sup is None:
                continue
            if svtype in ("BND", "TRA"):
                mc, mp = parse_bnd_alt(alt)
                if mc is None:
                    mc, end = info.get("CHR2"), info.get("END")
                    mp = int(end) if (end and end.isdigit()) else None
                if mc is None or mp is None:
                    continue
                ca, pa, cb, pb = canon(chrom, pos, mc, mp)
            else:
                end = info.get("END")
                if not (end and end.isdigit()):
                    continue
                ca, pa, cb, pb = canon(chrom, pos, chrom, int(end))
            out.append((ca, pa, cb, pb, sup))
    return out


def best_match(ca, pa, cb, pb, records, tol):
    """Support of the NEAREST caller record whose BOTH ends fall within tol bp
    (sum of per-end distances minimised). Nearest, not max, so that in dense
    loci where a caller emits several breakends a few bp apart, the count is
    attributed to the corresponding breakpoint rather than the strongest nearby
    one. Ties in distance are broken by higher support."""
    best_sup = None
    best_dist = None
    for (rca, rpa, rcb, rpb, sup) in records:
        if rca == ca and rcb == cb and abs(rpa - pa) <= tol and abs(rpb - pb) <= tol:
            dist = abs(rpa - pa) + abs(rpb - pb)
            if (best_dist is None or dist < best_dist
                    or (dist == best_dist and (best_sup is None or sup > best_sup))):
                best_dist, best_sup = dist, sup
    return best_sup


def main():
    ap = argparse.ArgumentParser(description="Layer per-caller read support onto mm_annotated.tsv.")
    ap.add_argument("--annotated", required=True)
    ap.add_argument("--sniffles", default=None)
    ap.add_argument("--cutesv", default=None)
    ap.add_argument("--severus", default=None)
    ap.add_argument("--output", required=True)
    ap.add_argument("--tol", type=int, default=25, help="bp tolerance per breakpoint (default 25).")
    args = ap.parse_args()

    snf = load_caller(args.sniffles, "SUPPORT")
    cut = load_caller(args.cutesv, "RE")
    sev = load_caller(args.severus, "DV")
    sys.stderr.write(f"loaded support records: sniffles={len(snf)} cutesv={len(cut)} severus={len(sev)}\n")

    with open(args.annotated, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        in_cols = list(reader.fieldnames or [])
        rows = [dict(r) for r in reader]

    new_cols = ["support_sniffles", "support_cutesv", "support_severus"]
    out_cols = in_cols + [c for c in new_cols if c not in in_cols]

    n_pop = 0
    for r in rows:
        ca_, pa_ = r.get("chrom_a", ""), r.get("pos_a", "")
        cb_, pb_ = r.get("chrom_b", ""), r.get("pos_b", "")
        if not (ca_ and pa_ and cb_ and pb_ and pa_.lstrip("-").isdigit() and pb_.lstrip("-").isdigit()):
            for c in new_cols:
                r[c] = ""
            continue
        ca, pa, cb, pb = canon(ca_, int(pa_), cb_, int(pb_))
        # Gate on the callers the merge already established for this junction
        # (from SUPP_VEC). Coordinate-matching only fetches that caller's read
        # count; it never credits a caller absent from this row, which prevents
        # a different nearby breakend in a dense locus from bleeding across.
        listed = {x.strip().lower() for x in (r.get("callers") or "").split(",") if x.strip()}
        def gated(name, recs):
            if listed and name not in listed:
                return None
            return best_match(ca, pa, cb, pb, recs, args.tol)
        s = gated("sniffles", snf)
        c = gated("cutesv", cut)
        v = gated("severus", sev)
        r["support_sniffles"] = str(s) if s is not None else ""
        r["support_cutesv"]   = str(c) if c is not None else ""
        r["support_severus"]  = str(v) if v is not None else ""
        present = [x for x in (s, c, v) if x is not None]
        if present:
            r["support_reads"] = str(max(present))
            n_pop += 1

    with open(args.output, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=out_cols, delimiter="\t",
                           extrasaction="ignore", restval="", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    sys.stderr.write(f"populated support_reads on {n_pop}/{len(rows)} rows -> {args.output}\n")


if __name__ == "__main__":
    main()
