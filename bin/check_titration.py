#!/usr/bin/env python3
"""
check_titration.py
==================

Track a fixed set of known findings across sequencing depths, to decide where
a run can stop.

The question this answers is not "how much depth do we get" but "at what point
does a finding we know is there stop being called". Coverage curves saturate
long before that: on this cohort every panel region cleared 15x by 16 hours,
yet a 7-read duplication call at 24 hours is a 5-read call at 16, and whether
that survives is a different question from whether the region was covered.

Findings are declared once, in TRUTH below, taken from the accredited assay.
Each is matched against the evidence file its own caller writes, so a fusion is
checked against the annotated SV table and a duplication against the focal-dup
table. Negatives are checked too: "no hotspot ITD" is a finding, and a run
length that starts producing one has failed in the other direction.

Reads directly from each timepoint's report.zip, so nothing needs re-running.

Usage:
  check_titration.py --results-glob 'results_*hr' --out titration_check.tsv
"""

from __future__ import annotations

import argparse
import csv
import io
import re
import sys
import zipfile
from pathlib import Path
from typing import Dict, List, Optional

__version__ = "0.1.0"


# -----------------------------------------------------------------------------
# What we expect to see. Sample keys are matched as a prefix, so BC13 matches
# BC13_16hr, BC13_24hr and so on.
# -----------------------------------------------------------------------------
TRUTH = [
    # sample, label,                    kind,        matcher
    ("BC13", "NPM1 frameshift",         "variant",   {"gene": "NPM1", "csq": "frameshift"}),
    ("BC13", "FLT3 p.Val592Ala",        "variant",   {"gene": "FLT3", "hgvsp": "Val592Ala"}),
    ("BC13", "no hotspot FLT3-ITD",     "absent_dup", {"gene": "FLT3"}),

    ("BC14", "UBTF-TD exon 13",         "dup",       {"gene": "UBTF"}),
    ("BC14", "FLT3-ITD, juxtamembrane", "dup",       {"gene": "FLT3"}),

    ("BC15", "DEK::NUP214 t(6;9)",      "fusion",    {"genes": ("DEK", "NUP214")}),
    ("BC15", "FLT3 p.Asn841Ile",        "variant",   {"gene": "FLT3", "hgvsp": "Asn841Ile"}),
    ("BC15", "no hotspot FLT3-ITD",     "absent_dup", {"gene": "FLT3"}),
]


def read_zip_tsv(zf: zipfile.ZipFile, pattern: str, sample: str) -> List[dict]:
    """Rows of the first member matching pattern for this sample."""
    rx = re.compile(pattern)
    for name in zf.namelist():
        base = name.rsplit("/", 1)[-1]
        if base.startswith(sample) and rx.search(base):
            with zf.open(name) as fh:
                text = io.TextIOWrapper(fh, encoding="utf-8", errors="replace")
                return [r for r in csv.DictReader(text, delimiter="\t")]
    return []


def num(v, default=0.0) -> float:
    try:
        return float(str(v).strip())
    except (TypeError, ValueError):
        return default


# -----------------------------------------------------------------------------
# Matchers. Each returns (found, evidence_string).
# -----------------------------------------------------------------------------
def match_variant(rows, spec):
    for r in rows:
        if (r.get("gene") or "").upper() != spec["gene"].upper():
            continue
        if "hgvsp" in spec and spec["hgvsp"].lower() not in (r.get("hgvsp") or "").lower():
            continue
        if "csq" in spec and spec["csq"] not in (r.get("consequence") or "").lower():
            continue
        return True, (f"VAF {r.get('tumor_af_pct','?')}%, "
                      f"{r.get('ALT_COUNT','?')}/{r.get('DP','?')} reads")
    return False, ""


def dup_gene(r) -> str:
    """Gene for a duplication row, from either schema.

    The caller writes gene and confidence columns. The report bundle rewrites
    the table into a consensus schema with neither: the graded call carries its
    label in status, an ungraded one carries the literal "negative", and the
    original fields are folded into raw_calls. Both are read here so the check
    works against a bundle or a raw results tree.
    """
    if r.get("gene"):
        return r["gene"].upper()
    st = (r.get("status") or "").upper()
    if st and st != "NEGATIVE":
        return st.split("-")[0]
    m = re.search(r"gene=([A-Za-z0-9]+)", r.get("raw_calls") or "")
    return m.group(1).upper() if m else ""


def dup_graded(r) -> bool:
    """Whether this row is a graded call rather than a positional artefact."""
    if r.get("confidence"):
        return r["confidence"].lower() in ("high", "moderate")
    st = (r.get("status") or "").strip().lower()
    if st in ("negative", "no_itd", "no-itd", ""):
        return False
    m = re.search(r"confidence=(\w+)", r.get("raw_calls") or "")
    return m.group(1).lower() in ("high", "moderate") if m else True


def dup_evidence(r) -> str:
    length = r.get("dup_length_bp") or r.get("length_bp") or "?"
    sup = r.get("support_reads") or ""
    span = r.get("spanning_reads") or ""
    if not sup:
        m = re.search(r"support_reads=(\d+)", r.get("raw_calls") or "")
        sup = m.group(1) if m else "?"
        m = re.search(r"spanning_reads=(\d+)", r.get("raw_calls") or "")
        span = m.group(1) if m else "?"
    dom = r.get("hotspot") or r.get("domain") or ""
    conf = r.get("confidence") or ""
    if not conf:
        m = re.search(r"confidence=(\w+)", r.get("raw_calls") or "")
        conf = m.group(1) if m else ""
    return (f"{length} bp, {sup}/{span} reads"
            + (f", {dom}" if dom else "") + (f", {conf}" if conf else ""))


def match_dup(rows, spec):
    """A graded duplication. Position matters more than support: an off-hotspot
    call is an artefact however many reads it carries, which is why the caller
    grades on position and why only high or moderate counts here."""
    best = None
    for r in rows:
        if dup_gene(r) != spec["gene"].upper() or not dup_graded(r):
            continue
        if best is None or num(r.get("support_reads") or r.get("vaf_pct_mean")) > \
                           num(best.get("support_reads") or best.get("vaf_pct_mean")):
            best = r
    return (True, dup_evidence(best)) if best is not None else (False, "")


def match_absent_dup(rows, spec):
    """A correct negative. Any graded call for this gene is a failure."""
    same = [r for r in rows if dup_gene(r) == spec["gene"].upper()]
    hits = [r for r in same if dup_graded(r)]
    if hits:
        return False, f"{len(hits)} graded call(s) — false positive"
    return True, (f"{len(same)} call(s), none graded" if same else "no calls")


def match_fusion(rows, spec):
    a, b = (g.upper() for g in spec["genes"])
    best = None
    for r in rows:
        ga = (r.get("gene_a") or "").upper()
        gb = (r.get("gene_b") or "").upper()
        if {a, b} <= {ga, gb} or (a in ga and b in gb) or (b in ga and a in gb):
            if best is None or num(r.get("support_reads")) > num(best.get("support_reads")):
                best = r
    if best is None:
        return False, ""
    named = (best.get("known_pair") or best.get("known_mm_pair") or "").strip()
    return True, (f"{best.get('support_reads','?')} reads, "
                  f"{best.get('n_callers','?')} caller(s)"
                  + (f", named {named}" if named else ", unnamed"))


MATCHERS = {"variant": match_variant, "dup": match_dup,
            "absent_dup": match_absent_dup, "fusion": match_fusion}

SOURCES = {"variant": r"clinical.*\.tsv$", "dup": r"flt3_consensus\.tsv$",
           "absent_dup": r"flt3_consensus\.tsv$", "fusion": r"annotated\.tsv$"}


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-glob", default="results_*hr",
                    help="Glob for the per-timepoint result directories "
                         "[results_*hr]. Each must contain report.zip.")
    ap.add_argument("--out", type=Path, default=None,
                    help="Write the matrix as TSV as well as printing it.")
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args = ap.parse_args()

    zips = {}
    for d in sorted(Path(".").glob(args.results_glob)):
        z = d / "report.zip"
        if not z.is_file():
            continue
        m = re.search(r"(\d+)hr", d.name)
        if m:
            zips[int(m.group(1))] = z
    if not zips:
        sys.stderr.write(f"ERROR: no report.zip under {args.results_glob}\n")
        return 1

    times = sorted(zips)
    sys.stderr.write(f"timepoints: {', '.join(str(t) + 'h' for t in times)}\n\n")

    # cache: (time, sample, kind) -> rows
    cache: Dict[tuple, List[dict]] = {}

    def rows_for(t, sample, kind):
        key = (t, sample, kind)
        if key not in cache:
            with zipfile.ZipFile(zips[t]) as zf:
                # sample keys are prefixes; find the full id in this bundle
                full = sample
                for n in zf.namelist():
                    b = n.rsplit("/", 1)[-1]
                    if b.startswith(sample + "_"):
                        full = b.split(".")[0].split("_report")[0]
                        full = "_".join(full.split("_")[:2])
                        break
                cache[key] = read_zip_tsv(zf, SOURCES[kind], full)
        return cache[key]

    results = []
    w = max(len(lbl) for _, lbl, _, _ in TRUTH) + 2
    header = f"{'':<6}{'FINDING':<{w}}" + "".join(f"{str(t)+'h':<10}" for t in times)
    print(header)
    print("-" * len(header))

    for sample, label, kind, spec in TRUTH:
        line = f"{sample:<6}{label:<{w}}"
        row = {"sample": sample, "finding": label, "kind": kind}
        for t in times:
            found, ev = MATCHERS[kind](rows_for(t, sample, kind), spec)
            line += f"{('yes' if found else 'NO'):<10}"
            row[f"{t}h"] = "yes" if found else "no"
            row[f"{t}h_evidence"] = ev
        print(line)
        results.append(row)

    # first depth at which every finding holds
    print()
    for t in times:
        n = sum(1 for r in results if r[f"{t}h"] == "yes")
        mark = "  <- all findings recovered" if n == len(results) else ""
        print(f"  {t:>3}h  {n}/{len(results)}{mark}")

    print("\nEvidence per finding:")
    for r in results:
        print(f"\n  {r['sample']}  {r['finding']}")
        for t in times:
            ev = r.get(f"{t}h_evidence") or "-"
            print(f"    {str(t)+'h':<5} {r[f'{t}h']:<4} {ev}")

    if args.out:
        cols = ["sample", "finding", "kind"]
        for t in times:
            cols += [f"{t}h", f"{t}h_evidence"]
        with open(args.out, "w", newline="") as fh:
            wr = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                                lineterminator="\n")
            wr.writeheader()
            wr.writerows(results)
        sys.stderr.write(f"\nwrote {args.out}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
