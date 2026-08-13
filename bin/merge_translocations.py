#!/usr/bin/env python3
"""
merge_translocations.py

Collapse duplicate translocation (TRA) calls in an mm_annotated.tsv where the
same junction is reported multiple times - as reciprocal mates (Severus _1/_2
with A/B swapped) or as near-identical calls from different callers a few bp
apart. Two TRA rows are united when, after orientation normalisation, they share
the same chromosome pair AND BOTH breakpoints lie within --max-dist bp
(default 100). Non-TRA rows pass through unchanged.

Ig-aware second pass (added): immunoglobulin-partner translocations are special.
The Ig locus (IGH/IGK/IGL) is large (IGH ~1.5 Mb) and the Ig-side breakpoint of
the *same* event lands at different J/switch positions across callers, so two
genuine FGFR3::IGH calls have Ig-side coordinates far more than --max-dist apart
and never collapse under the positional rule - producing many redundant rows.
The second pass therefore collapses any partner::Ig translocations that share
the same unordered gene pair and a partner-side breakpoint within
--ig-partner-tol bp (default 1000), treating the whole Ig locus as one anchor
and ignoring the Ig-side coordinate. Non-Ig::non-Ig pairs keep the positional
rule unchanged. Ig loci are detected by gene-symbol prefix, so this is
reference-agnostic.

Standard-library only; runs in any python3 (incl. awgs_sv).

Merge semantics per cluster:
  - orientation normalised: each row's two ends are ordered canonically
    (chrom-sorted), so a mate's swapped A/B collapses onto its partner;
  - representative row = most callers, then most support_reads, then first;
  - positions/genes taken from the representative (gene stays with its end;
    a panel-gene annotation is preferred over OFF_PANEL);
  - callers = union across members; n_callers recomputed from that union;
  - known_mm_pair / known_freq = first non-empty across members;
  - filter = PASS if any member is PASS;
  - two columns added: n_merged, merged_sv_ids.

Single-linkage clustering (transitive) is used; clusters stay small so chaining
is not a concern.

Usage:
  python3 merge_translocations.py --input <sample>.mm_annotated.tsv [...] \
      [--outdir DIR] [--max-dist 100] [--ig-partner-tol 1000] [--keep-all-sv]
"""

# [cytoband-partner-annotation applied]
import argparse
import csv
import os
import sys

CALLER_ORDER = {"Sniffles": 0, "CuteSV": 1, "Severus": 2, "nanomonsv": 3}


def chrom_sort_key(chrom):
    name = chrom[3:] if chrom.lower().startswith("chr") else chrom
    if name.isdigit():
        return (0, int(name), "")
    return (1, {"X": 0, "Y": 1, "M": 2, "MT": 2}.get(name.upper(), 99), name.upper())


def to_int(value, default=None):
    try:
        return int(str(value).strip())
    except (ValueError, TypeError):
        return default


def is_dispersed_anchor(gene):
    """True if the gene is a dispersed-breakpoint translocation anchor: an
    immunoglobulin locus (IGH/IGK/IGL, incl. @-suffixed or _locus forms) or MYC.
    These loci carry the same event's breakpoint at scattered positions across a
    wide region, so redundant calls must be collapsed on the *partner* side
    rather than by requiring both breakpoints to coincide. Reference-agnostic;
    keyed on gene symbol so it is independent of coordinates/assembly."""
    g = (gene or "").strip().upper().rstrip("@")
    if g == "MYC":
        return True
    return g.startswith("IGH") or g.startswith("IGK") or g.startswith("IGL")


def canonical_ends(row):
    """Return ((chromA,posA,geneA),(chromB,posB,geneB)) chrom-sorted so that
    reciprocal mates normalise identically."""
    a = (row["chrom_a"], to_int(row["pos_a"], 0), row.get("gene_a", ""),
         (row.get("gene_a_source", "") or "coordinate"))
    b = (row["chrom_b"], to_int(row["pos_b"], 0), row.get("gene_b", ""),
         (row.get("gene_b_source", "") or "coordinate"))
    ka = (chrom_sort_key(a[0]), a[1])
    kb = (chrom_sort_key(b[0]), b[1])
    return (a, b) if ka <= kb else (b, a)


def callers_set(row):
    raw = (row.get("callers") or "").strip()
    return {c.strip() for c in raw.split(",") if c.strip()}


_SOURCE_RANK = {"panel": 0, "cytoband": 1, "coordinate": 2}


def best_end(members, end_index):
    """Pick the best (gene, source) for end 0 or 1 across clustered members.

    Prefers the strongest provenance (panel > cytoband > coordinate). This
    replaces the old string test against "OFF_PANEL": off-panel partners are
    now labelled by cytoband, so the merge must rank on the provenance token
    rather than a magic gene name. Ties are broken by taking the end from the
    representative member (most callers / most support), so an adjacent-band
    disagreement resolves to the highest-support call rather than input order.
    """
    best = None
    best_rank = None
    ranked = sorted(members, key=_rep_key, reverse=True)
    for m in ranked:
        end = m["_ends"][end_index]
        gene = end[2].strip()
        source = (end[3] if len(end) > 3 else "coordinate").strip() or "coordinate"
        if not gene:
            continue
        rank = _SOURCE_RANK.get(source, 3)
        if best is None or rank < best_rank:
            best, best_rank = (gene, source), rank
    return best if best is not None else ("", "coordinate")


def best_gene(members, end_index):
    """Backwards-compatible shim: return only the gene label."""
    return best_end(members, end_index)[0]


def first_nonempty(members, col):
    for m in members:
        v = (m.get(col) or "").strip()
        if v:
            return v
    return ""


def _rep_key(m):
    return (to_int(m.get("n_callers"), 0), to_int(m.get("support_reads"), 0))


def merge_cluster(members):
    # representative: most callers, then most support_reads, then first.
    rep = max(members, key=_rep_key)

    end1, end2 = rep["_ends"]
    callers = set()
    for m in members:
        callers |= callers_set(m)
    callers_sorted = sorted(callers, key=lambda c: (CALLER_ORDER.get(c, 99), c))

    # max, not sum: the same junction seen by multiple callers is overlapping
    # read sets, so summing would double/triple-count the same molecules.
    supports = [to_int(m.get("support_reads")) for m in members]
    supports = [s for s in supports if s is not None]
    support_reads = str(max(supports)) if supports else ""

    def _max_col(col):
        vals = [to_int(m.get(col)) for m in members]
        vals = [v for v in vals if v is not None]
        return str(max(vals)) if vals else ""
    sup_snf = _max_col("support_sniffles")
    sup_cut = _max_col("support_cutesv")
    sup_sev = _max_col("support_severus")

    filt = "PASS" if any((m.get("filter") or "").strip() == "PASS" for m in members) \
        else (rep.get("filter") or "")

    return {
        "sample": rep.get("sample", ""),
        "sv_id": rep.get("sv_id", ""),
        "sv_type": "TRA",
        "filter": filt,
        "chrom_a": end1[0], "pos_a": str(end1[1]),
        "gene_a": best_end(members, 0)[0],
        "gene_a_source": best_end(members, 0)[1],
        "chrom_b": end2[0], "pos_b": str(end2[1]),
        "gene_b": best_end(members, 1)[0],
        "gene_b_source": best_end(members, 1)[1],
        "known_mm_pair": first_nonempty(members, "known_mm_pair"),
        "known_freq": first_nonempty(members, "known_freq"),
        "callers": ",".join(callers_sorted),
        "n_callers": str(len(callers_sorted)),
        "supp_vec": rep.get("supp_vec", ""),
        "support_reads": support_reads,
        "support_sniffles": sup_snf,
        "support_cutesv": sup_cut,
        "support_severus": sup_sev,
        "n_merged": str(len(members)),
        "merged_sv_ids": ",".join(m.get("sv_id", "") for m in members),
    }


def cluster_tra(tra_rows, max_dist):
    """Single-linkage union-find on canonical-pair + both-ends-within-threshold."""
    n = len(tra_rows)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        parent[find(x)] = find(y)

    for i in range(n):
        e1i, e2i = tra_rows[i]["_ends"]
        for j in range(i + 1, n):
            e1j, e2j = tra_rows[j]["_ends"]
            if e1i[0] != e1j[0] or e2i[0] != e2j[0]:
                continue
            if abs(e1i[1] - e1j[1]) <= max_dist and abs(e2i[1] - e2j[1]) <= max_dist:
                union(i, j)

    clusters = {}
    for i in range(n):
        clusters.setdefault(find(i), []).append(tra_rows[i])
    return clusters


def ig_aware_union(clusters, ig_partner_tol):
    """Second pass over positional clusters: collapse partner::anchor translocations
    whose anchor-side breakpoints are scattered across a dispersed anchor (Ig
    locus or MYC).

    The anchor side is identified by is_dispersed_anchor (IGH/IGK/IGL or MYC).
    Eligible clusters are those with exactly one anchor end. They are grouped by
    (Ig-side chromosome, partner chromosome) and then single-linkage clustered on
    the *partner* breakpoint position: two clusters merge if their partner
    positions are within --ig-partner-tol bp. Partner gene name is NOT used, so
    off-panel/unnamed partners (e.g. OFF_PANEL::MYC, IGH_locus::OFF_PANEL) still
    collapse on position. The Ig-side coordinate is ignored entirely. Clusters
    that are not partner::Ig (both Ig, or neither) are left exactly as the
    positional pass produced them.

    Returns a new clusters-style dict (values are member lists)."""
    cluster_list = list(clusters.values())
    n = len(cluster_list)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        parent[find(x)] = find(y)

    tol = max(1, ig_partner_tol)

    # Collect eligible clusters with their partner anchor (ig_chrom, partner_chrom, partner_pos).
    groups = {}
    for ci, members in enumerate(cluster_list):
        rep = max(members, key=_rep_key)
        (ca, pa, ga, _), (cb, pb, gb, _) = rep["_ends"]
        a_ig, b_ig = is_dispersed_anchor(ga), is_dispersed_anchor(gb)
        if a_ig == b_ig:                       # both Ig or neither -> leave alone
            continue
        if a_ig:
            ig_c, partner_c, partner_p = ca, cb, pb
        else:
            ig_c, partner_c, partner_p = cb, ca, pa
        groups.setdefault((ig_c, partner_c), []).append((partner_p, ci))

    # Single-linkage on partner position within each (ig_chrom, partner_chrom) group.
    for members in groups.values():
        members.sort()
        for k in range(1, len(members)):
            if members[k][0] - members[k - 1][0] <= tol:
                union(members[k][1], members[k - 1][1])

    merged = {}
    for ci in range(n):
        merged.setdefault(find(ci), []).extend(cluster_list[ci])
    return merged


def process_file(path, outdir, max_dist, keep_all_sv, ig_partner_tol):
    stem = os.path.basename(path)
    # default: translocations-only file; --keep-all-sv: full table w/ TRA merged
    suffix = ".sv_tra_merged.tsv" if keep_all_sv else ".translocations.tsv"
    base = stem
    for s_ in (".mm_annotated.tsv", "_mm_annotated.tsv", ".tsv"):
        if base.endswith(s_):
            base = base[: -len(s_)]
            break
    out_name = base + suffix
    os.makedirs(outdir, exist_ok=True)
    out_path = os.path.join(outdir, out_name)

    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        in_cols = list(reader.fieldnames or [])
        rows = [dict(r) for r in reader]

    out_cols = in_cols + [c for c in ("n_merged", "merged_sv_ids") if c not in in_cols]

    # Index TRA rows; precompute canonical ends.
    tra_idx = []
    for i, r in enumerate(rows):
        if (r.get("sv_type") or "").strip() == "TRA":
            r["_ends"] = canonical_ends(r)
            tra_idx.append(i)

    tra_rows = [rows[i] for i in tra_idx]
    clusters = cluster_tra(tra_rows, max_dist)
    clusters = ig_aware_union(clusters, ig_partner_tol)

    # Map each TRA row's original position -> its cluster's representative output,
    # emitted once at the earliest member position.
    pos_of = {id(r): tra_idx[k] for k, r in enumerate(tra_rows)}
    merged_by_first = {}
    member_first = {}
    for members in clusters.values():
        first_pos = min(pos_of[id(m)] for m in members)
        merged_by_first[first_pos] = merge_cluster(members)
        for m in members:
            member_first[pos_of[id(m)]] = first_pos

    out_rows = []
    n_tra_in = len(tra_rows)
    for i, r in enumerate(rows):
        if (r.get("sv_type") or "").strip() == "TRA":
            if member_first.get(i) == i:        # representative emission point
                out_rows.append(merged_by_first[i])
            # else: a non-first cluster member -> skip
        elif keep_all_sv:
            rr = {k: r.get(k, "") for k in in_cols}
            rr["n_merged"] = "1"
            rr["merged_sv_ids"] = r.get("sv_id", "")
            out_rows.append(rr)
        # else: non-TRA rows are dropped (translocations-only output)

    with open(out_path, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=out_cols, delimiter="\t",
                           extrasaction="ignore", restval="", lineterminator="\n")
        w.writeheader()
        for rr in out_rows:
            w.writerow(rr)

    n_tra_out = sum(1 for rr in out_rows if rr.get("sv_type") == "TRA")
    scope = "full SV table" if keep_all_sv else "translocations only"
    sys.stderr.write(
        f"[{stem}] TRA {n_tra_in} -> {n_tra_out} merged; "
        f"{len(out_rows)} rows written ({scope}); {out_path}\n"
    )
    return out_path


def main():
    ap = argparse.ArgumentParser(description="Unite near-identical translocation calls in mm_annotated.tsv.")
    ap.add_argument("-i", "--input", required=True, nargs="+", help="One or more *.mm_annotated.tsv files.")
    ap.add_argument("-o", "--outdir", default=None, help="Output dir (default: alongside each input).")
    ap.add_argument("--max-dist", type=int, default=100, help="Max bp between breakpoints to unite non-Ig pairs (default 100).")
    ap.add_argument("--ig-partner-tol", type=int, default=2000,
                    help="bp tolerance on the partner side when collapsing partner::Ig "
                         "translocations across the Ig locus (default 1000).")
    ap.add_argument("--keep-all-sv", action="store_true", help="Keep all SV types in output (TRA merged, others passthrough). Default: translocations only.")
    args = ap.parse_args()
    for path in args.input:
        if not os.path.isfile(path):
            sys.stderr.write(f"ERROR: not found: {path}\n"); sys.exit(1)
        outdir = args.outdir or os.path.dirname(os.path.abspath(path))
        process_file(path, outdir, args.max_dist, args.keep_all_sv, args.ig_partner_tol)


if __name__ == "__main__":
    main()
