#!/usr/bin/env python3
"""
Filter Clair3/VEP somatic candidates to the adaptive-WGS acute-leukaemia panel.

Panel membership and per-region labels are derived from the hg38 panel BED at
runtime (BED col-4 labels), NOT from a hardcoded gene dictionary. This is the
design decision carried over from the acute leukaemia pipeline, and it is what lets one
script serve both the AML and the ALL panel: the panel BED handed in at runtime
determines membership, so the AML run keeps AML genes and the ALL run keeps ALL
genes with no code path in common.

For each input *_somatic_candidates.tsv (VEP-annotated Clair3/ClairS-TO output):
  1. PANEL-GENE membership : keep variants whose gene symbol maps to a BED
     col-4 label (directly, via a small alias map, or via an IG-locus prefix
     when --include-ig is set).
  2. Consequence filter     : keep only protein-altering primary consequences.
  3. Population-AF flag      : pop_af_max above --max-pop-af is soft-flagged
                              COMMON_POLYMORPHISM (kept, not dropped).

Outputs (filenames unchanged for downstream compatibility):
  <stem>.al_filtered.tsv    in-panel + reportable rows, with Filter column
  <stem>.al_clinical.tsv    the Filter == PASS subset
  al_filter_summary.tsv     one row per input file

No sample-specific variant, gene-pair, or finding is hardcoded; only column
names and the panel BED are read.

Usage:
  filter_al_somatic_candidates.py \\
      --input  <sample>_somatic_candidates.tsv [...] \\
      --panel-bed AML_panel_hg38.bed \\
      --outdir al_filtered
"""
import argparse
import csv
import os
import sys

MISSING = "-1"

# ----------------------------------------------------------------------------
# Panel membership and labels are derived from the hg38 panel BED at runtime.
# Nothing about the gene set is hardcoded here, so the filter does not drift
# when the panel is revised. The BED is the single source of truth: a variant
# is in-panel iff its gene symbol maps to a BED col-4 label.
#
# The only hand-maintained data is ALIAS_TO_SYMBOL: cases where VEP emits a
# gene symbol that is not literally present in any BED label. These are true
# synonyms, not membership decisions. Extend only for genuine naming
# mismatches, never to add or remove panel genes.
# ----------------------------------------------------------------------------
# Gene synonyms and excluded loci are loaded from TSV assets rather than
# written here. Nothing about the gene set, and no clinical prior, is
# hardcoded: panel membership comes from the BED, synonyms from
# assets/gene_aliases.tsv, and excluded loci from assets/excluded_loci.tsv.
ALIAS_TO_SYMBOL = {}
IG_PREFIXES = ()


def load_aliases(path):
    """alias -> symbol. Naming mismatches only, never membership decisions."""
    out = {}
    if not path:
        return out
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            a = (row.get("alias") or "").strip().upper()
            s = (row.get("symbol") or "").strip().upper()
            if a and s:
                out[a] = s
    return out


def load_excluded_loci(path):
    """Label prefixes excluded from the SNV report.

    The Ig and TR loci are on the panel as fusion anchors, but undergo somatic
    hypermutation and V(D)J recombination, so variant calls across them are
    recombination artefacts. They stay fully active for fusion detection on the
    T2T track; this exclusion applies to the SNV report only, and --include-ig
    disables it.
    """
    if not path:
        return ()
    out = []
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            pfx = (row.get("prefix") or "").strip().upper()
            if pfx:
                out.append(pfx)
    return tuple(out)


REPORTABLE_CONSEQUENCES = {
    "frameshift_variant",
    "inframe_deletion",
    "inframe_insertion",
    "missense_variant",
    "start_lost",
    "stop_gained",
    "stop_lost",
    "stop_retained_variant",
}

# Final column order; any unexpected columns are appended after these.
PREFERRED_COLS = [
    "chrom", "pos", "ref", "alt", "qual", "variant_type",
    "gene", "panel_label", "transcript", "hgvsc", "hgvsp", "biotype", "canonical",
    "consequence", "csq_primary", "impact", "exon_rank", "exon_total",
    "rs_id", "pop_af_max", "pop_af_max_source", "clinvar_sig",
    "tumor_af", "tumor_af_pct", "REF_COUNT", "ALT_COUNT", "DP", "Filter",
]

MISSING_TOKENS = {"", "-1", ".", "nan", "NA"}


def is_missing(value):
    return value is None or str(value).strip() in MISSING_TOKENS


def safe_float(value):
    """Return float(value) or None if missing/non-numeric."""
    if is_missing(value):
        return None
    try:
        return float(value)
    except (ValueError, TypeError):
        return None


def clean_cell(value):
    if value is None:
        return MISSING
    s = str(value).strip()
    return MISSING if s in ("", ".") else s


def _label_from_bed_name(raw_name):
    """Convert a BED col-4 label to the downstream panel_label style:
    compound separators '/' and '+' become '_'. E.g. 'FGFR3/NSD2' ->
    'FGFR3_NSD2', 'TP53+TNFSF12' -> 'TP53_TNFSF12'. This preserves the
    label convention the dashboard already consumes."""
    return raw_name.strip().replace("/", "_").replace("+", "_")


def load_panel(bed_path):
    """Build {gene_symbol_upper: panel_label} from the hg38 panel BED.

    Each BED col-4 label may be compound ('FCRL5/FCRL4', 'TP53+TNFSF12').
    Every component symbol is indexed to the same normalized panel_label, so a
    VEP call on any single component (e.g. TNFSF12, or FCRL4) resolves to the
    region's label.
    """
    symbol_to_label = {}
    with open(bed_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if len(fields) < 4:
                continue
            raw_name = fields[3].strip()
            if not raw_name:
                continue
            label = _label_from_bed_name(raw_name)
            # Split compound region names into their component symbols.
            for component in raw_name.replace("+", "/").split("/"):
                component = component.strip().upper()
                if component:
                    symbol_to_label[component] = label
    if not symbol_to_label:
        sys.stderr.write(
            f"ERROR: no panel regions parsed from BED: {bed_path}\n")
        sys.exit(1)
    return symbol_to_label


def resolve_panel_label(gene_symbol, panel, include_ig):
    """Return the panel_label for a gene symbol, or None if off-panel.

    Resolution order: direct BED match -> alias -> BED match -> IG-locus
    prefix (only when include_ig). Membership is entirely BED-driven; the
    alias map only rewrites a symbol into one the BED already knows.
    """
    g = str(gene_symbol).strip().upper()
    if not g or g in MISSING_TOKENS:
        return None
    if g in panel:
        return panel[g]
    alias = ALIAS_TO_SYMBOL.get(g)
    if alias and alias.upper() in panel:
        return panel[alias.upper()]
    if include_ig:
        for prefix in IG_PREFIXES:
            if g.startswith(prefix):
                # Prefer an actual BED label for this locus if one exists.
                for sym, lab in panel.items():
                    if sym.startswith(prefix):
                        return lab
                return f"{prefix}_locus"
    return None


def split_exon(value):
    s = str(value).strip()
    if s in ("", "-1", "."):
        return MISSING, MISSING
    if "/" in s:
        rank, total = s.split("/", 1)
        return (rank.strip() or MISSING), (total.strip() or MISSING)
    return s, MISSING


def af_to_pct(value):
    f = safe_float(value)
    if f is None:
        return MISSING
    return f"{f * 100:.2f}"


def filter_one_file(input_path, outdir, max_pop_af, include_ig, panel,
                    keep_off_panel):
    stem = os.path.basename(input_path)
    for suffix in (".tsv", ".txt"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]

    with open(input_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        in_cols = list(reader.fieldnames or [])
        rows = [dict(r) for r in reader]

    n_input = len(rows)
    kept = []
    n_off_panel = 0
    for row in rows:
        # 1. Normalise missing values.
        for key in list(row.keys()):
            row[key] = clean_cell(row[key])

        # 2. Ensure count columns exist.
        for key in ("REF_COUNT", "ALT_COUNT", "DP"):
            row.setdefault(key, MISSING)
            if is_missing(row[key]):
                row[key] = MISSING

        # 2a. Exon split (and drop original 'exon').
        if "exon" in row:
            row["exon_rank"], row["exon_total"] = split_exon(row.pop("exon"))

        # 2b. AF percentage.
        if "tumor_af" in row:
            row["tumor_af_pct"] = af_to_pct(row["tumor_af"])

        # 2c. Drop domains.
        row.pop("domains", None)

        # 3. Panel membership + primary consequence.
        label = resolve_panel_label(row.get("gene", ""), panel, include_ig)
        primary = str(row.get("consequence", "")).strip().split("&")[0]

        # 4. Filters. Consequence is always enforced. Off-panel rows are
        #    dropped unless --keep-off-panel is set (audit mode), in which
        #    case they are retained with panel_label = OFF_PANEL.
        if primary not in REPORTABLE_CONSEQUENCES:
            continue
        if label is None:
            n_off_panel += 1
            if not keep_off_panel:
                continue
            label = "OFF_PANEL"
        row["panel_label"] = label
        row["csq_primary"] = primary

        # 5. Soft polymorphism flag.
        pop_af = safe_float(row.get("pop_af_max"))
        row["Filter"] = "COMMON_POLYMORPHISM" if (pop_af is not None and pop_af > max_pop_af) else "PASS"

        kept.append(row)

    # Sort by panel gene, chrom, numeric position.
    def sort_key(r):
        return (r.get("panel_label", ""), r.get("chrom", ""), safe_float(r.get("pos")) or 0.0)
    kept.sort(key=sort_key)

    # Column order.
    all_cols = set()
    for r in kept:
        all_cols.update(r.keys())
    ordered = [c for c in PREFERRED_COLS if c in all_cols]
    ordered += [c for c in all_cols if c not in ordered]

    os.makedirs(outdir, exist_ok=True)
    filtered_path = os.path.join(outdir, f"{stem}.al_filtered.tsv")
    clinical_path = os.path.join(outdir, f"{stem}.al_clinical.tsv")

    def write_tsv(path, records):
        with open(path, "w", newline="") as out:
            w = csv.DictWriter(out, fieldnames=ordered, delimiter="\t", lineterminator="\n",
                               extrasaction="ignore", restval=MISSING)
            w.writeheader()
            for rec in records:
                w.writerow(rec)

    write_tsv(filtered_path, kept)
    clinical = [r for r in kept if r["Filter"] == "PASS"]
    write_tsv(clinical_path, clinical)

    n_pass = len(clinical)
    n_common = sum(1 for r in kept if r["Filter"] == "COMMON_POLYMORPHISM")
    genes_hit = sorted({r["panel_label"] for r in clinical})
    sys.stderr.write(
        f"[{stem}] input={n_input}  in-panel+reportable={len(kept)}  "
        f"PASS={n_pass}  COMMON_POLYMORPHISM={n_common}  off-panel-reportable={n_off_panel}\n"
    )
    sys.stderr.write(f"[{stem}] PASS genes: {', '.join(genes_hit) if genes_hit else '(none)'}\n")
    return {"sample": stem, "input": n_input, "kept": len(kept),
            "pass": n_pass, "common_polymorphism": n_common,
            "off_panel_reportable": n_off_panel}


def parse_args():
    ap = argparse.ArgumentParser(
        description="Filter Clair3/VEP somatic candidates to the acute-leukaemia panel "
                    "(membership derived from the hg38 panel BED).")
    ap.add_argument("-i", "--input", required=True, nargs="+",
                    help="One or more *_somatic_candidates.tsv files.")
    ap.add_argument("--panel-bed", required=True,
                    help="hg38 panel BED (col-4 labels define the in-panel gene "
                         "set). Pass ${params.panel_bed_hg38}.")
    ap.add_argument("-o", "--outdir", default="al_filtered",
                    help="Output directory (default: ./al_filtered).")
    ap.add_argument("--max-pop-af", type=float, default=0.01,
                    help="pop_af_max above this is flagged COMMON_POLYMORPHISM "
                         "(default: 0.01).")
    ap.add_argument("--aliases", default=None,
                    help="TSV of gene synonyms (alias, symbol). Optional.")
    ap.add_argument("--excluded-loci", default=None,
                    help="TSV of label prefixes excluded from the SNV report "
                         "(prefix, reason). Optional.")
    ap.add_argument("--include-ig", action="store_true",
                    help="Include IGH/IGK/IGL-locus genes (off by default: "
                         "IG-locus SNVs are dominated by somatic hypermutation).")
    ap.add_argument("--keep-off-panel", action="store_true",
                    help="Retain reportable variants that are off-panel, labeled "
                         "OFF_PANEL, instead of dropping them. Audit aid for panel "
                         "transitions; off by default.")
    return ap.parse_args()


def main():
    args = parse_args()
    if not os.path.isfile(args.panel_bed):
        sys.stderr.write(f"ERROR: panel BED not found: {args.panel_bed}\n")
        sys.exit(1)
    panel = load_panel(args.panel_bed)

    # Module-level tables are populated from the asset TSVs. Both are optional:
    # without an alias file, only literal symbol matches resolve; without an
    # excluded-loci file, nothing is excluded on locus grounds.
    global ALIAS_TO_SYMBOL, IG_PREFIXES
    ALIAS_TO_SYMBOL = load_aliases(args.aliases)
    IG_PREFIXES = load_excluded_loci(args.excluded_loci)

    sys.stderr.write(
        f"Loaded {len(panel)} gene symbols from panel BED: {args.panel_bed}; "
        f"{len(ALIAS_TO_SYMBOL)} aliases; "
        f"{len(IG_PREFIXES)} excluded locus prefixes\n")

    summary_rows = []
    for path in args.input:
        if not os.path.isfile(path):
            sys.stderr.write(f"ERROR: input not found: {path}\n")
            sys.exit(1)
        summary_rows.append(
            filter_one_file(path, args.outdir, args.max_pop_af,
                            args.include_ig, panel, args.keep_off_panel))

    os.makedirs(args.outdir, exist_ok=True)
    summary_path = os.path.join(args.outdir, "al_filter_summary.tsv")
    with open(summary_path, "w", newline="") as out:
        cols = ["sample", "input", "kept", "pass", "common_polymorphism",
                "off_panel_reportable"]
        w = csv.DictWriter(out, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in summary_rows:
            w.writerow(r)
    sys.stderr.write(f"\nWrote cohort summary: {summary_path}\n")


if __name__ == "__main__":
    main()
