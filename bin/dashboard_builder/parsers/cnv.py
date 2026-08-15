"""Parse CNV-related outputs.

Per the agreed run layout, a sample's CNV outputs look like::

  <sample>/
    cnv_consensus/
      <sample>_cnv_clinical.tsv          # the clinical CNV table
    cnvkit_plots/
      <sample>.final-scatter.png         # genome-wide scatter (PNG)
      <sample>.final-diagram.pdf         # genome-wide diagram (PDF)
      overview/
        <sample>_genome_wide.png         # full-genome overview
        <sample>_gene_summary_heatmap.png  # panel gene heatmap
      combined/
        <sample>_combined_chr<N>.png     # 23 per-chrom multi-panel views
      per_chromosome/
        <sample>_chr<N>_<GENES>.png      # focused per-chrom views
      per_gene/
        <sample>_gene_<GENE>.png         # one PNG per panel gene

NOTE on filename convention: the cnvkit `batch --scatter --diagram` outputs use
``<sample>.final-...`` (dot separator), not ``<sample>_final-...`` (underscore).
An earlier version of this parser used the underscore form and therefore never
matched the real files. We now glob both forms to be permissive.

The CNV clinical TSV columns we expect (from the pipeline's cnv_annotate.py)::

  gene tier call chromosome arm start end cnvkit_log2 cn_estimate
  zscore LOO_FP_pct cnvkit_confidence concordance arm_level_event
  clinical_significance

Returns a dict with:
  - 'clinical_table' : pass-through DataTable input {columns, rows, n}
  - 'scatter_png'    : str | None (relative to sample_dir)
  - 'diagram_pdf'    : str | None
  - 'overview'       : list[{label, path}]
  - 'combined'       : list[{chrom, path}]
  - 'per_chrom'      : list[{label, path}]    # focused per-chromosome plots
  - 'per_gene_priority' : list[{gene, path}]  # the 8 clinically important genes
  - 'per_gene_other'    : list[{gene, path}]  # everything else (alphabetical)
"""

from pathlib import Path
import re

import pandas as pd


# Clinically important CNV genes for the myeloid panel. These surface at the
# top of the per-gene plot grid in the dashboard's CNV tab. Order here is
# preserved in the rendered grid.
#
# Naming notes:
#   - "PTPN" is the BED-file label used by this panel for the PTPN11 locus
#     (per myeloid_driver_genes.tsv row 87, "BED label likely means PTPN11").
#     Plot filenames follow the BED label.
#   - IKZF interpreted as IKZF1 (IKZF2/IKZF3 are not on the myeloid panel).
#   - ATM is not currently on this panel; included here so it surfaces
#     automatically if the pipeline starts emitting a plot for it. Until
#     then, the gallery footer notes it as a missing priority gene.
PRIORITY_GENES = [
    "KMT2A", "IKZF1", "CDKN2A", "CDKN2B",
    "PTPN", "ETV6", "TET2", "TP53",
    "NF1", "ATM", "RUNX1",
]


# Regex extractors for the per-chrom and per-gene plot filenames.
# We extract the chromosome label / gene name from the filename so the
# template can present a clean caption alongside each image.
_PER_CHROM_RE = re.compile(r"_chr([0-9XY]+)_", re.IGNORECASE)
_GENE_RE      = re.compile(r"_gene_([A-Za-z0-9-]+)\.png$")
_COMBINED_RE  = re.compile(r"_combined_chr([0-9XY]+)\.png$", re.IGNORECASE)


def _chrom_sort_key(label):
    """Sort chromosomes numerically (1, 2, ..., 22) then X, Y."""
    label = str(label or "").upper()
    if label == "X": return (1, 23)
    if label == "Y": return (1, 24)
    if label.isdigit():
        return (0, int(label))
    return (2, label)


def _rel(p, sample_dir):
    try:
        return str(p.relative_to(sample_dir))
    except ValueError:
        return str(p)


def parse(sample_dir, sample):
    sample_dir = Path(sample_dir)
    plots_dir = sample_dir / "cnvkit_plots"

    # ---- Clinical CNV table ----
    clinical_path = sample_dir / "cnv_consensus" / f"{sample}_cnv_clinical.tsv"
    clinical_table = None
    if clinical_path.exists():
        try:
            df = pd.read_csv(
                clinical_path, sep="\t", dtype=str,
                keep_default_na=False, na_values=[""],
            ).fillna("")
            clinical_table = {
                "columns": list(df.columns),
                "rows":    df.to_dict(orient="records"),
                "n":       len(df),
            }
        except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError):
            clinical_table = None

    # ---- Genome-wide scatter PNG / diagram PDF ----
    # Be permissive about dot vs underscore separator.
    scatter_png = None
    for cand in (plots_dir.glob(f"{sample}.final-scatter.png"),
                 plots_dir.glob(f"{sample}_final-scatter.png")):
        for p in cand:
            scatter_png = _rel(p, sample_dir)
            break
        if scatter_png:
            break

    diagram_pdf = None
    for cand in (plots_dir.glob(f"{sample}.final-diagram.pdf"),
                 plots_dir.glob(f"{sample}_final-diagram.pdf")):
        for p in cand:
            diagram_pdf = _rel(p, sample_dir)
            break
        if diagram_pdf:
            break

    # ---- overview/ : whole-panel context plots (typically 1-3 PNGs) ----
    overview = []
    overview_dir = plots_dir / "overview"
    if overview_dir.is_dir():
        for p in sorted(overview_dir.glob("*.png")):
            # Friendly label from the filename suffix after the sample name.
            stem = p.stem
            if stem.startswith(sample + "_"):
                label = stem[len(sample) + 1:].replace("_", " ").strip()
            else:
                label = stem.replace("_", " ")
            overview.append({"label": label, "path": _rel(p, sample_dir)})

    # ---- combined/ : per-chrom multi-panel views ----
    combined = []
    combined_dir = plots_dir / "combined"
    if combined_dir.is_dir():
        for p in sorted(combined_dir.glob("*.png")):
            m = _COMBINED_RE.search(p.name)
            chrom = m.group(1) if m else p.stem
            combined.append({"chrom": chrom, "path": _rel(p, sample_dir)})
        combined.sort(key=lambda x: _chrom_sort_key(x["chrom"]))

    # ---- per_chromosome/ : focused per-chrom views (with gene list in name) ----
    per_chrom = []
    per_chrom_dir = plots_dir / "per_chromosome"
    if per_chrom_dir.is_dir():
        for p in sorted(per_chrom_dir.glob("*.png")):
            m = _PER_CHROM_RE.search(p.name)
            chrom = m.group(1) if m else "?"
            # The full filename without sample prefix and .png suffix makes a
            # readable label (e.g. "chr11 HRAS WT1 SF1 EED KMT2A +1").
            stem = p.stem
            if stem.startswith(sample + "_"):
                label = stem[len(sample) + 1:].replace("_", " ")
            else:
                label = stem.replace("_", " ")
            per_chrom.append({
                "chrom": chrom,
                "label": label,
                "path":  _rel(p, sample_dir),
            })
        per_chrom.sort(key=lambda x: _chrom_sort_key(x["chrom"]))

    # ---- per_gene/ : split into priority + other ----
    per_gene_all = {}
    per_gene_dir = plots_dir / "per_gene"
    if per_gene_dir.is_dir():
        for p in sorted(per_gene_dir.glob("*.png")):
            m = _GENE_RE.search(p.name)
            if not m:
                continue
            gene = m.group(1)
            per_gene_all[gene] = _rel(p, sample_dir)

    per_gene_priority = []
    for g in PRIORITY_GENES:
        if g in per_gene_all:
            per_gene_priority.append({"gene": g, "path": per_gene_all[g]})

    priority_set = set(PRIORITY_GENES)
    per_gene_other = [
        {"gene": g, "path": per_gene_all[g]}
        for g in sorted(per_gene_all)
        if g not in priority_set
    ]

    # ---- Primary genome-wide scatter (added 2026-05-24) ----
    # The clean genome-wide scatter is the primary view for arm-level
    # CNV events. The parser pulls it out of the generic overview list
    # into its own slot so the template can render it at full width as
    # the first thing in the CNV tab.
    primary_genome_scatter = None
    clean_scatter_path = plots_dir / "overview" / f"{sample}_genome_scatter_clean.png"
    if clean_scatter_path.exists():
        primary_genome_scatter = _rel(clean_scatter_path, sample_dir)
        # Remove from the generic overview list to avoid duplicate rendering
        overview = [item for item in overview
                    if item.get("path") != primary_genome_scatter]

    return {
        "clinical_table":     clinical_table,
        "scatter_png":        scatter_png,
        "primary_genome_scatter": primary_genome_scatter,
        "diagram_pdf":        diagram_pdf,
        "overview":           overview,
        "combined":           combined,
        "per_chrom":          per_chrom,
        "per_gene_priority":  per_gene_priority,
        "per_gene_other":     per_gene_other,
        # Backward-compat aliases so any code still referencing the old keys
        # doesn't crash. The template will be updated to use the new keys.
        "per_chrom_pngs":     [item["path"] for item in per_chrom],
        "per_gene_pngs":      [item["path"] for item in per_gene_priority]
                              + [item["path"] for item in per_gene_other],
    }
