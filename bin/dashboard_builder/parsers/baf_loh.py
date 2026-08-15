"""Parse the cohort BAF / LOH screen output for a single sample.

Unlike the other parsers in this package, the underlying artefact is
cohort-scoped rather than per-sample: the screen writes one table covering every
sample in the run, because it normalises heterozygous site density for each
panel region against the cohort median for that same region. This parser is
therefore given the run directory as well as the sample directory, locates the
cohort table, and filters it to the sample being reported.

Expected layout::

  <run_dir>/
    hg38/baf_loh/
      cohort.baf_screen.tsv
      figures/
        <sample>.genome_baf_cn.png
        <sample>.region_baf.png
        cohort_baf_deflection_heatmap.png

Screen table columns (from bin/baf_loh_screen.py)::

  sample region chrom start end span_mb n_het het_per_mb het_density_ratio
  median_dp frac_phased median_baf frac_central band_lo band_hi baf_deflection
  depletion_score bimodality flag flag_reason cn cn_event cn_call
  tumour_fraction cn_note

Interpretation note, which the template is expected to surface: a flag and a
copy number call must be read together. Allelic imbalance is produced by
deletion, by gain and by copy-neutral loss of heterozygosity alike, so a
flagged region only indicates copy-neutral LOH when the copy number at that
region is neutral *and* the copy number call is trustworthy. Where the ichorCNA
fit for the sample was degenerate, every copy number call for that sample is
unreliable and the parser surfaces a single sample-level warning rather than
relying on the per-row note being read.

Returns:
  {
    'table':          {columns, rows, n},   # DataTable input, this sample only
    'summary': {
      'n_loh_likely':   int,
      'n_equivocal':    int,
      'n_no_loh':       int,
      'n_unassessable': int,
      'n_regions':      int,
      'tumour_fraction': float | None,
      'top_flagged':    [{region, frac_central, bimodality, cn_call}, ... up to 8]
    },
    'cn_warning':     str | None,           # sample-level ichorCNA fit caveat
    'cn_available':   bool,
    'figures':        [{label, path}],      # paths relative to the sample dir
    'cohort_figure':  {label, path} | None
  }
  or None when no screen output is present.
"""

from pathlib import Path

import pandas as pd


# Flags in the order they should be presented: actionable first.
FLAG_ORDER = ["LOH_LIKELY", "EQUIVOCAL", "UNASSESSABLE", "NO_LOH"]

# Columns surfaced in the report table. The full table carries intermediate
# statistics that are useful for debugging but only clutter a clinical view.
DISPLAY_COLUMNS = [
    "region", "chrom", "start", "end", "n_het", "median_dp",
    "median_baf", "frac_central", "bimodality", "depletion_score",
    "flag", "cn", "cn_call", "cn_note",
]

# Regions whose copy number call cannot be trusted even when present.
UNTRUSTWORTHY_CN_CALLS = {"GAP", "NO_CN"}


def _find_screen_table(run_dir):
    """
    Locate the cohort screen table beneath a run directory.

    The published location is checked first; a recursive search is used as a
    fallback so that the parser still works against directories assembled by
    hand rather than by the workflow.
    """
    if run_dir is None:
        return None
    run_dir = Path(run_dir)

    # Two published layouts are checked directly before falling back to a
    # search: the pipeline results tree, where the screen is published under
    # the hg38 track, and the report bundle, where cohort-level files sit at
    # the bundle root beside filter_summary.tsv.
    for candidate in (
        run_dir / "hg38" / "baf_loh" / "cohort.baf_screen.tsv",
        run_dir / "baf_loh" / "cohort.baf_screen.tsv",
    ):
        if candidate.exists():
            return candidate

    matches = sorted(run_dir.glob("**/cohort.baf_screen.tsv"))
    return matches[0] if matches else None


def _relative_or_none(path, base):
    """Return ``path`` relative to ``base``, or None if it is not beneath it."""
    try:
        return str(Path(path).relative_to(base))
    except ValueError:
        return None


def _collect_figures(screen_path, sample, sample_dir):
    """
    Locate the per-sample and cohort figures beside the screen table.

    Paths are returned relative to the sample directory where possible, since
    the report is written into that directory and loaded over file:// as well
    as http://. Where a figure lies outside the sample directory a path
    relative to the report is still emitted, which the template renders as-is.
    """
    figures = []
    cohort_figure = None

    # Per-sample figures live beside the screen table in the results tree, and
    # under the sample's own directory in the report bundle. Both are checked so
    # the parser works against either layout.
    figure_dirs = [
        screen_path.parent / "figures",
        Path(sample_dir) / "baf_loh",
    ]

    for label, filename in (
        ("Genome-wide BAF and copy number", "{0}.genome_baf_cn.png".format(sample)),
        ("Per-region B-allele frequency", "{0}.region_baf.png".format(sample)),
    ):
        for figure_dir in figure_dirs:
            path = figure_dir / filename
            if path.exists():
                relative = _relative_or_none(path, Path(sample_dir))
                figures.append({"label": label, "path": relative or str(path)})
                break

    figure_dir = screen_path.parent / "figures"
    if not figure_dir.is_dir():
        figure_dir = screen_path.parent

    heatmap = figure_dir / "cohort_baf_deflection_heatmap.png"
    if heatmap.exists():
        relative = _relative_or_none(heatmap, Path(sample_dir))
        cohort_figure = {
            "label": "Cohort allelic imbalance by panel region",
            "path": relative or str(heatmap),
        }

    return figures, cohort_figure


def _flag_sort_key(flag):
    """Sort key placing actionable flags first, unknown flags last."""
    try:
        return FLAG_ORDER.index(flag)
    except ValueError:
        return len(FLAG_ORDER)


def parse(sample_dir, sample, run_dir=None):
    """
    Build the BAF / LOH context for one sample.

    ``run_dir`` is required to locate the cohort table. When it is not supplied
    the sample directory's parents are walked, which covers both the flat and
    the subdirectory publication layouts.
    """
    sample_dir = Path(sample_dir)

    screen_path = _find_screen_table(run_dir)
    if screen_path is None:
        # Walk upwards from the sample directory as a fallback. Two levels
        # covers <run>/<sample>/ and <run>/<sample>/<subdir>/.
        for parent in list(sample_dir.parents)[:3]:
            screen_path = _find_screen_table(parent)
            if screen_path is not None:
                break
    if screen_path is None:
        return None

    try:
        df = pd.read_csv(screen_path, sep="\t", dtype=str,
                         keep_default_na=False, na_values=[""])
    except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError):
        return None

    if "sample" not in df.columns:
        return None

    df = df[df["sample"] == sample].fillna("")
    if df.empty:
        return None

    # Sample-level copy number caveat. The note is identical on every row for a
    # given sample, so it is lifted out and presented once.
    cn_warning = None
    if "cn_note" in df.columns:
        warnings = [note for note in df["cn_note"].unique()
                    if note.startswith("unreliable ichorCNA fit")]
        if warnings:
            cn_warning = warnings[0]

    cn_available = False
    if "cn_call" in df.columns:
        cn_available = bool(set(df["cn_call"].unique()) - UNTRUSTWORTHY_CN_CALLS - {""})

    tumour_fraction = None
    if "tumour_fraction" in df.columns:
        values = pd.to_numeric(df["tumour_fraction"], errors="coerce").dropna()
        if len(values):
            tumour_fraction = float(values.iloc[0])

    flags = df["flag"] if "flag" in df.columns else pd.Series(dtype=str)
    counts = flags.value_counts().to_dict()

    # Order rows so that flagged regions appear first; the template renders the
    # whole set and leaves filtering to the DataTable.
    display = df.copy()
    display["_flag_rank"] = display["flag"].map(_flag_sort_key)
    numeric_fc = pd.to_numeric(display.get("frac_central"), errors="coerce")
    display["_fc"] = numeric_fc.fillna(1.0)
    display = display.sort_values(["_flag_rank", "_fc"]).drop(columns=["_flag_rank", "_fc"])

    columns = [c for c in DISPLAY_COLUMNS if c in display.columns]
    rows = display[columns].to_dict(orient="records")

    # A short list for the summary card, mirroring the coverage parser's
    # low-coverage examples.
    top_flagged = []
    flagged = display[display["flag"] == "LOH_LIKELY"].head(8)
    for record in flagged.to_dict(orient="records"):
        top_flagged.append({
            "region": record.get("region", ""),
            "frac_central": record.get("frac_central", ""),
            "bimodality": record.get("bimodality", ""),
            "cn_call": record.get("cn_call", ""),
        })

    figures, cohort_figure = _collect_figures(screen_path, sample, sample_dir)

    return {
        "table": {"columns": columns, "rows": rows, "n": len(rows)},
        "summary": {
            "n_loh_likely": int(counts.get("LOH_LIKELY", 0)),
            "n_equivocal": int(counts.get("EQUIVOCAL", 0)),
            "n_no_loh": int(counts.get("NO_LOH", 0)),
            "n_unassessable": int(counts.get("UNASSESSABLE", 0)),
            "n_regions": int(len(df)),
            "tumour_fraction": tumour_fraction,
            "top_flagged": top_flagged,
        },
        "cn_warning": cn_warning,
        "cn_available": cn_available,
        "figures": figures,
        "cohort_figure": cohort_figure,
    }
