#!/usr/bin/env python3
"""Build static HTML dashboard for a targeted-seq pipeline run directory.

Usage:
    python build.py /path/to/run_directory
    python build.py /path/to/run_directory --subdir clinical

Produces, alongside the existing files:
    <run_dir>/cohort_index.html
    <run_dir>/assets/                       (copied from the builder's vendored set)
    <run_dir>/<sample>/<subdir>/<sample>_report.html (one per sample subdir)

If --subdir is given (e.g. ``clinical``), every per-sample file lookup, every
per-sample write (the report HTML, annotation caches, IGV hash-router patch),
and every URL in the cohort index treats ``<sample>/<subdir>/`` as the sample's
working directory instead of ``<sample>/`` directly. This accommodates pipeline
layouts (notably the nf-core port) that publish per-sample outputs into a named
subdirectory under each sample.

The builder is otherwise read-only against existing files in the run directory,
with one documented exception: <sample>_igv_report.html is patched in place to
add an idempotent hash-router <script> at the bottom of <body>. This lets the
dashboard select a variant in the IGV report by URL fragment ("#row_<uid>"),
which is the only mechanism that works under both http:// and file:// loads.
The injection is delimited by HTML comment sentinels and can be re-run safely.
"""

import argparse
import logging
import re
import shutil
import sys
from datetime import datetime
from pathlib import Path

from jinja2 import Environment, FileSystemLoader, select_autoescape

# Make `parsers` importable when running this script directly.
sys.path.insert(0, str(Path(__file__).resolve().parent))

from parsers import hsmetrics as p_hsmetrics
from parsers import variants as p_variants
from parsers import flt3 as p_flt3
from parsers import coverage as p_coverage
from parsers import cnv as p_cnv
from parsers import igv as p_igv
from parsers import genebe as p_genebe
from parsers import mobidetails as p_mobidetails
from parsers import oncokb as p_oncokb
from parsers import cancervar as p_cancervar
from parsers import baf_loh as p_baf_loh
from parsers import translocations as p_translocations
from parsers import ichor as p_ichor
from parsers import qc as p_qc


BUILDER_VERSION = "0.5.0-triage+exon"

# Directories under a run dir that are NOT samples.
# baf_loh holds the cohort-level BAF/LOH screen table and heatmap rather than
# a sample. It sits at the bundle root for the same reason filter_summary.tsv
# does: the screen judges each region against the cohort, so one table covers
# the whole run and there is nothing per-sample to report from it.
NON_SAMPLE_DIRS = {"pipeline_info", "assets", "baf_loh"}

# File descriptions for the Files tab. Suffix -> human-readable.
# Both the production naming convention (underscores throughout) and the
# nf-core port's dot-separated convention are listed so the Files tab
# labels them correctly regardless of which pipeline produced the run.
FILE_DESCRIPTIONS = [
    ("_hsmetrics.txt",                       "Picard CollectHsMetrics output"),
    ("_somaticseq_clinical_final.tsv",       "Curated clinical variants"),
    (".somaticseq.clinical.final.tsv",       "Curated clinical variants"),
    ("_somaticseq_filtered.tsv",             "Full annotated variants"),
    (".somaticseq.filtered.tsv",             "Full annotated variants"),
    ("_flt3_consensus.tsv",                  "FLT3-ITD consensus calls"),
    ("_exon_coverage.tsv",                   "Per-exon coverage (duplicates included)"),
    ("_igv_report.html",                     "IGV visualization (clinical variants)"),
    ("_focal_dup_igv.html",                  "IGV visualization (focal duplications)"),
    ("_NV_fastp.html",                       "fastp QC report"),
    ("_fastp.html",                          "fastp QC report"),
    ("_dashboard.html",                      "Pre-existing per-sample dashboard"),
    ("_genebe_cache.json",                   "GeneBe annotation cache (audit)"),
    ("_mobidetails_cache.json",              "MobiDetails lookup cache (audit, not surfaced in UI)"),
    ("_oncokb_cache.json",                   "OncoKB annotation cache (audit)"),
    ("_cancervar_cache.json",                "CancerVar (AMP/ASCO/CAP) annotation cache (audit)"),
]


def describe_file(filename):
    """Return a brief description for a per-sample file based on its suffix."""
    for suffix, desc in FILE_DESCRIPTIONS:
        if filename.endswith(suffix):
            return desc
    return ""


def discover_samples(run_dir):
    """Return sorted list of sample subdirectories.

    A subdirectory is a sample if it is not in NON_SAMPLE_DIRS and is a directory.
    We do not require a particular sentinel file here so missing-file behavior
    is handled per-parser later.
    """
    samples = []
    for entry in sorted(run_dir.iterdir()):
        if not entry.is_dir():
            continue
        if entry.name in NON_SAMPLE_DIRS:
            continue
        samples.append(entry)
    return samples


# Exon token from a panel-BED column-4 label, e.g.
#   "Target=1;ProbeIdx=12;GNB1_Ex_11"      -> "11"
#   "926537_53648240_ARID1A_Ex_1_1"        -> "1"
#   "WT1_Ex_1A"                            -> "1A"
# Captures digits plus an optional single trailing letter (1A/1B/8A exist in
# the myeloid panel); the probe-index suffix (_1) after that is ignored.
_BED_EXON_RE = re.compile(r"_Ex_?([0-9]+[A-Za-z']?)", re.IGNORECASE)


def _load_panel_bed_exons(path):
    """Build {chrom: [(start, end, exon_token), ...]} from a panel BED.

    Returns None if the path is falsy or unreadable so callers can skip the
    BED step cleanly. BED is 0-based half-open; we keep coordinates as-is and
    test with start <= pos < end after converting the variant's 1-based POS.
    """
    if not path:
        return None
    p = Path(path)
    if not p.exists():
        logging.warning("--panel-bed not found, skipping BED exon fallback: %s", path)
        return None
    index = {}
    try:
        with open(p) as fh:
            for line in fh:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 4:
                    continue
                chrom = parts[0]
                try:
                    start = int(parts[1])
                    end = int(parts[2])
                except ValueError:
                    continue
                m = _BED_EXON_RE.search(parts[3])
                if not m:
                    continue
                index.setdefault(chrom, []).append((start, end, m.group(1)))
    except OSError as exc:
        logging.warning("could not read --panel-bed (%s): %s", path, exc)
        return None
    logging.info("Loaded panel BED exon index: %d contigs", len(index))
    return index


def _bed_exon_for(bed_index, chrom, pos_1based):
    """Return the exon token whose BED interval contains the 1-based position."""
    if not bed_index or not chrom:
        return ""
    try:
        pos0 = int(str(pos_1based)) - 1
    except (TypeError, ValueError):
        return ""
    for start, end, exon in bed_index.get(chrom, ()):
        if start <= pos0 < end:
            return exon
    return ""


def _report_exon(row, bed_index):
    """Effective reporting exon for one variant row.

    Precedence: VariantValidator (VV_Exon) -> panel BED -> VEP EXON. Each source
    is treated as missing when blank or the pipeline's "-1" sentinel. VEP's EXON
    is "<n>/<total>" (e.g. "8/52"); we keep only the exon number.
    """
    vv = str(row.get("VV_Exon", "") or "").strip()
    if vv and vv != "-1":
        return vv
    bed = _bed_exon_for(bed_index, row.get("Chr", ""), row.get("Start", ""))
    if bed:
        return bed
    vep = str(row.get("EXON", "") or "").strip()
    if vep and vep != "-1":
        return vep.split("/")[0]
    return ""


def _apply_report_exon(parsed, bed_index):
    """Add an EXON_REPORT key to every row of a parsed-variant dict, in place."""
    if not parsed or not parsed.get("rows"):
        return
    for row in parsed["rows"]:
        row["EXON_REPORT"] = _report_exon(row, bed_index)
    if "EXON_REPORT" not in parsed["columns"]:
        parsed["columns"].append("EXON_REPORT")


def collect_sample_context(sample_dir, build_time, subdir="",
                           annotate_genebe=False, genebe_user=None, genebe_key=None,
                           annotate_mobidetails=False,
                           annotate_oncokb=False, oncokb_token=None,
                           annotate_cancervar=False, panel_bed_index=None):
    """Run every parser on one sample directory; return a context dict for the template.

    ``subdir`` (default ``""``) is an optional path component appended to
    ``sample_dir`` to obtain the directory where per-sample files actually
    live. With ``subdir=""`` the effective directory equals ``sample_dir``
    (v0.4.4 behavior). With ``subdir="clinical"`` everything reads from and
    writes to ``sample_dir / "clinical"`` instead. The sample identifier
    (used in filenames and as the display name) is always derived from
    ``sample_dir.name``, regardless of subdir.
    """
    sample = sample_dir.name
    effective_dir = (sample_dir / subdir) if subdir else sample_dir

    ctx = {
        "sample": sample,
        "build_time": build_time,
        "builder_version": BUILDER_VERSION,
        "hsmetrics": None,
        "clinical": None,
        "filtered": None,
        "flt3": None,
        "coverage": None,
        "cnv": None,
        "igv": None,
        "genebe": {},        # chr:pos:ref:alt -> annotation dict (empty unless --annotate-genebe)
        "mobidetails": {},   # chr:pos:ref:alt -> {url, mobidetails_id, ...} (empty unless --annotate-mobidetails)
        "oncokb": {},        # chr:pos:ref:alt -> OncoKB annotation dict (empty unless --annotate-oncokb)
        "cancervar": {},     # chr:pos:ref:alt -> CancerVar annotation dict (empty unless --annotate-cancervar)
        "ichor_sol_pages": [],
        "files": {"fastp": None, "igv_report": None, "focal_dup_igv": None,
                  "existing_dashboard": None, "listing": []},
    }

    # --- HsMetrics ---
    hsmetrics_path = effective_dir / f"{sample}_hsmetrics.txt"
    try:
        ctx["hsmetrics"] = p_hsmetrics.parse(hsmetrics_path)
    except Exception as exc:
        logging.warning("[%s] hsmetrics parse failed: %s", sample, exc)

    # --- Clinical variants ---
    # Production names this <sample>_somaticseq_clinical_final.tsv; the nf-core
    # port's ORGANIZE_OUTPUT publishes <sample>.somaticseq.clinical.final.tsv.
    # Try both, in that order.
    clinical_path = None
    for cand in (f"{sample}_somaticseq_clinical_final.tsv",
                 f"{sample}.somaticseq.clinical.final.tsv"):
        p = effective_dir / cand
        if p.exists():
            clinical_path = p
            break
    if clinical_path is None:
        # Pick the canonical name so the WARNING log is informative.
        clinical_path = effective_dir / f"{sample}_somaticseq_clinical_final.tsv"
    try:
        ctx["clinical"] = p_variants.parse(clinical_path)
        _apply_report_exon(ctx["clinical"], panel_bed_index)
    except Exception as exc:
        logging.warning("[%s] clinical variants parse failed: %s", sample, exc)

    # --- Filtered variants ---
    # Same naming-convention story as the clinical TSV above.
    filtered_path = None
    for cand in (f"{sample}_somaticseq_filtered.tsv",
                 f"{sample}.somaticseq.filtered.tsv"):
        p = effective_dir / cand
        if p.exists():
            filtered_path = p
            break
    if filtered_path is None:
        filtered_path = effective_dir / f"{sample}_somaticseq_filtered.tsv"
    try:
        ctx["filtered"] = p_variants.parse(filtered_path)
        _apply_report_exon(ctx["filtered"], panel_bed_index)
    except Exception as exc:
        logging.warning("[%s] filtered variants parse failed: %s", sample, exc)

    # --- FLT3 ---
    flt3_path = effective_dir / f"{sample}_flt3_consensus.tsv"
    try:
        ctx["flt3"] = p_flt3.parse(flt3_path)
    except Exception as exc:
        logging.warning("[%s] flt3 parse failed: %s", sample, exc)

    # --- Exon coverage ---
    cov_path = effective_dir / f"{sample}_exon_coverage.tsv"
    try:
        ctx["coverage"] = p_coverage.parse(cov_path)
    except Exception as exc:
        logging.warning("[%s] coverage parse failed: %s", sample, exc)

    # --- CNV ---
    try:
        ctx["cnv"] = p_cnv.parse(effective_dir, sample)

        # The BAF / LOH screen writes one table for the whole run rather than
        # one per sample, so the parser is given the run directory and filters
        # to this sample itself. The parser walks upwards from the directory it
        # is given if the table does not sit directly beneath it.
        ctx["baf_loh"] = p_baf_loh.parse(effective_dir, sample,
                                         run_dir=sample_dir.parent)
    except Exception as exc:
        logging.warning("[%s] cnv parse failed: %s", sample, exc)
        ctx["cnv"] = {"clinical_table": None, "scatter_png": None, "diagram_pdf": None,
                      "per_chrom_pngs": [], "per_gene_pngs": []}

    # --- IGV lookup + idempotent hash-router injection ---
    # The hash-router script lets the parent dashboard select a variant by
    # setting iframe.src to "...#row_<uid>", which works under file:// where
    # reading the iframe's DOM is blocked by cross-origin policy. We modify
    # the IGV report file in place (idempotent via sentinel comments) since
    # writing a 6+ MB patched copy alongside the original is wasteful.
    # Focal-duplication IGV page. Separate from the clinical SNV report because
    # the duplications are called on a different table and reviewed on their own
    # terms: what matters is whether the reads at that position actually carry a
    # tandem duplication, or whether the aligner put them there.
    fd_igv = effective_dir / f"{sample}_focal_dup_igv.html"
    if fd_igv.exists():
        ctx["files"]["focal_dup_igv"] = fd_igv.name

    # Rasterised ichorCNA solution pages, if the bundle produced them.
    sol_dir = effective_dir / "cnv" / "all_sols"
    if not sol_dir.is_dir():
        sol_dir = effective_dir / "all_sols"
    pages = []
    if sol_dir.is_dir():
        for img in sorted(sol_dir.glob(f"{sample}.sol*.png")):
            n = "".join(ch for ch in img.stem.split("sol")[-1] if ch.isdigit())
            pages.append({"name": n or img.stem,
                          "label": f"Solution page {n or '?'}",
                          "path": str(img.relative_to(effective_dir))})
    ctx["ichor_sol_pages"] = pages

    igv_path = effective_dir / f"{sample}_igv_report.html"
    if igv_path.exists():
        ctx["files"]["igv_report"] = igv_path.name
        try:
            ctx["igv"] = p_igv.extract_lookup(igv_path)
        except Exception as exc:
            logging.warning("[%s] igv lookup extraction failed: %s", sample, exc)
        try:
            if p_igv.inject_hash_router(igv_path):
                logging.info("[%s] Injected hash-router into IGV report.", sample)
        except Exception as exc:
            logging.warning("[%s] IGV hash-router injection failed: %s", sample, exc)

    # --- fastp HTML ---
    for cand in (f"{sample}_NV_fastp.html", f"{sample}_fastp.html"):
        if (effective_dir / cand).exists():
            ctx["files"]["fastp"] = cand
            break

    # --- Existing pipeline-side dashboard (QC tab embed) ---
    existing = effective_dir / f"{sample}_dashboard.html"
    if existing.exists():
        ctx["files"]["existing_dashboard"] = existing.name

    # --- Optional: GeneBe build-time annotation of clinical variants ---
    if annotate_genebe and ctx["clinical"] and ctx["clinical"].get("rows"):
        try:
            ctx["genebe"] = p_genebe.annotate(
                clinical_rows=ctx["clinical"]["rows"],
                sample_dir=effective_dir,
                sample=sample,
                api_user=genebe_user,
                api_key=genebe_key,
            )
        except Exception as exc:
            logging.warning("[%s] GeneBe annotation crashed: %s", sample, exc)
            ctx["genebe"] = {}

    # --- Optional: MobiDetails build-time lookup of clinical variants ---
    if annotate_mobidetails and ctx["clinical"] and ctx["clinical"].get("rows"):
        try:
            md_cache = p_mobidetails.annotate(
                clinical_rows=ctx["clinical"]["rows"],
                sample_dir=effective_dir,
                sample=sample,
            )
            # Embed only the hits in the page JSON; not-found entries stay on disk only.
            ctx["mobidetails"] = p_mobidetails.filter_for_frontend(md_cache)
        except Exception as exc:
            logging.warning("[%s] MobiDetails annotation crashed: %s", sample, exc)
            ctx["mobidetails"] = {}

    # --- Optional: OncoKB build-time annotation of clinical variants ---
    if annotate_oncokb and ctx["clinical"] and ctx["clinical"].get("rows"):
        try:
            ctx["oncokb"] = p_oncokb.annotate(
                clinical_rows=ctx["clinical"]["rows"],
                sample_dir=effective_dir,
                sample=sample,
                oncokb_token=oncokb_token,
            )
        except Exception as exc:
            logging.warning("[%s] OncoKB annotation crashed: %s", sample, exc)
            ctx["oncokb"] = {}

    # --- Optional: CancerVar build-time annotation of clinical variants ---
    # AMP/ASCO/CAP tier classification + 12 CBP criteria + OPAI score.
    # Anonymous endpoint, no token required.
    if annotate_cancervar and ctx["clinical"] and ctx["clinical"].get("rows"):
        try:
            cv_cache = p_cancervar.annotate(
                clinical_rows=ctx["clinical"]["rows"],
                sample_dir=effective_dir,
                sample=sample,
            )
            # Negative cache entries stay on disk (so we don't re-query them)
            # but we strip them out of the JSON embedded in the HTML.
            ctx["cancervar"] = p_cancervar.filter_for_frontend(cv_cache)
        except Exception as exc:
            logging.warning("[%s] CancerVar annotation crashed: %s", sample, exc)
            ctx["cancervar"] = {}

    # --- Files tab listing (everything in the sample dir we recognize) ---
    listing = []
    for entry in sorted(effective_dir.iterdir()):
        if entry.is_file():
            listing.append((entry.name, describe_file(entry.name)))
    ctx["files"]["listing"] = listing

    # Structural variants, copy number and adaptive-sampling QC.
    # Each parser returns None when its inputs are absent, and the
    # matching template renders an empty state rather than failing.
    ctx["translocations"] = p_translocations.parse(effective_dir, sample)
    ctx["ichor"] = p_ichor.parse(effective_dir, sample)
    ctx["qc"] = p_qc.parse(effective_dir, sample)
    return ctx


def copy_assets(builder_dir, run_dir):
    """Copy the vendored assets/ directory into the run dir.

    We overwrite (idempotent) so the run dir always has the assets that match
    the builder version that produced the HTML.
    """
    src = builder_dir / "assets"
    dst = run_dir / "assets"
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst)
    logging.info("Copied assets to %s", dst)


def format_num_filter(value, ndigits=2):
    """Jinja2 filter: format a number to N decimals, or empty string if not numeric."""
    if value is None or value == "":
        return ""
    try:
        return f"{float(value):.{ndigits}f}"
    except (TypeError, ValueError):
        return str(value)


def build(run_dir, subdir="",
          annotate_genebe=False, genebe_user=None, genebe_key=None,
          annotate_mobidetails=False,
          annotate_oncokb=False, oncokb_token=None,
          annotate_cancervar=False, panel_bed=None):
    run_dir = Path(run_dir).resolve()
    if not run_dir.is_dir():
        raise SystemExit(f"Not a directory: {run_dir}")

    # Normalize and validate subdir. Forbid leading/trailing slashes, ``..``,
    # and absolute paths; ``Path`` would otherwise let us escape sample_dir.
    subdir = subdir.strip("/") if subdir else ""
    if subdir:
        if ".." in Path(subdir).parts or Path(subdir).is_absolute():
            raise SystemExit(f"--subdir must be a simple relative path, got: {subdir!r}")

    builder_dir = Path(__file__).resolve().parent

    env = Environment(
        loader=FileSystemLoader(str(builder_dir / "templates")),
        autoescape=select_autoescape(["html"]),
        trim_blocks=True,
        lstrip_blocks=True,
    )
    env.filters["format_num"] = format_num_filter
    macros = env.get_template("macros.html.j2").module

    build_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Panel BED exon index (used as the fallback exon source when a variant has
    # no VV_Exon). Loaded once and shared across samples.
    panel_bed_index = _load_panel_bed_exons(panel_bed)

    # Per-sample report sits at <run_dir>/<sample>/<subdir>/<sample>_report.html
    # The depth below run_dir is 1 (sample dir) + however many parts in subdir.
    # ``assets/`` always lives at <run_dir>/assets, so the relative prefix needs
    # exactly that many ``../`` segments.
    sample_depth = 1 + (len(Path(subdir).parts) if subdir else 0)
    sample_assets_prefix = "/".join([".."] * sample_depth) + "/assets"
    # Relative link from a sample report back to <run_dir>/cohort_index.html.
    # Same depth as assets_prefix; just swap the final segment.
    sample_cohort_link = "/".join([".."] * sample_depth) + "/cohort_index.html"
    # For the cohort_index template, the URL fragment between <sample>/ and
    # <sample>_report.html is either empty (no subdir) or "<subdir>/".
    report_subdir = (subdir + "/") if subdir else ""

    # 1. Copy assets first so per-sample reports can already reference them.
    copy_assets(builder_dir, run_dir)

    # 2. Walk samples.
    sample_dirs = discover_samples(run_dir)
    if not sample_dirs:
        logging.warning("No sample subdirectories found in %s", run_dir)

    sample_contexts = []
    sample_template = env.get_template("sample_report.html.j2")

    for sample_dir in sample_dirs:
        ctx = collect_sample_context(
            sample_dir, build_time,
            subdir=subdir,
            annotate_genebe=annotate_genebe,
            genebe_user=genebe_user,
            genebe_key=genebe_key,
            annotate_mobidetails=annotate_mobidetails,
            annotate_oncokb=annotate_oncokb,
            oncokb_token=oncokb_token,
            annotate_cancervar=annotate_cancervar,
            panel_bed_index=panel_bed_index,
        )
        out_dir = (sample_dir / subdir) if subdir else sample_dir
        if not out_dir.is_dir():
            logging.warning("[%s] expected subdir %s/ not found; skipping sample report",
                            sample_dir.name, subdir)
            continue
        out_path = out_dir / f"{sample_dir.name}_report.html"
        with open(out_path, "w", encoding="utf-8") as fh:
            fh.write(sample_template.render(ctx=ctx, macros=macros,
                                            assets_prefix=sample_assets_prefix,
                                            cohort_link=sample_cohort_link))
        logging.info("Wrote %s", out_path)
        sample_contexts.append(ctx)

    # 3. Cohort index.
    cohort_template = env.get_template("cohort_index.html.j2")
    cohort_ctx = {
        "run_name": run_dir.name,
        "build_time": build_time,
        "samples": sample_contexts,
        "report_subdir": report_subdir,
    }
    cohort_path = run_dir / "cohort_index.html"
    with open(cohort_path, "w", encoding="utf-8") as fh:
        fh.write(cohort_template.render(ctx=cohort_ctx, assets_prefix="assets"))
    logging.info("Wrote %s", cohort_path)

    return cohort_path


def main():
    parser = argparse.ArgumentParser(description="Build static dashboard for a targeted-seq run directory.")
    parser.add_argument("run_dir", help="Path to the run directory (e.g. /path/to/20260521_133954)")
    parser.add_argument(
        "--subdir", default="",
        help="Optional path component under each sample directory where per-sample "
             "files live (e.g. 'clinical' for the nf-core pipeline layout). When set, "
             "every file lookup, the in-place IGV hash-router patch, the per-sample "
             "report write, and cohort_index URLs all use <sample>/<subdir>/ as the "
             "sample's working directory. Default is empty (v0.4.4 behavior)."
    )
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")
    parser.add_argument(
        "--annotate-genebe", action="store_true",
        help="Annotate clinical variants via the GeneBe REST API at build time (network required)."
    )
    parser.add_argument(
        "--genebe-user", default=None,
        help="Email for GeneBe Basic auth (recommended; raises rate limits). "
             "Used only when --annotate-genebe is set."
    )
    parser.add_argument(
        "--genebe-key", default=None,
        help="API key from your GeneBe profile page. "
             "Used only when --annotate-genebe is set."
    )
    parser.add_argument(
        "--annotate-mobidetails", action="store_true",
        help="Resolve clinical variants to MobiDetails record IDs at build time "
             "via the /api/variant/exists endpoint (network required, no API key). "
             "Writes <sample>_mobidetails_cache.json with which variants are known "
             "to MD. Note: this does NOT add UI links -- MD has no reliable "
             "anonymous deep link, so the cache is for audit purposes only. "
             "The user-facing path is the Copy VV_HGVS dropdown."
    )
    parser.add_argument(
        "--annotate-oncokb", action="store_true",
        help="Annotate clinical variants via OncoKB at build time. Requires "
             "--oncokb-token (free academic token from "
             "https://www.oncokb.org/account/register). Cached per-sample."
    )
    parser.add_argument(
        "--oncokb-token", default=None,
        help="OncoKB Bearer token. Used only when --annotate-oncokb is set."
    )
    parser.add_argument(
        "--annotate-cancervar", action="store_true",
        help="Annotate clinical variants via CancerVar at build time (AMP/ASCO/"
             "CAP 2017 tier classification + 12 CBP criteria + OPAI score). "
             "Anonymous endpoint, no token required. Cached per-sample. When "
             "present, the Reporting tab's tier column is pre-filled with the "
             "CancerVar tier for newly-selected variants."
    )
    parser.add_argument(
        "--panel-bed", default=None,
        help="Panel BED whose column-4 labels carry exon names (e.g. "
             "GNB1_Ex_11). Used as the fallback exon source for the Reporting "
             "tab and the per-variant Copy Exon action when a variant has no "
             "VariantValidator VV_Exon. Optional; without it, exon falls back "
             "to VEP's EXON column."
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s  %(levelname)-7s %(message)s",
    )

    cohort_path = build(
        args.run_dir,
        subdir=args.subdir,
        annotate_genebe=args.annotate_genebe,
        genebe_user=args.genebe_user,
        genebe_key=args.genebe_key,
        annotate_mobidetails=args.annotate_mobidetails,
        annotate_oncokb=args.annotate_oncokb,
        oncokb_token=args.oncokb_token,
        annotate_cancervar=args.annotate_cancervar,
        panel_bed=args.panel_bed,
    )
    print(f"\nDashboard built. Open:\n  {cohort_path}\n")


if __name__ == "__main__":
    main()
