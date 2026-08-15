"""Optional build-time CancerVar annotation of clinical variants.

CancerVar (https://cancervar.wglab.org) implements the AMP/ASCO/CAP 2017
guideline for clinical interpretation of somatic cancer variants. It returns,
per variant:

  * a four-level Tier classification (Tier I strong / Tier II potential /
    Tier III unknown / Tier IV benign), which maps directly onto the
    AMP/ASCO/CAP framework that pathologists report against;
  * 12 CBP (Clinical-Based Prediction) criteria scores, each 0/1/2;
  * an OPAI ("Oncogenic Pathogenicity AI") score in [0, 1], derived from a
    semi-supervised deep-learning model trained on curated cancer variants.

Endpoint (anonymous, no API key):
    GET https://cancervar.wglab.org/api_new.php
        ?queryType=position&build=hg38&chr=<C>&pos=<P>&ref=<R>&alt=<A>

A successful response is a single flat JSON object, e.g.::

    {"Cancervar":"10#Tier_II_potential","Build":"hg38","Chromosome":2,
     "Position":25234373,"Ref_allele":"C","Alt_allele":"T","Gene":"DNMT3A",
     "CBP_1":1,...,"CBP_12":1,"OPAI":0.54}

Unknown / unparseable variants return an empty body (zero bytes) with HTTP
200, which we treat as "no annotation" and still record in the cache to avoid
re-querying every build.

This module is OPT-IN. The builder only calls it when --annotate-cancervar
is passed. Annotations are cached on disk at::

    <sample_dir>/<sample>_cancervar_cache.json

keyed by ``chr:pos:ref:alt`` (matching the GeneBe/OncoKB cache key scheme).
A negative entry is stored as ``{"_found": False, "_fetched_at": ...}``.

Reference:
    Li Q, Ren Z, Cao K, Li MM, Zhou Y, Wang K. CancerVar: an Artificial
    Intelligence empowered platform for clinical interpretation of somatic
    mutations in cancer. Sci Adv. 2022.
"""

from pathlib import Path
import json
import logging
import time
from datetime import datetime


CANCERVAR_URL = "https://cancervar.wglab.org/api_new.php"
TIMEOUT_S = 30
SLEEP_BETWEEN_REQUESTS_S = 0.5  # be polite -- the API is on older PHP/Apache infra

# CBP criteria labels per Li et al 2022 (Sci Adv) supplementary methods.
# Each criterion is scored 0 / 1 / 2 by the CancerVar engine, where 2 is
# strong evidence supporting the criterion. The display label is shown in
# the dashboard's detail panel next to the per-criterion score.
CBP_DESCRIPTIONS = {
    "CBP_1":  "FDA-approved or well-established therapeutic association (same cancer type)",
    "CBP_2":  "FDA-approved or well-established therapeutic association (different cancer type)",
    "CBP_3":  "Evidence from well-powered studies with consensus (investigational biomarker)",
    "CBP_4":  "Evidence from multiple small studies (investigational biomarker, same cancer type)",
    "CBP_5":  "Preclinical studies only",
    "CBP_6":  "No evidence for oncogenic function",
    "CBP_7":  "Listed in cancer mutation hotspot databases",
    "CBP_8":  "Observed in functional domain with established cancer mechanism",
    "CBP_9":  "Population allele frequency (0 = common, 1 = rare, 2 = absent/very rare)",
    "CBP_10": "Predicted damaging by computational tools",
    "CBP_11": "Listed as somatic mutation in cancer databases (COSMIC, cBioPortal)",
    "CBP_12": "Literature-reported as somatic mutation in cancer",
}

# Map the slug returned by the API to a clean human label we surface in the UI.
# The live API has been observed to emit the tier names with varying case
# ("Tier_III_Uncertain" vs "Tier_III_unknown" in some downstream docs), so the
# helper normalises the third token before lookup. Anything we don't recognise
# falls back to the raw slug.
TIER_LABELS = {
    "Tier_I":   "Tier I",
    "Tier_II":  "Tier II",
    "Tier_III": "Tier III",
    "Tier_IV":  "Tier IV",
}

TIER_DESCRIPTIONS = {
    "Tier_I":   "Variants of Strong Clinical Significance (FDA-approved or well-established)",
    "Tier_II":  "Variants of Potential Clinical Significance (investigational/preclinical evidence)",
    "Tier_III": "Variants of Unknown Clinical Significance",
    "Tier_IV":  "Benign or Likely Benign Variants",
}


def _tier_root(slug):
    """Reduce a tier slug like ``"Tier_III_Uncertain"`` to its leading
    Roman-numeral root ``"Tier_III"`` so labels resolve regardless of how the
    server cased the third token. Returns ``""`` if the slug is missing or
    doesn't match the expected pattern."""
    if not slug or not isinstance(slug, str):
        return ""
    parts = slug.split("_")
    if len(parts) < 2 or parts[0] != "Tier":
        return ""
    return "Tier_" + parts[1]


def _variant_key(chrom, pos, ref, alt):
    return f"{chrom}:{pos}:{ref}:{alt}"


def _strip_chr(chrom):
    s = str(chrom or "")
    return s[3:] if s.startswith("chr") else s


def _parse_cancervar_field(value):
    """Split CancerVar's ``"10#Tier_II_potential"`` into (score, slug).

    Returns ``(None, None)`` if the field is missing or malformed. The score
    is an integer the CancerVar engine emits internally; we don't use it for
    classification (the slug already carries that), but we surface it for audit.
    """
    if not value or not isinstance(value, str) or "#" not in value:
        return (None, None)
    score_part, _, slug = value.partition("#")
    try:
        score = int(score_part)
    except (TypeError, ValueError):
        score = None
    return (score, slug or None)


def _summarise(record):
    """Reduce a raw CancerVar response to the fields the dashboard renders."""
    out = {"_found": True}

    score, slug = _parse_cancervar_field(record.get("Cancervar"))
    out["tier_score"]       = score                                  # internal CancerVar score, e.g. 10
    out["tier_slug"]        = slug                                   # e.g. "Tier_II_potential" (raw, preserved)
    tier_root = _tier_root(slug)
    out["tier_label"]       = TIER_LABELS.get(tier_root, slug or "") # e.g. "Tier II"  (used to prefill the report)
    out["tier_description"] = TIER_DESCRIPTIONS.get(tier_root, "")

    out["gene"]             = record.get("Gene")
    out["build"]            = record.get("Build")
    out["opai"]             = record.get("OPAI")             # 0.0 - 1.0

    # Per-criterion 0/1/2 scores, in their canonical order. We carry the labels
    # too so the frontend can render the table without duplicating them.
    cbp = []
    for k in ("CBP_1", "CBP_2", "CBP_3", "CBP_4", "CBP_5", "CBP_6",
              "CBP_7", "CBP_8", "CBP_9", "CBP_10", "CBP_11", "CBP_12"):
        if k in record:
            cbp.append({
                "id":          k,
                "score":       record.get(k),
                "description": CBP_DESCRIPTIONS.get(k, ""),
            })
    out["cbp"] = cbp
    return out


def annotate(clinical_rows, sample_dir, sample):
    """Annotate clinical variants via CancerVar byPosition.

    Returns ``{chr:pos:ref:alt -> annotation dict}``. The dict carries
    ``_found: True`` for hits and ``_found: False`` for negative cache entries.
    Per-variant failures (network, parse) are logged and don't abort the pass.

    No token / API key required -- the endpoint is anonymous.
    """
    try:
        import requests
    except ImportError:
        logging.warning(
            "[%s] CancerVar annotation requested but 'requests' is not installed.",
            sample,
        )
        return {}

    cache_path = Path(sample_dir) / f"{sample}_cancervar_cache.json"
    cache = {}
    if cache_path.exists():
        try:
            with open(cache_path, "r", encoding="utf-8") as fh:
                cache = json.load(fh)
        except (OSError, json.JSONDecodeError) as exc:
            logging.warning(
                "[%s] Could not read CancerVar cache, ignoring: %s", sample, exc,
            )
            cache = {}

    rows = list(clinical_rows)
    to_query = []
    for r in rows:
        key = _variant_key(
            r.get("Chr", ""), r.get("Start", ""),
            r.get("Ref", ""), r.get("Alt", ""),
        )
        if key in cache:
            continue
        try:
            chrom = _strip_chr(r.get("Chr", ""))
            pos = int(r.get("Start"))
            ref = str(r.get("Ref", ""))
            alt = str(r.get("Alt", ""))
        except (TypeError, ValueError):
            continue
        # CancerVar's byPosition endpoint only accepts SNVs and simple indels
        # in the standard chr/pos/ref/alt form -- there's no separate "end"
        # column. We pass what we have; the server returns an empty body for
        # variants it can't resolve, which we cache as a negative entry.
        to_query.append((key, chrom, pos, ref, alt))

    if not to_query:
        logging.info(
            "[%s] All %d clinical variants present in CancerVar cache.",
            sample, len(rows),
        )
        return cache

    logging.info(
        "[%s] Annotating %d new clinical variant(s) via CancerVar...",
        sample, len(to_query),
    )

    timestamp = datetime.utcnow().strftime("%Y-%m-%d %H:%M UTC")
    headers = {"Accept": "application/json"}
    n_hit = 0
    n_miss = 0
    n_err = 0

    for key, chrom, pos, ref, alt in to_query:
        params = {
            "queryType": "position",
            "build":     "hg38",
            "chr":       chrom,
            "pos":       pos,
            "ref":       ref,
            "alt":       alt,
        }
        try:
            response = requests.get(
                CANCERVAR_URL,
                params=params,
                headers=headers,
                timeout=TIMEOUT_S,
            )
        except requests.RequestException as exc:
            logging.warning(
                "[%s] CancerVar request failed for %s: %s",
                sample, key, exc,
            )
            n_err += 1
            time.sleep(SLEEP_BETWEEN_REQUESTS_S)
            continue

        if response.status_code != 200:
            logging.warning(
                "[%s] CancerVar HTTP %s for %s: %s",
                sample, response.status_code, key, response.text[:200],
            )
            n_err += 1
            time.sleep(SLEEP_BETWEEN_REQUESTS_S)
            continue

        # CancerVar returns an empty body for variants it can't classify.
        # We treat that as a negative cache entry so a future build doesn't
        # re-query them. Genuinely empty bodies have content length 0.
        body = response.text or ""
        if not body.strip():
            cache[key] = {"_found": False, "_fetched_at": timestamp}
            n_miss += 1
            time.sleep(SLEEP_BETWEEN_REQUESTS_S)
            continue

        try:
            data = response.json()
        except ValueError as exc:
            logging.warning(
                "[%s] CancerVar returned non-JSON for %s: %s",
                sample, key, exc,
            )
            n_err += 1
            time.sleep(SLEEP_BETWEEN_REQUESTS_S)
            continue

        if isinstance(data, dict) and "Cancervar" in data:
            ann = _summarise(data)
            ann["_fetched_at"] = timestamp
            cache[key] = ann
            n_hit += 1
        else:
            # Recognised JSON shape but no Cancervar field -- treat as miss.
            cache[key] = {"_found": False, "_fetched_at": timestamp}
            n_miss += 1

        time.sleep(SLEEP_BETWEEN_REQUESTS_S)

    logging.info(
        "[%s] CancerVar done: %d annotated, %d not classified, %d errors "
        "(cache total: %d).",
        sample, n_hit, n_miss, n_err, len(cache),
    )

    try:
        with open(cache_path, "w", encoding="utf-8") as fh:
            json.dump(cache, fh, indent=2, sort_keys=True)
        logging.info(
            "[%s] CancerVar cache written: %s entries -> %s",
            sample, len(cache), cache_path,
        )
    except OSError as exc:
        logging.warning(
            "[%s] Could not write CancerVar cache: %s", sample, exc,
        )

    return cache


def filter_for_frontend(cache):
    """Return only positive cache entries (``_found: True``) for template embed.

    The negative cache entries are useful on disk -- they prevent re-querying
    the same unresolvable variants on every build -- but they have nothing for
    the frontend to render, so we strip them out before embedding the JSON in
    the sample report HTML.
    """
    return {k: v for k, v in (cache or {}).items() if v.get("_found")}
