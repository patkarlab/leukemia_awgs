"""Optional build-time variant annotation via the GeneBe REST API.

GeneBe (https://genebe.net) provides ACMG classification, ClinVar status, gnomAD
allele frequencies and more. We POST batches of clinical variants and embed the
returned annotations into the per-sample dashboard.

Endpoint:
    POST https://api.genebe.net/cloud/api-public/v1/variants?genome=hg38
    Body: [{"chr":"2","pos":25234373,"ref":"C","alt":"T"}, ...]
    Optional Basic auth: -u email:api_key (higher rate limits)

This module is OPT-IN. The builder only calls it when --annotate-genebe is passed.

Annotations are cached on disk at:
    <sample_dir>/<sample>_genebe_cache.json

so that re-running the builder does not re-hit the network for variants we've
already seen. Cache entries are keyed by chr:pos:ref:alt. The cache is loaded,
updated for any missing variants, and rewritten.

Returns a dict keyed by chr:pos:ref:alt mapping to {acmg_classification,
acmg_criteria, clinvar_classification, clinvar_disease, gnomad_exome_af,
gnomad_genome_af, _fetched_at}. Missing keys mean "no annotation available".
"""

from __future__ import annotations

import json
import logging
import time
from datetime import datetime
from pathlib import Path
from typing import Iterable, Optional


GENEBE_BATCH_URL = "https://api.genebe.net/cloud/api-public/v1/variants"
BATCH_SIZE = 100
TIMEOUT_S = 60


def _variant_key(chrom: str, pos, ref: str, alt: str) -> str:
    return f"{chrom}:{pos}:{ref}:{alt}"


def _strip_chr(chrom: str) -> str:
    return str(chrom or "").replace("chr", "", 1) if str(chrom or "").startswith("chr") else str(chrom or "")


def _summarise_variant(variant_record: dict) -> dict:
    """Reduce a GeneBe API response record to the fields we display.

    Field names below are taken from a live GeneBe response (May 2026) -- they
    match the keys returned by /v1/variants on hg38.
    """
    out = {}
    out["acmg_classification"]    = variant_record.get("acmg_classification")
    out["acmg_criteria"]          = variant_record.get("acmg_criteria")
    out["acmg_score"]             = variant_record.get("acmg_score")
    out["clinvar_classification"] = variant_record.get("clinvar_classification")
    out["clinvar_disease"]        = variant_record.get("clinvar_disease")
    out["clinvar_review_status"]  = variant_record.get("clinvar_review_status")
    out["gnomad_exome_af"]        = variant_record.get("gnomad_exomes_af")
    out["gnomad_genome_af"]       = variant_record.get("gnomad_genomes_af")
    out["revel_score"]            = variant_record.get("revel_score")
    out["alphamissense_prediction"] = variant_record.get("alphamissense_prediction")
    out["spliceai_max_score"]     = variant_record.get("spliceai_max_score")
    out["effect"]                 = variant_record.get("effect")
    out["gene_symbol"]            = variant_record.get("gene_symbol")
    return out


def annotate(
    clinical_rows: Iterable[dict],
    sample_dir: Path,
    sample: str,
    api_user: Optional[str] = None,
    api_key: Optional[str] = None,
) -> dict:
    """Annotate the clinical-variants set via the GeneBe batch endpoint.

    Returns {chr:pos:ref:alt -> annotation dict}. On network/auth/parse failure,
    logs a warning and returns whatever has been cached so far.
    """
    # Lazy import so the builder is importable without these libs when
    # --annotate-genebe is not used.
    try:
        import requests
    except ImportError:
        logging.warning("[%s] GeneBe annotation requested but 'requests' is not installed.", sample)
        return {}

    cache_path = Path(sample_dir) / f"{sample}_genebe_cache.json"
    cache: dict = {}
    if cache_path.exists():
        try:
            with open(cache_path, "r", encoding="utf-8") as fh:
                cache = json.load(fh)
        except (OSError, json.JSONDecodeError) as exc:
            logging.warning("[%s] Could not read GeneBe cache, ignoring: %s", sample, exc)
            cache = {}

    rows = list(clinical_rows)
    to_query = []
    for r in rows:
        key = _variant_key(r.get("Chr", ""), r.get("Start", ""), r.get("Ref", ""), r.get("Alt", ""))
        if key not in cache:
            to_query.append((key, r))

    if not to_query:
        logging.info("[%s] All %d clinical variants present in GeneBe cache.", sample, len(rows))
        return cache

    logging.info("[%s] Annotating %d new clinical variants via GeneBe...", sample, len(to_query))

    auth = None
    if api_user and api_key:
        auth = (api_user, api_key)

    timestamp = datetime.utcnow().strftime("%Y-%m-%d %H:%M UTC")

    for i in range(0, len(to_query), BATCH_SIZE):
        chunk = to_query[i : i + BATCH_SIZE]
        payload = []
        for key, r in chunk:
            try:
                payload.append({
                    "chr":  _strip_chr(r.get("Chr", "")),
                    "pos":  int(r.get("Start")),
                    "ref":  str(r.get("Ref", "")),
                    "alt":  str(r.get("Alt", "")),
                })
            except (TypeError, ValueError):
                logging.warning("[%s] Skipping malformed variant for GeneBe: %s", sample, key)
                continue

        if not payload:
            continue

        try:
            response = requests.post(
                GENEBE_BATCH_URL,
                params={"genome": "hg38"},
                json=payload,
                auth=auth,
                timeout=TIMEOUT_S,
                headers={"Accept": "application/json", "Content-Type": "application/json"},
            )
        except requests.RequestException as exc:
            logging.warning("[%s] GeneBe request failed (chunk %d): %s", sample, i // BATCH_SIZE, exc)
            continue

        if response.status_code != 200:
            logging.warning(
                "[%s] GeneBe returned HTTP %s on chunk %d: %s",
                sample, response.status_code, i // BATCH_SIZE, response.text[:200],
            )
            continue

        try:
            data = response.json()
        except ValueError as exc:
            logging.warning("[%s] GeneBe returned non-JSON: %s", sample, exc)
            continue

        variants_out = data.get("variants") if isinstance(data, dict) else data
        if not isinstance(variants_out, list):
            logging.warning("[%s] GeneBe response shape unexpected; skipping.", sample)
            continue

        # Each returned variant may have re-normalised chr/pos/ref/alt (e.g. left-trim).
        # We index input -> output by position within the chunk.
        for (key, _), out in zip(chunk, variants_out):
            ann = _summarise_variant(out)
            ann["_fetched_at"] = timestamp
            cache[key] = ann

        # Gentle pacing — unauthenticated requests are rate-limited.
        time.sleep(0.5)

    # Persist updated cache (best-effort).
    try:
        with open(cache_path, "w", encoding="utf-8") as fh:
            json.dump(cache, fh, indent=2, sort_keys=True)
        logging.info("[%s] GeneBe cache written: %s entries -> %s", sample, len(cache), cache_path)
    except OSError as exc:
        logging.warning("[%s] Could not write GeneBe cache: %s", sample, exc)

    return cache
