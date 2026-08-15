"""Optional build-time variant annotation via the MobiDetails REST API.

MobiDetails (https://mobidetails.chu-montpellier.fr) hosts curated genomic
variant records with HGVS, transcript, splicing and population-frequency
data. We use the `/api/variant/exists/{variant_ghgvs}` endpoint to resolve
an HGVS genomic string to a stable MobiDetails record URL at build time, so
the dashboard can deep-link directly to the variant page (skipping the
create/annotate step the live "Open MobiDetails" button triggers).

Endpoint (no API key required, public access for lookup only):

    GET https://mobidetails.chu-montpellier.fr/api/variant/exists/{variant_ghgvs}
    -> 200 with {"mobidetails_id": int, "url": str}            (found)
    -> 200 with {"mobidetails_warning": "...does not exist..."} (not found)

This module is OPT-IN. The builder only calls it when --annotate-mobidetails
is passed. Annotations are cached on disk at:

    <sample_dir>/<sample>_mobidetails_cache.json

keyed by chr:pos:ref:alt (matching the GeneBe cache key scheme) so that
re-runs do not re-hit the network. Cache entries are
{mobidetails_id, url, hgvs_g, _fetched_at} on hit, or
{not_found: true, hgvs_g, _fetched_at} when MD has no record yet.

Variants with empty or "-1" VV_HGVSg are skipped silently — they cannot be
resolved by genomic HGVS lookup.
"""

from pathlib import Path
import json
import logging
import time
from datetime import datetime
from urllib.parse import quote


MD_EXISTS_URL = "https://mobidetails.chu-montpellier.fr/api/variant/exists/{ghgvs}"
TIMEOUT_S = 30
# Politeness pause between requests. The endpoint has no documented rate limit
# but a small delay keeps us a good citizen for unauthenticated lookups.
SLEEP_BETWEEN_REQUESTS_S = 0.2


def _variant_key(chrom, pos, ref, alt):
    return f"{chrom}:{pos}:{ref}:{alt}"


def _is_useful(value):
    """True iff a TSV cell carries real data (not blank, not the -1 sentinel)."""
    if value is None:
        return False
    s = str(value).strip()
    return s != "" and s != "-1" and s != "-1.0"


def annotate(clinical_rows, sample_dir, sample):
    """Resolve clinical variants to MobiDetails record URLs.

    Returns {chr:pos:ref:alt -> annotation dict}. Each annotation either has
    {mobidetails_id, url, hgvs_g, _fetched_at} (success) or
    {not_found: True, hgvs_g, _fetched_at} (MD has no record yet) or is omitted
    entirely (HGVSg missing or request failed). On any failure mode we log and
    return what we have so far -- the dashboard JS falls back gracefully to
    the create_g deep-link for any unresolved variant.
    """
    try:
        import requests
    except ImportError:
        logging.warning(
            "[%s] MobiDetails annotation requested but 'requests' is not installed.",
            sample,
        )
        return {}

    cache_path = Path(sample_dir) / f"{sample}_mobidetails_cache.json"
    cache = {}
    if cache_path.exists():
        try:
            with open(cache_path, "r", encoding="utf-8") as fh:
                cache = json.load(fh)
        except (OSError, json.JSONDecodeError) as exc:
            logging.warning(
                "[%s] Could not read MobiDetails cache, ignoring: %s",
                sample, exc,
            )
            cache = {}

    rows = list(clinical_rows)
    to_query = []
    n_skipped_no_hgvs = 0
    for r in rows:
        key = _variant_key(
            r.get("Chr", ""), r.get("Start", ""),
            r.get("Ref", ""), r.get("Alt", ""),
        )
        if key in cache:
            continue
        hgvs_g = r.get("VV_HGVSg", "")
        if not _is_useful(hgvs_g):
            n_skipped_no_hgvs += 1
            continue
        to_query.append((key, str(hgvs_g)))

    if n_skipped_no_hgvs:
        logging.info(
            "[%s] %d clinical variant(s) had no VV_HGVSg -- skipped for MobiDetails.",
            sample, n_skipped_no_hgvs,
        )

    if not to_query:
        logging.info(
            "[%s] All %d clinical variants present in MobiDetails cache (or unresolvable).",
            sample, len(rows),
        )
        return cache

    logging.info(
        "[%s] Looking up %d clinical variant(s) in MobiDetails...",
        sample, len(to_query),
    )

    timestamp = datetime.utcnow().strftime("%Y-%m-%d %H:%M UTC")
    n_hit = 0
    n_miss = 0
    n_err = 0

    for key, hgvs_g in to_query:
        url = MD_EXISTS_URL.format(ghgvs=quote(hgvs_g, safe=""))
        try:
            response = requests.get(
                url,
                timeout=TIMEOUT_S,
                headers={"Accept": "application/json"},
            )
        except requests.RequestException as exc:
            logging.warning(
                "[%s] MobiDetails request failed for %s: %s",
                sample, hgvs_g, exc,
            )
            n_err += 1
            time.sleep(SLEEP_BETWEEN_REQUESTS_S)
            continue

        if response.status_code != 200:
            logging.warning(
                "[%s] MobiDetails returned HTTP %s for %s: %s",
                sample, response.status_code, hgvs_g, response.text[:200],
            )
            n_err += 1
            time.sleep(SLEEP_BETWEEN_REQUESTS_S)
            continue

        try:
            data = response.json()
        except ValueError as exc:
            logging.warning(
                "[%s] MobiDetails returned non-JSON for %s: %s",
                sample, hgvs_g, exc,
            )
            n_err += 1
            time.sleep(SLEEP_BETWEEN_REQUESTS_S)
            continue

        if isinstance(data, dict) and "mobidetails_id" in data:
            cache[key] = {
                "mobidetails_id": data.get("mobidetails_id"),
                "url":            data.get("url"),
                "hgvs_g":         hgvs_g,
                "_fetched_at":    timestamp,
            }
            n_hit += 1
        elif isinstance(data, dict) and "mobidetails_warning" in data:
            cache[key] = {
                "not_found":   True,
                "warning":     data.get("mobidetails_warning"),
                "hgvs_g":      hgvs_g,
                "_fetched_at": timestamp,
            }
            n_miss += 1
        else:
            logging.warning(
                "[%s] MobiDetails response shape unexpected for %s: %s",
                sample, hgvs_g, str(data)[:200],
            )
            n_err += 1

        time.sleep(SLEEP_BETWEEN_REQUESTS_S)

    logging.info(
        "[%s] MobiDetails done: %d found, %d not in MD, %d errors (cache total: %d).",
        sample, n_hit, n_miss, n_err, len(cache),
    )

    try:
        with open(cache_path, "w", encoding="utf-8") as fh:
            json.dump(cache, fh, indent=2, sort_keys=True)
        logging.info(
            "[%s] MobiDetails cache written: %s entries -> %s",
            sample, len(cache), cache_path,
        )
    except OSError as exc:
        logging.warning(
            "[%s] Could not write MobiDetails cache: %s",
            sample, exc,
        )

    return cache


def filter_for_frontend(cache):
    """Strip not-found entries from the cache before embedding in HTML.

    The dashboard JS only needs {url} for hits -- not-found entries waste
    bytes in the inline JSON and would have no chip to render anyway, since
    the JS already falls back to the create_g deep-link.
    """
    out = {}
    for key, ann in (cache or {}).items():
        if not isinstance(ann, dict):
            continue
        if ann.get("not_found"):
            continue
        url = ann.get("url")
        if not url:
            continue
        out[key] = {
            "url":            url,
            "mobidetails_id": ann.get("mobidetails_id"),
            "_fetched_at":    ann.get("_fetched_at"),
        }
    return out
