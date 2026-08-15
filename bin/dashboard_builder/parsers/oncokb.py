"""Optional build-time OncoKB annotation of clinical variants.

OncoKB (https://www.oncokb.org) is a precision-oncology knowledge base curated
by Memorial Sloan Kettering. The REST API requires an authentication token --
register a free academic account at https://www.oncokb.org/account/register
and copy the token from your profile page, then pass it via --oncokb-token.

Endpoint:
    GET https://www.oncokb.org/api/v1/annotate/mutations/byGenomicChange
        ?genomicLocation=<chr>,<start>,<end>,<ref>,<alt>&referenceGenome=GRCh38
    Headers: Authorization: Bearer <token>, Accept: application/json

Returns oncogenic classification (Oncogenic/Likely Oncogenic/Likely Neutral/
Inconclusive/Unknown), mutation effect, hotspot membership, the variant's
highest FDA/CLIA evidence level, and (when present) actionable treatments
with drug names + sensitivity/resistance level.

This module is OPT-IN. The builder only calls it when --annotate-oncokb is
passed AND --oncokb-token is non-empty.

Annotations are cached on disk at:
    <sample_dir>/<sample>_oncokb_cache.json

keyed by chr:pos:ref:alt (matching the GeneBe/MobiDetails cache key scheme).
"""

from pathlib import Path
import json
import logging
import time
from datetime import datetime


ONCOKB_URL = "https://www.oncokb.org/api/v1/annotate/mutations/byGenomicChange"
TIMEOUT_S = 30
SLEEP_BETWEEN_REQUESTS_S = 0.2  # be polite even with a valid token


def _variant_key(chrom, pos, ref, alt):
    return f"{chrom}:{pos}:{ref}:{alt}"


def _strip_chr(chrom):
    s = str(chrom or "")
    return s[3:] if s.startswith("chr") else s


def _summarise(record):
    """Reduce an OncoKB response to the fields we display in the dashboard.

    Field names are taken from a live OncoKB v1 response. We deliberately keep
    a small surface; the full record is verbose and contains nested objects.
    """
    out = {}
    out["oncogenic"]                  = record.get("oncogenic")
    out["mutation_effect_known"]      = (record.get("mutationEffect") or {}).get("knownEffect")
    out["mutation_effect_description"] = (record.get("mutationEffect") or {}).get("description")
    out["highest_sensitive_level"]    = record.get("highestSensitiveLevel")
    out["highest_resistance_level"]   = record.get("highestResistanceLevel")
    out["highest_diagnostic_implication_level"] = record.get("highestDiagnosticImplicationLevel")
    out["highest_prognostic_implication_level"] = record.get("highestPrognosticImplicationLevel")
    out["hotspot"]                    = record.get("hotspot")
    out["vus"]                        = record.get("vus")
    # treatments: list of {drugs:[{drugName}], level, levelAssociatedCancerType, pmids, abstracts}
    treatments_out = []
    for t in (record.get("treatments") or []):
        treatments_out.append({
            "drugs":   [d.get("drugName") for d in (t.get("drugs") or []) if d.get("drugName")],
            "level":   t.get("level"),
            "cancer_type": (t.get("levelAssociatedCancerType") or {}).get("mainType"),
            "pmids":   t.get("pmids", []),
        })
    out["treatments"] = treatments_out
    # geneSummary / variantSummary are short prose blurbs OncoKB writes; they
    # are well-suited to embed in a tooltip / expand panel.
    out["gene_summary"]               = record.get("geneSummary")
    out["variant_summary"]            = record.get("variantSummary")
    out["tumor_type_summary"]         = record.get("tumorTypeSummary")
    return out


def annotate(clinical_rows, sample_dir, sample, oncokb_token):
    """Annotate clinical variants via OncoKB byGenomicChange.

    Returns {chr:pos:ref:alt -> annotation dict}. On any failure mode (network,
    auth, parse) we log and return whatever we have cached so far. Per-variant
    request failures don't abort the whole annotation pass.
    """
    try:
        import requests
    except ImportError:
        logging.warning(
            "[%s] OncoKB annotation requested but 'requests' is not installed.",
            sample,
        )
        return {}

    if not oncokb_token:
        logging.warning(
            "[%s] OncoKB annotation requested but --oncokb-token was not provided. "
            "Register at https://www.oncokb.org/account/register for a free token.",
            sample,
        )
        return {}

    cache_path = Path(sample_dir) / f"{sample}_oncokb_cache.json"
    cache = {}
    if cache_path.exists():
        try:
            with open(cache_path, "r", encoding="utf-8") as fh:
                cache = json.load(fh)
        except (OSError, json.JSONDecodeError) as exc:
            logging.warning(
                "[%s] Could not read OncoKB cache, ignoring: %s", sample, exc,
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
        # OncoKB's byGenomicChange takes start AND end -- for a SNV they are
        # the same; for indels OncoKB normalises internally so passing
        # start = pos and end = pos + max(len(ref), len(alt)) - 1 is a safe
        # default that matches their docs.
        end = pos + max(len(ref), len(alt)) - 1
        to_query.append((key, chrom, pos, end, ref, alt))

    if not to_query:
        logging.info(
            "[%s] All %d clinical variants present in OncoKB cache.",
            sample, len(rows),
        )
        return cache

    logging.info(
        "[%s] Annotating %d new clinical variant(s) via OncoKB...",
        sample, len(to_query),
    )

    timestamp = datetime.utcnow().strftime("%Y-%m-%d %H:%M UTC")
    headers = {
        "Authorization": f"Bearer {oncokb_token}",
        "Accept": "application/json",
    }
    n_hit = 0
    n_err = 0

    for key, chrom, pos, end, ref, alt in to_query:
        params = {
            "genomicLocation": f"{chrom},{pos},{end},{ref},{alt}",
            "referenceGenome": "GRCh38",
        }
        try:
            response = requests.get(
                ONCOKB_URL,
                params=params,
                headers=headers,
                timeout=TIMEOUT_S,
            )
        except requests.RequestException as exc:
            logging.warning(
                "[%s] OncoKB request failed for %s: %s",
                sample, key, exc,
            )
            n_err += 1
            time.sleep(SLEEP_BETWEEN_REQUESTS_S)
            continue

        if response.status_code == 401:
            logging.warning(
                "[%s] OncoKB returned 401 Unauthorized. Check that the token "
                "is correct and your account is active. Aborting OncoKB pass.",
                sample,
            )
            break

        if response.status_code != 200:
            logging.warning(
                "[%s] OncoKB HTTP %s for %s: %s",
                sample, response.status_code, key, response.text[:200],
            )
            n_err += 1
            time.sleep(SLEEP_BETWEEN_REQUESTS_S)
            continue

        try:
            data = response.json()
        except ValueError as exc:
            logging.warning(
                "[%s] OncoKB returned non-JSON for %s: %s",
                sample, key, exc,
            )
            n_err += 1
            time.sleep(SLEEP_BETWEEN_REQUESTS_S)
            continue

        if isinstance(data, dict):
            ann = _summarise(data)
            ann["_fetched_at"] = timestamp
            cache[key] = ann
            n_hit += 1
        else:
            logging.warning(
                "[%s] OncoKB response shape unexpected for %s.", sample, key,
            )
            n_err += 1

        time.sleep(SLEEP_BETWEEN_REQUESTS_S)

    logging.info(
        "[%s] OncoKB done: %d annotated, %d errors (cache total: %d).",
        sample, n_hit, n_err, len(cache),
    )

    try:
        with open(cache_path, "w", encoding="utf-8") as fh:
            json.dump(cache, fh, indent=2, sort_keys=True)
        logging.info(
            "[%s] OncoKB cache written: %s entries -> %s",
            sample, len(cache), cache_path,
        )
    except OSError as exc:
        logging.warning(
            "[%s] Could not write OncoKB cache: %s", sample, exc,
        )

    return cache
