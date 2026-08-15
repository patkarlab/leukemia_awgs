"""Parse FLT3-ITD consensus output.

flt3_consensus.tsv columns:
  sample, status, n_tools, tools, length_bp, length_range, pos_hg38,
  vaf_pct_min, vaf_pct_max, vaf_pct_mean, ar_min, ar_max, ar_mean,
  hgvsc, hgvsp, domain, inserted_seq, raw_calls

Returns:
  {
    'rows': [list of dicts],
    'n_positive': int rows with status indicating ITD detected (any non-empty status that is not 'negative')
  }
or None if the file is missing/unreadable.

Note: when the file has only a header (no rows), the sample is FLT3-ITD negative.
"""

from pathlib import Path

import pandas as pd


def parse(path):
    path = Path(path)
    if not path.exists():
        return None

    try:
        df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, na_values=[""])
    except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError):
        return None

    df = df.fillna("")
    rows = df.to_dict(orient="records")

    n_positive = 0
    if "status" in df.columns:
        # Treat anything that is non-empty and not "negative"/"no_itd" as a positive event.
        statuses = df["status"].str.lower()
        n_positive = int((statuses.notna() & (statuses != "") & (~statuses.isin(["negative", "no_itd", "no-itd"]))).sum())

    return {"rows": rows, "n_positive": n_positive, "n_total": len(rows)}
