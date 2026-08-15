"""Parse per-exon coverage TSV.

Schema (from pipeline):
  Gene, Exon, Chr, Start, End, Length_bp, Mean_Coverage, Pct_100x, Pct_250x, Pct_500x, Flag

The pipeline's clinical convention is to INCLUDE duplicates in coverage metrics,
so this file already reflects that. No filtering is applied here.

Returns:
  {
    'columns': [list of column names],
    'rows':    [list of dicts],
    'n':       int,
    'summary': {
      'mean_of_per_exon_means': float,    # the headline "Mean coverage" value
      'n_exons':                int,
      'n_low_lt_100':           int,       # exons with mean < 100x
      'n_low_lt_250':           int,
      'low_lt_100_examples':    [{Gene, Exon, Mean_Coverage}, ... up to 8]
    } or None
  }
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
    cols = list(df.columns)
    rows = df.to_dict(orient="records")

    summary = None
    if "Mean_Coverage" in df.columns and len(df):
        # Coerce numerically. Pandas read everything as str (keep_default_na=False);
        # convert and treat blanks/non-numeric as NaN.
        numeric = pd.to_numeric(df["Mean_Coverage"], errors="coerce")
        n_exons = int(numeric.notna().sum())
        if n_exons > 0:
            mean_of_means = float(numeric.dropna().mean())
            n_lt_100 = int((numeric < 100).sum())
            n_lt_250 = int((numeric < 250).sum())

            # Sample of the worst-coverage exons (Gene/Exon/Mean), sorted ascending.
            low_df = (
                df.assign(_cov=numeric)
                  .loc[numeric < 100]
                  .sort_values("_cov", na_position="last")
                  .head(8)
            )
            low_examples = [
                {
                    "Gene": r.get("Gene", ""),
                    "Exon": r.get("Exon", ""),
                    "Mean_Coverage": r.get("Mean_Coverage", ""),
                }
                for r in low_df.to_dict(orient="records")
            ]

            summary = {
                "mean_of_per_exon_means": mean_of_means,
                "n_exons": n_exons,
                "n_low_lt_100": n_lt_100,
                "n_low_lt_250": n_lt_250,
                "low_lt_100_examples": low_examples,
            }

    return {"columns": cols, "rows": rows, "n": len(df), "summary": summary}
