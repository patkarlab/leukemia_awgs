"""Parse Picard CollectHsMetrics output.

Picard files have two blocks:

  ## METRICS CLASS picard.analysis.directed.HsMetrics
  <header line>
  <data line>

  ## HISTOGRAM java.lang.Integer
  <header line>
  <rows>

We return a dict with:
  - metrics: dict of column -> value (numeric where possible, else str)
  - histogram: list of dicts (one per row) — only included if present

Returns None if the file cannot be parsed.
"""

from pathlib import Path


def _coerce(value):
    """Best-effort numeric coercion. Returns float, int, or original string."""
    if value == "" or value == "?":
        return None
    try:
        f = float(value)
        # If it round-trips to an int and there is no decimal, return int.
        if f.is_integer() and "." not in value and "e" not in value.lower():
            return int(f)
        return f
    except ValueError:
        return value


def parse(path):
    """Parse a Picard HsMetrics file.

    Parameters
    ----------
    path : pathlib.Path or str
        Path to the *_hsmetrics.txt file.

    Returns
    -------
    dict or None
        {'metrics': {...}, 'histogram': [{'depth': int, 'count': int}, ...]}
    """
    path = Path(path)
    if not path.exists():
        return None

    try:
        with open(path, "r") as fh:
            lines = fh.read().splitlines()
    except OSError:
        return None

    metrics = {}
    histogram = []

    # Find the METRICS CLASS block. The two lines after it are header + data.
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("## METRICS CLASS"):
            # Next non-empty line is the header, the line after is the data.
            j = i + 1
            while j < len(lines) and not lines[j].strip():
                j += 1
            if j + 1 < len(lines):
                header = lines[j].split("\t")
                data = lines[j + 1].split("\t")
                for col, val in zip(header, data):
                    metrics[col] = _coerce(val)
            i = j + 2
            continue

        if line.startswith("## HISTOGRAM"):
            j = i + 1
            while j < len(lines) and not lines[j].strip():
                j += 1
            if j >= len(lines):
                break
            hist_header = lines[j].split("\t")
            j += 1
            while j < len(lines):
                row = lines[j].strip()
                if not row or row.startswith("#"):
                    j += 1
                    continue
                parts = row.split("\t")
                if len(parts) < 2:
                    j += 1
                    continue
                # Picard column 1 = coverage_or_base_quality, column 2 = high_quality_coverage_count
                # Column 3 (unfiltered_baseq_count) is the base-quality histogram, kept aside.
                depth = _coerce(parts[0])
                count = _coerce(parts[1]) if len(parts) > 1 else None
                histogram.append({"depth": depth, "count": count})
                j += 1
            i = j
            continue
        i += 1

    if not metrics:
        return None

    return {"metrics": metrics, "histogram": histogram, "hist_header": hist_header if histogram else []}
