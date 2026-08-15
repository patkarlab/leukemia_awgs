"""
qc.py - dashboard parser for adaptive-sampling QC.

Replaces the hybrid-capture QC parsers inherited from the reference pipeline
(hsmetrics, fastp, Picard), none of which apply to nanopore adaptive sampling.
The QC that does apply is what QC_ONTARGET already emits per sample:

    <effective_dir>/qc/<sample>.region_coverage.tsv
    <effective_dir>/qc/<sample>.region_coverage.png
    <effective_dir>/qc/<sample>.readlen_hist.png
    <effective_dir>/qc/<sample>.qscore_hist.png
    <effective_dir>/qc/<sample>.readlen_qscore.tsv

The three PNGs are inlined as base64 so the report survives being moved or
mailed. The region coverage TSV is parsed into a sortable table with its own
header preserved, so a panel change that adds or renames a column shows up in
the report without a code change.

readlen_qscore.tsv is per-read and can run to millions of rows, so it is never
rendered as a table. It is listed as a downloadable artefact, and summary
statistics are computed from it in a single streaming pass when its columns can
be identified.

Standard library only.
"""

import base64
import csv
import os


# Per-read tables above this many rows are summarised, never rendered.
MAX_TABLE_ROWS = 2000

# Candidate column names for the per-read summary. Matching is case-insensitive
# and substring-based so minor header changes upstream do not break the summary.
READLEN_KEYS = ("read_length", "readlen", "length")
QSCORE_KEYS = ("mean_q", "qscore", "meanq", "quality")

# Column names in the region coverage table that should sort numerically.
NUMERIC_HINTS = (
    "start",
    "end",
    "depth",
    "coverage",
    "cov",
    "bases",
    "length",
    "size",
    "reads",
    "pct",
    "percent",
    "fraction",
    "mean",
    "median",
)

PLOTS = [
    (
        "region_coverage",
        ".region_coverage.png",
        "Per-region on-target coverage",
        "Mean depth per panel region. Wide breakpoint windows at the "
        "immunoglobulin and WWOX/MAF loci read lower than focal gene-body "
        "windows because depth is averaged over a much larger interval; "
        "interpret each region against its own window size, not against the "
        "panel-wide median.",
    ),
    (
        "readlen_hist",
        ".readlen_hist.png",
        "Read-length distribution",
        "On-target read lengths. The shearing target sits in the 8-12 kb "
        "window; a left-shifted mode indicates over-fragmentation.",
    ),
    (
        "qscore_hist",
        ".qscore_hist.png",
        "Per-read mean quality",
        "Distribution of per-read mean Q.",
    ),
]


def _inline_png(path):
    """Return a base64 data URI for a PNG, or None if unreadable."""
    try:
        with open(path, "rb") as handle:
            encoded = base64.b64encode(handle.read()).decode("ascii")
    except OSError:
        return None
    return "data:image/png;base64," + encoded


def _qc_dir(effective_dir, sample):
    """Locate the QC directory for this sample."""
    direct = os.path.join(effective_dir, "qc")
    if os.path.isdir(direct):
        return direct

    # Only accept a directory that holds THIS sample's coverage table. A
    # directory named qc/ belonging to another sample must not be picked up.
    target = "%s.region_coverage.tsv" % sample
    for root, dirs, _files in os.walk(effective_dir):
        for name in dirs:
            candidate = os.path.join(root, name)
            try:
                contents = os.listdir(candidate)
            except OSError:
                continue
            if target in contents:
                return candidate
    return None


def _is_numeric_column(name, values):
    """Decide whether a column sorts numerically.

    Uses the header name as a hint and then confirms against the data, so a
    column called "mean_depth" holding "NA" everywhere is not forced numeric.
    """
    lowered = name.lower()
    hinted = any(hint in lowered for hint in NUMERIC_HINTS)
    if not hinted:
        return False
    parsed = 0
    for value in values:
        text = (value or "").strip()
        if not text:
            continue
        try:
            float(text)
            parsed += 1
        except ValueError:
            return False
    return parsed > 0


def _read_table(path, max_rows=MAX_TABLE_ROWS):
    """Read a TSV into a header/rows structure with per-column numeric flags."""
    with open(path, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            return None
        header = [h.strip() for h in header]
        rows = []
        truncated = False
        for raw in reader:
            if not raw or not any(c.strip() for c in raw):
                continue
            if len(rows) >= max_rows:
                truncated = True
                break
            rows.append(list(raw) + [""] * (len(header) - len(raw)))

    columns = []
    for i, name in enumerate(header):
        values = [row[i] for row in rows]
        columns.append(
            {
                "key": name,
                "label": name.replace("_", " "),
                "numeric": _is_numeric_column(name, values),
            }
        )

    cell_rows = []
    for row in rows:
        cells = []
        for i, column in enumerate(columns):
            value = (row[i] or "").strip()
            order = None
            if column["numeric"] and value:
                try:
                    order = float(value)
                except ValueError:
                    order = None
            cells.append({"text": value, "order": order})
        cell_rows.append(cells)

    return {
        "columns": columns,
        "rows": cell_rows,
        "n_rows": len(cell_rows),
        "truncated": truncated,
    }


def _summarise_reads(path):
    """Compute read-count, N50 and mean Q from the per-read table.

    Streams the file so a multi-million-row table costs memory proportional to
    the read-length column alone. Returns None when the columns cannot be
    identified, in which case the file is still offered as a download.
    """
    try:
        handle = open(path, newline="")
    except OSError:
        return None

    with handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = [h.strip().lower() for h in next(reader)]
        except StopIteration:
            return None

        def find(keys):
            for i, name in enumerate(header):
                if any(key in name for key in keys):
                    return i
            return None

        len_i = find(READLEN_KEYS)
        q_i = find(QSCORE_KEYS)

        data_rows = []
        for row in reader:
            if row and any(c.strip() for c in row):
                data_rows.append(row)
            if len(data_rows) > 2:
                break

        # A single data row means this file is already a summary, not a
        # per-read table. Computing an N50 from it would report one read.
        if len(data_rows) == 1:
            values = list(data_rows[0]) + [""] * (len(header) - len(data_rows[0]))
            return {
                "mode": "summary",
                "pairs": [
                    {"key": k, "value": (values[i] or "").strip()}
                    for i, k in enumerate(header)
                    if k.strip()
                ],
            }

        if len_i is None:
            return None

        lengths = []
        q_total = 0.0
        q_count = 0
        for row in data_rows + list(reader):
            if len(row) <= len_i:
                continue
            try:
                lengths.append(int(float(row[len_i])))
            except (ValueError, TypeError):
                continue
            if q_i is not None and len(row) > q_i:
                try:
                    q_total += float(row[q_i])
                    q_count += 1
                except (ValueError, TypeError):
                    pass

    if not lengths:
        return None

    lengths.sort(reverse=True)
    total = sum(lengths)
    running = 0
    n50 = lengths[-1]
    for value in lengths:
        running += value
        if running >= total / 2.0:
            n50 = value
            break

    return {
        "mode": "computed",
        "n_reads": len(lengths),
        "total_bases": total,
        "read_n50": n50,
        "read_mean": total / float(len(lengths)),
        "read_max": lengths[0],
        "mean_q": (q_total / q_count) if q_count else None,
    }


DEPTH_HINTS = ("mean_depth", "median_depth", "depth", "mean_cov", "coverage")


def _depth_summary(coverage):
    """Median, mean and range of per-region depth.

    The panel deliberately mixes focal gene-body windows with megabase
    breakpoint windows, so a panel-wide mean understates the focal regions.
    The median across regions is the more representative single number, and
    the range is carried alongside so neither is read as uniform coverage.
    """
    if not coverage or not coverage.get("columns"):
        return None

    index = None
    for hint in DEPTH_HINTS:
        for i, column in enumerate(coverage["columns"]):
            if column["key"].lower() == hint:
                index = i
                break
        if index is not None:
            break
    if index is None:
        for i, column in enumerate(coverage["columns"]):
            if column["numeric"] and any(
                h in column["key"].lower() for h in ("depth", "cov")
            ):
                index = i
                break
    if index is None:
        return None

    values = []
    for row in coverage["rows"]:
        if index >= len(row):
            continue
        order = row[index].get("order")
        if order is not None:
            values.append(order)
    if not values:
        return None

    values.sort()
    n = len(values)
    median = values[n // 2] if n % 2 else (values[n // 2 - 1] + values[n // 2]) / 2.0

    return {
        "column": coverage["columns"][index]["key"],
        "n_regions": n,
        "median": median,
        "mean": sum(values) / float(n),
        "min": values[0],
        "max": values[-1],
        "n_below_5": sum(1 for v in values if v < 5),
        "n_below_10": sum(1 for v in values if v < 10),
    }


def _normalise_reads(reads):
    """Expose n_reads, read_n50 and mean_q whichever shape the file had.

    A one-row summary file and a per-read table both end up with the same
    three keys, so the templates do not have to know which was on disk.
    """
    if not reads:
        return reads
    if reads.get("mode") != "summary":
        return reads

    lookup = {}
    for pair in reads.get("pairs", []):
        lookup[pair["key"].strip().lower()] = pair["value"]

    def find(*candidates):
        for candidate in candidates:
            if candidate in lookup:
                return lookup[candidate]
        for key, value in lookup.items():
            if any(c in key for c in candidates):
                return value
        return None

    reads["n_reads"] = find("n_reads", "num_reads", "read_count", "reads")
    reads["read_n50"] = find("read_n50", "n50")
    reads["mean_q"] = find("mean_q", "meanq", "qscore")
    return reads


def parse(effective_dir, sample):
    """Collect QC plots, the coverage table and per-read summary statistics."""
    def not_found(reason):
        return {
            "found": False,
            "reason": reason,
            "searched": str(effective_dir),
            "plots": [],
            "coverage": None,
            "reads": None,
            "files": [],
        }

    qc_dir = _qc_dir(effective_dir, sample)
    if not qc_dir:
        return not_found("no qc/ directory under this sample")

    plots = []
    for key, suffix, title, caption in PLOTS:
        path = os.path.join(qc_dir, "%s%s" % (sample, suffix))
        if not os.path.isfile(path):
            continue
        data_uri = _inline_png(path)
        if not data_uri:
            continue
        plots.append(
            {
                "key": key,
                "title": title,
                "caption": caption,
                "src": data_uri,
                "filename": os.path.basename(path),
            }
        )

    coverage = None
    coverage_path = os.path.join(qc_dir, "%s.region_coverage.tsv" % sample)
    if os.path.isfile(coverage_path):
        coverage = _read_table(coverage_path)
        if coverage:
            coverage["filename"] = os.path.basename(coverage_path)

    reads = None
    reads_path = os.path.join(qc_dir, "%s.readlen_qscore.tsv" % sample)
    if os.path.isfile(reads_path):
        reads = _summarise_reads(reads_path)

    files = sorted(
        name
        for name in os.listdir(qc_dir)
        if os.path.isfile(os.path.join(qc_dir, name))
    )

    if not plots and not coverage and not reads:
        return not_found("qc/ exists but holds no recognised products: %s" % qc_dir)

    reads = _normalise_reads(reads)

    return {
        "found": True,
        "searched": str(effective_dir),
        "depth": _depth_summary(coverage),
        "qc_dir": qc_dir,
        "plots": plots,
        "coverage": coverage,
        "reads": reads,
        "files": files,
    }
