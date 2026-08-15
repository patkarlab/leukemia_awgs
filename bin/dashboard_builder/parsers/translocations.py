"""
translocations.py - dashboard parser for merged, annotated structural variants.

Reads the per-sample annotated translocation table produced by
merge_translocations.py + annotate_mm_translocations.py:

    <effective_dir>/translocations/<sample>.mm_annotated.tsv

The schema is read from the file's own header row; nothing about gene names,
gene pairs, expected karyotypes or breakpoint coordinates is encoded here. The
only column knowledge used is which column *names* should be treated as numeric
for sorting, and which two coordinate pairs describe the breakpoints, so that
the template can build the paired IGV view.

Return contract
---------------
parse() returns either None (no table found) or a dict:

    {
      "path":       str,          # path actually read, for the Files tab
      "n_events":   int,
      "n_pass":     int,          # rows whose filter column is PASS
      "columns":    [ {"key","label","numeric","default_hidden"} , ... ],
      "rows":       [ {"cells": [...], "event_id", "igv": {...}}, ... ],
      "has_igv":    bool,         # whether breakpoint coordinates were resolvable
    }

Each cell is {"text": str, "order": float|None}. ``order`` is emitted into the
DataTables ``data-order`` attribute so that numeric and range-valued columns
sort correctly regardless of how they are formatted for display. Blank cells
sort below every populated value rather than alphabetically among them.

Standard library only. The dashboard environment has pandas, but a flat TSV of
a few dozen rows does not need it, and csv avoids the dtype coercion that turns
an empty support column into NaN and a sample identifier into a float.
"""

import csv
import json
import os
import re


# Columns rendered as numbers. Anything not listed sorts as text.
NUMERIC_COLUMNS = {
    "pos_a",
    "pos_b",
    "n_callers",
    "support_reads",
    "support_sniffles",
    "support_cutesv",
    "support_severus",
}

# Columns present in the table but hidden by default in the rendered view.
# They stay available through the DataTables column-visibility button.
DEFAULT_HIDDEN = {
    "sample",
    "supp_vec",
    "gene_a_source",
    "gene_b_source",
}

# Human-facing column labels. Any column absent from this map is displayed
# with its raw header name, so a schema change degrades to an ugly header
# rather than a missing column.
COLUMN_LABELS = {
    "sample": "Sample",
    "sv_id": "SV ID",
    "sv_type": "Type",
    "filter": "Filter",
    "chrom_a": "Chr A",
    "pos_a": "Pos A",
    "gene_a": "Gene / band A",
    "chrom_b": "Chr B",
    "pos_b": "Pos B",
    "gene_b": "Gene / band B",
    "gene_a_source": "A source",
    "gene_b_source": "B source",
    "known_mm_pair": "Known pair",
    "known_freq": "Reported freq",
    "callers": "Callers",
    "n_callers": "N callers",
    "supp_vec": "Support vector",
    "support_reads": "Supporting reads",
    "support_sniffles": "Sniffles",
    "support_cutesv": "CuteSV",
    "support_severus": "Severus",
}

# Column used as the default sort key, descending, when present.
DEFAULT_SORT_COLUMN = "support_reads"

# sv_type values that denote a rearrangement between two loci. Everything else
# in the annotated table (DEL, INS, INV, DUP) is a different evidence class and
# is hidden from this tab by default.
TRANSLOCATION_TYPES = {"TRA", "BND", "CTX"}


def _find_table(effective_dir, sample):
    """Locate the annotated translocation TSV for one sample.

    Looks first at the canonical bundle location, then falls back to a
    recursive search so the parser also works when pointed straight at a
    results tree rather than a report bundle.
    """
    direct = os.path.join(
        effective_dir, "translocations", "%s.mm_annotated.tsv" % sample
    )
    if os.path.isfile(direct):
        return direct

    # Bound to the sample identifier. There is deliberately no "any annotated
    # table" fallback: a directory that does not hold this sample's file must
    # render empty rather than silently show another sample's rearrangements.
    for root, _dirs, files in os.walk(effective_dir):
        for name in files:
            if name == "%s.mm_annotated.tsv" % sample:
                return os.path.join(root, name)
    return None


def _load_igv_manifest(effective_dir, sample):
    """Load the per-breakpoint IGV manifest written by igv_snapshots.py.

    The manifest maps each event identifier to the two standalone HTML pages
    holding its breakpoints. Returns an empty dict when IGV snapshots have not
    been generated, in which case the table renders without IGV buttons rather
    than with buttons that lead nowhere.
    """
    candidates = [
        os.path.join(
            effective_dir,
            "igv",
            "translocations",
            "%s.translocations.manifest.json" % sample,
        ),
        os.path.join(
            effective_dir, "igv", "%s.translocations.manifest.json" % sample
        ),
    ]

    manifest_path = None
    for candidate in candidates:
        if os.path.isfile(candidate):
            manifest_path = candidate
            break

    if manifest_path is None:
        for root, _dirs, files in os.walk(effective_dir):
            for name in files:
                if name == "%s.translocations.manifest.json" % sample:
                    manifest_path = os.path.join(root, name)
                    break
            if manifest_path:
                break

    if not manifest_path:
        return {}

    try:
        with open(manifest_path) as handle:
            manifest = json.load(handle)
    except (OSError, ValueError):
        return {}

    base = os.path.dirname(manifest_path)
    lookup = {}
    for event in manifest.get("events", []):
        event_id = event.get("event_id")
        if not event_id:
            continue
        entry = {}
        for side in ("a", "b"):
            html = (event.get(side) or {}).get("html")
            if not html:
                continue
            absolute = html if os.path.isabs(html) else os.path.join(base, html)
            if not os.path.isfile(absolute):
                continue
            try:
                entry[side] = os.path.relpath(absolute, effective_dir)
            except ValueError:
                entry[side] = absolute
        if entry:
            lookup[event_id] = entry
    return lookup


def _numeric_order(value):
    """Return a float sort key for a plain numeric cell, or None if blank."""
    text = (value or "").strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _frequency_order(value):
    """Return a float sort key for a reported-frequency string.

    Frequencies in the annotation dictionary are written as human ranges
    ("5-10%"), open bounds ("<1%") or approximations ("~15%"). Sorting those
    as text puts "<1%" above "5-10%". This extracts every number in the string
    and returns their mean, halving the value when the string is an upper
    bound, so the column orders by magnitude. Returns None when no number is
    present, which sorts the cell below all populated values.
    """
    text = (value or "").strip()
    if not text:
        return None
    numbers = [float(n) for n in re.findall(r"\d+(?:\.\d+)?", text)]
    if not numbers:
        return None
    mean = sum(numbers) / len(numbers)
    if text.lstrip().startswith("<"):
        return mean / 2.0
    return mean


def _order_for(key, value):
    """Dispatch a cell to the right sort-key function."""
    if key in NUMERIC_COLUMNS:
        return _numeric_order(value)
    if key == "known_freq":
        return _frequency_order(value)
    return None


def _locus(chrom, pos, flank):
    """Build an IGV locus string from a breakpoint, or None if incomplete."""
    chrom = (chrom or "").strip()
    pos = (pos or "").strip()
    if not chrom or not pos:
        return None
    try:
        centre = int(float(pos))
    except ValueError:
        return None
    start = max(1, centre - flank)
    end = centre + flank
    return "%s:%d-%d" % (chrom, start, end)


def parse(effective_dir, sample, igv_flank=5000):
    """Parse the annotated translocation table for one sample.

    ``igv_flank`` sets the half-width of the locus strings handed to the
    template for the paired IGV view. It should match the flanking value the
    IGV reports were generated with so that the requested locus falls inside
    the region the report actually embedded.
    """
    def not_found(reason):
        return {
            "found": False,
            "reason": reason,
            "searched": str(effective_dir),
            "n_events": 0,
            "n_translocations": 0,
            "n_other": 0,
            "n_pass": 0,
            "columns": [],
            "rows": [],
            "has_igv": False,
            "default_sort_index": None,
        }

    path = _find_table(effective_dir, sample)
    if not path:
        return not_found("no %s.mm_annotated.tsv under this directory" % sample)

    with open(path, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            return not_found("table is empty: %s" % path)
        header = [h.strip() for h in header]
        raw_rows = [row for row in reader if row and any(c.strip() for c in row)]

    columns = [
        {
            "key": key,
            "label": COLUMN_LABELS.get(key, key),
            "numeric": key in NUMERIC_COLUMNS,
            "default_hidden": key in DEFAULT_HIDDEN,
        }
        for key in header
    ]

    index = {key: i for i, key in enumerate(header)}
    igv_lookup = _load_igv_manifest(effective_dir, sample)

    def field(row, key):
        i = index.get(key)
        if i is None or i >= len(row):
            return ""
        return (row[i] or "").strip()

    rows = []
    n_pass = 0
    n_translocations = 0
    has_igv = False

    for raw in raw_rows:
        # Pad short rows so a truncated trailing column does not shift cells.
        padded = list(raw) + [""] * (len(header) - len(raw))

        cells = []
        for i, key in enumerate(header):
            value = (padded[i] or "").strip()
            cells.append({"text": value, "order": _order_for(key, value)})

        if field(padded, "filter").upper() == "PASS":
            n_pass += 1

        locus_a = _locus(
            field(padded, "chrom_a"), field(padded, "pos_a"), igv_flank
        )
        locus_b = _locus(
            field(padded, "chrom_b"), field(padded, "pos_b"), igv_flank
        )

        event_id = field(padded, "sv_id") or "event_%d" % (len(rows) + 1)
        pages = igv_lookup.get(event_id, {})
        if pages.get("a") or pages.get("b"):
            has_igv = True

        # The annotated table is annotation over the whole merged callset, so
        # it carries deletions, insertions, inversions and duplications
        # alongside the rearrangements. A row counts as a rearrangement if it
        # is typed TRA or if its two breakpoints sit on different contigs;
        # either condition alone would miss cases depending on how a caller
        # typed the record.
        sv_type = field(padded, "sv_type").upper()
        chrom_a = field(padded, "chrom_a")
        chrom_b = field(padded, "chrom_b")
        is_translocation = sv_type in TRANSLOCATION_TYPES or (
            bool(chrom_a) and bool(chrom_b) and chrom_a != chrom_b
        )
        if is_translocation:
            n_translocations += 1

        rows.append(
            {
                "cells": cells,
                "event_id": event_id,
                "sv_type": sv_type,
                "is_translocation": is_translocation,
                "igv": {
                    "locus_a": locus_a,
                    "locus_b": locus_b,
                    "html_a": pages.get("a"),
                    "html_b": pages.get("b"),
                    "label_a": field(padded, "gene_a")
                    or field(padded, "chrom_a"),
                    "label_b": field(padded, "gene_b")
                    or field(padded, "chrom_b"),
                    "point_a": "%s:%s"
                    % (field(padded, "chrom_a"), field(padded, "pos_a")),
                    "point_b": "%s:%s"
                    % (field(padded, "chrom_b"), field(padded, "pos_b")),
                },
            }
        )

    default_sort_index = None
    if DEFAULT_SORT_COLUMN in index:
        default_sort_index = index[DEFAULT_SORT_COLUMN]

    return {
        "found": True,
        "searched": str(effective_dir),
        "path": path,
        "filename": os.path.basename(path),
        "n_events": len(rows),
        "n_translocations": n_translocations,
        "n_other": len(rows) - n_translocations,
        "n_pass": n_pass,
        "columns": columns,
        "rows": rows,
        "has_igv": has_igv,
        "default_sort_index": default_sort_index,
        "igv_flank": igv_flank,
    }
