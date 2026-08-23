#!/usr/bin/env python3
"""
igv_snapshots.py

Generate self-contained igv-reports HTML snapshots for the two evidence
classes this pipeline reports on.

Modes
-----
somatic
    One report per sample covering every row of the clinical SNV table
    (<sample>.clinical.tsv, the output of the somatic filter). Rendered
    against the hg38 BAM, because that is the reference the somatic calls
    were made on. This is the report the Variants tab links into.

translocations
    Two reports per event, one for each breakpoint, rendered against the T2T
    BAM. Each is a standalone page, which is what allows the dashboard to show
    both partners of a rearrangement side by side in a single view: the two
    pages are loaded into two iframes. Nothing is shared between them, so
    there is no cross-frame scripting, no srcdoc escaping, and no dependency
    on igv-reports' internal row router.

    A manifest (<sample>.translocations.manifest.json) records which page
    belongs to which breakpoint of which event, so the dashboard resolves
    links by identifier rather than by reconstructing filenames.

Design notes
------------
- Standard library only. igv-reports' ``create_report`` must be on PATH; it is
  present in the awgs_sv environment.
- create_report is driven through a tab-delimited sites file rather than a
  VCF, so the variant table in each page shows exactly the columns written
  into that file, under their own names.
- No variant, gene pair, breakpoint coordinate, FISH finding or expected
  karyotype is encoded here. Every site comes from the upstream table. The
  only column knowledge used is which column names hold coordinates.
- Zero rows is a normal outcome, not an error: a sample may legitimately have
  no clinical SNVs or no annotated translocations. In that case a small
  placeholder page is written and the exit status is 0, so the report bundle
  and dashboard still compose cleanly.

Exit codes
----------
0   success, including the legitimate zero-sites case
2   usage or input error (missing input, missing BAM, create_report failure)
"""

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile


# Column names holding coordinates, per input type. Lookup is by name, so a
# reordered upstream table does not break the mapping.
SOMATIC_CHROM = "chrom"
SOMATIC_POS = "pos"

TRANSLOCATION_SIDES = (
    ("a", "chrom_a", "pos_a", "gene_a"),
    ("b", "chrom_b", "pos_b", "gene_b"),
)
TRANSLOCATION_ID = "sv_id"

# Columns carried into the somatic sites file, in display order, when present.
SOMATIC_INFO_COLUMNS = [
    "gene",
    "panel_label",
    "variant_type",
    "consequence",
    "csq_primary",
    "impact",
    "hgvsc",
    "hgvsp",
    "transcript",
    "exon_rank",
    "exon_total",
    "rs_id",
    "clinvar_sig",
    "pop_af_max",
    "tumor_af_pct",
    "REF_COUNT",
    "ALT_COUNT",
    "DP",
    "qual",
    "Filter",
]

# Columns carried into each translocation breakpoint sites file, when present.
TRANSLOCATION_INFO_COLUMNS = [
    "sv_id",
    "sv_type",
    "filter",
    "breakpoint",
    "partner_locus",
    "partner_gene",
    "gene_a",
    "gene_b",
    "known_mm_pair",
    "known_pair",      # the acute-leukaemia annotator's name for the same field
    "entity",
    "tier",
    "match_quality",
    "anchor",
    "known_freq",
    "callers",
    "n_callers",
    "support_reads",
    "support_sniffles",
    "support_cutesv",
    "support_severus",
]

PLACEHOLDER_HTML = """<!DOCTYPE html>
<html lang="en">
<head><meta charset="utf-8"><title>{title}</title>
<style>body{{font-family:system-ui,-apple-system,"Segoe UI",sans-serif;
margin:3rem;color:#333}}p{{max-width:44rem;line-height:1.5}}</style></head>
<body><h2>{title}</h2><p>{message}</p></body></html>
"""


def eprint(*args):
    """Write to stderr; stdout is reserved for machine-readable summaries."""
    print(*args, file=sys.stderr)


def fail(message, code=2):
    eprint("ERROR: " + message)
    sys.exit(code)


def safe_name(value):
    """Reduce an identifier to characters that are safe in a filename."""
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("._-")
    return cleaned or "event"


def read_tsv(path):
    """Read a tab-delimited file with a header into (header, rows-as-dicts)."""
    with open(path, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = [h.strip() for h in next(reader)]
        except StopIteration:
            return [], []
        rows = []
        for raw in reader:
            if not raw or not any(cell.strip() for cell in raw):
                continue
            padded = list(raw) + [""] * (len(header) - len(raw))
            rows.append(
                {key: (padded[i] or "").strip() for i, key in enumerate(header)}
            )
    return header, rows


def write_placeholder(out_html, title, message):
    """Write the no-sites page and report success."""
    os.makedirs(os.path.dirname(os.path.abspath(out_html)) or ".", exist_ok=True)
    with open(out_html, "w") as handle:
        handle.write(PLACEHOLDER_HTML.format(title=title, message=message))
    eprint("No sites to render; wrote placeholder: %s" % out_html)


def ensure_bam_index(bam):
    """Index the BAM if no index is present. Non-fatal on failure."""
    for suffix in (".bai", ".csi"):
        if os.path.isfile(bam + suffix):
            return
    stem = os.path.splitext(bam)[0]
    for suffix in (".bai", ".csi"):
        if os.path.isfile(stem + suffix):
            return
    eprint("BAM index missing, indexing: %s" % bam)
    result = subprocess.run(
        ["samtools", "index", "-@", "4", bam],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if result.returncode != 0:
        eprint(
            "WARNING: samtools index failed; create_report may fail:\n%s"
            % result.stderr.decode("utf-8", "replace")
        )


def run_create_report(sites_path, header, bam, fasta, out_html, flanking, title):
    """Invoke create_report for one sites file.

    The sites file always carries coordinates in the first three columns, so
    --sequence/--begin/--end are fixed at 1/2/3 and remain correct regardless
    of what informational columns follow.
    """
    os.makedirs(os.path.dirname(os.path.abspath(out_html)) or ".", exist_ok=True)

    command = [
        "create_report",
        sites_path,
        "--sequence", "1",
        "--begin", "2",
        "--end", "3",
        "--fasta", fasta,
        "--tracks", bam,
        "--flanking", str(flanking),
        "--output", out_html,
    ]
    if title:
        command += ["--title", title]

    result = subprocess.run(
        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    if result.returncode != 0:
        eprint("create_report failed. Command was:")
        eprint("  " + " ".join(command))
        eprint(result.stdout.decode("utf-8", "replace"))
        eprint(result.stderr.decode("utf-8", "replace"))
        return False
    return True


def write_sites(rows, columns, path):
    """Write a sites TSV: CHROM, START, END, then the informational columns.

    The coordinate columns are capitalised because igv-reports carries the
    sites file's header straight into the tableJson it embeds, and the
    dashboard's lookup parser reads that header by name. Renaming these is
    enough to break the cross-link between a variant card and its IGV row.
    """
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["CHROM", "START", "END"] + columns)
        for row in rows:
            writer.writerow(
                [row["chrom"], row["start"], row["end"]]
                + [row.get(column, "") for column in columns]
            )


def site_bounds(position):
    """Convert a 1-based point position to a half-open single-base interval."""
    point = int(float(position))
    return max(0, point - 1), point


# ---------------------------------------------------------------------------
# Mode: somatic
# ---------------------------------------------------------------------------

def run_somatic(args):
    header, rows = read_tsv(args.sites_tsv)
    title = "%s clinical SNVs" % args.sample

    if not rows:
        write_placeholder(
            args.out_html,
            title,
            "This sample has no clinical SNV or indel calls. An empty result "
            "is an expected outcome, not a pipeline failure.",
        )
        return 0

    if SOMATIC_CHROM not in header or SOMATIC_POS not in header:
        fail(
            "clinical table lacks '%s'/'%s' columns: %s"
            % (SOMATIC_CHROM, SOMATIC_POS, args.sites_tsv)
        )

    ensure_bam_index(args.bam)

    # POSITION, REF and ALT are written explicitly and first. The dashboard
    # builds its lookup key as CHROM:POSITION:REF:ALT and matches it against
    # the variant row's Chr:Start:Ref:Alt, so all four must be present in the
    # rendered table and POSITION must be the 1-based coordinate -- not START,
    # which is 0-based for the browser's benefit. Omitting REF and ALT, as the
    # previous column set did, leaves the lookup empty and every variant card
    # reports "no IGV" with nothing to say why.
    info_columns = ["POSITION", "REF", "ALT"] + [
        c for c in SOMATIC_INFO_COLUMNS if c in header
    ]
    sites = []
    skipped = 0
    for row in rows:
        try:
            start, end = site_bounds(row[SOMATIC_POS])
        except (ValueError, KeyError):
            skipped += 1
            continue
        site = dict(row)
        site["chrom"] = row[SOMATIC_CHROM]
        site["start"] = start
        site["end"] = end
        site["POSITION"] = row[SOMATIC_POS]
        site["REF"] = row.get("ref", "")
        site["ALT"] = row.get("alt", "")
        sites.append(site)

    if skipped:
        eprint("Skipped %d row(s) with unparseable positions" % skipped)

    if not sites:
        write_placeholder(
            args.out_html, title, "No rows carried a usable coordinate."
        )
        return 0

    workdir = tempfile.mkdtemp(prefix="igv_somatic_")
    try:
        sites_path = os.path.join(workdir, "sites.tsv")
        write_sites(sites, info_columns, sites_path)
        ok = run_create_report(
            sites_path,
            info_columns,
            args.bam,
            args.fasta,
            args.out_html,
            args.flanking,
            title,
        )
    finally:
        shutil.rmtree(workdir, ignore_errors=True)

    if not ok:
        return 2
    print("somatic\t%s\t%d\t%s" % (args.sample, len(sites), args.out_html))
    return 0


# ---------------------------------------------------------------------------
# Mode: translocations
# ---------------------------------------------------------------------------

def select_events(rows, sv_types, interchromosomal_only, min_callers, max_events):
    """Choose which rows get breakpoint pages.

    The annotated table is annotation over the whole merged SV callset, so it
    holds deletions, insertions, inversions and duplications alongside the
    rearrangements. Rendering a page per breakpoint for all of them produces
    thousands of files and gigabytes of embedded alignment data, none of which
    is the evidence class this report is about.

    Selection is by SV type, optionally restricted to interchromosomal events,
    with a floor on caller count and a cap on the total. When the cap bites,
    events are kept in descending order of supporting reads so the
    best-evidenced survive.

    A row carrying a tier bypasses every filter here and is exempt from the
    cap. The actionable and defining calls in this assay are frequently
    intrachromosomal DEL or INV, so a type-and-chromosome filter drops exactly
    the events a reader needs the evidence for: PDGFRB::EBF1 (DEL plus its
    reciprocal INV pair) and P2RY8::CRLF2 (INV) had no pages at all. Graded
    rows are a handful per sample, so the page cost is negligible.

    Note on min_callers: the default is 1, deliberately. Single-caller
    rearrangements at low read support are not noise to be filtered away in
    this assay; they are the calls this panel exists to recover. Raising this
    is a decision to be made per run, not a default.

    Returns (selected_rows, summary_dict).
    """
    wanted = {t.strip().upper() for t in sv_types.split(",") if t.strip()}
    counts = {
        "total": len(rows),
        "by_type": {},
        "excluded_type": 0,
        "excluded_intrachromosomal": 0,
        "excluded_callers": 0,
        "excluded_cap": 0,
        "kept_graded": 0,
    }

    selected = []
    graded = []
    for row in rows:
        sv_type = (row.get("sv_type") or "").strip().upper()
        counts["by_type"][sv_type or "(blank)"] = (
            counts["by_type"].get(sv_type or "(blank)", 0) + 1
        )

        if (row.get("tier") or "").strip():
            counts["kept_graded"] += 1
            graded.append(row)
            continue

        if wanted and sv_type not in wanted:
            counts["excluded_type"] += 1
            continue

        chrom_a = (row.get("chrom_a") or "").strip()
        chrom_b = (row.get("chrom_b") or "").strip()
        if interchromosomal_only and chrom_a and chrom_b and chrom_a == chrom_b:
            counts["excluded_intrachromosomal"] += 1
            continue

        try:
            n_callers = int(float(row.get("n_callers") or 0))
        except ValueError:
            n_callers = 0
        if n_callers < min_callers:
            counts["excluded_callers"] += 1
            continue

        selected.append(row)

    def support_of(row):
        try:
            return float(row.get("support_reads") or 0)
        except ValueError:
            return 0.0

    if max_events and len(selected) > max_events:
        selected.sort(key=support_of, reverse=True)
        counts["excluded_cap"] = len(selected) - max_events
        selected = selected[:max_events]

    # Graded rows lead, and are never capped away.
    selected = graded + selected

    counts["selected"] = len(selected)
    return selected, counts


def report_selection(counts, sv_types, interchromosomal_only, min_callers):
    """Print why each row was kept or dropped, so the cost is never a surprise."""
    eprint("Event selection:")
    eprint("  rows in table            : %d" % counts["total"])
    by_type = ", ".join(
        "%s=%d" % (k, v)
        for k, v in sorted(counts["by_type"].items(), key=lambda kv: -kv[1])
    )
    eprint("  by sv_type               : %s" % by_type)
    eprint("  keeping sv_type          : %s" % sv_types)
    eprint("  interchromosomal only    : %s" % interchromosomal_only)
    eprint("  kept, graded (bypass)    : %d" % counts.get("kept_graded", 0))
    eprint("  minimum callers          : %d" % min_callers)
    eprint("  dropped, wrong sv_type   : %d" % counts["excluded_type"])
    eprint("  dropped, same chromosome : %d" % counts["excluded_intrachromosomal"])
    eprint("  dropped, below caller min: %d" % counts["excluded_callers"])
    if counts["excluded_cap"]:
        eprint("  dropped, over cap        : %d" % counts["excluded_cap"])
    eprint(
        "  SELECTED                 : %d events -> %d pages"
        % (counts["selected"], counts["selected"] * 2)
    )


def run_translocations(args):
    header, rows = read_tsv(args.sites_tsv)
    title = "%s translocations" % args.sample
    outdir = args.out_dir or os.path.dirname(os.path.abspath(args.out_html))
    os.makedirs(outdir, exist_ok=True)

    manifest_path = os.path.join(
        outdir, "%s.translocations.manifest.json" % args.sample
    )

    if not rows:
        write_placeholder(
            args.out_html,
            title,
            "This sample has no annotated structural variants. An empty "
            "result is an expected outcome, not a pipeline failure.",
        )
        with open(manifest_path, "w") as handle:
            json.dump(
                {
                    "sample": args.sample,
                    "flanking": args.flanking,
                    "events": [],
                },
                handle,
                indent=2,
            )
        return 0

    required = [column for _s, column, _p, _g in TRANSLOCATION_SIDES]
    required += [column for _s, _c, column, _g in TRANSLOCATION_SIDES]
    missing = [column for column in required if column not in header]
    if missing:
        fail(
            "annotated table lacks breakpoint columns %s: %s"
            % (", ".join(missing), args.sites_tsv)
        )

    ensure_bam_index(args.bam)

    selected, counts = select_events(
        rows,
        args.sv_types,
        args.interchromosomal_only,
        args.min_callers,
        args.max_events,
    )
    report_selection(
        counts, args.sv_types, args.interchromosomal_only, args.min_callers
    )

    if not selected:
        write_placeholder(
            args.out_html,
            title,
            "No rearrangements met the selection criteria. Other structural "
            "variant classes are present in the annotated table but are not "
            "rendered here.",
        )
        with open(manifest_path, "w") as handle:
            json.dump(
                {
                    "sample": args.sample,
                    "flanking": args.flanking,
                    "selection": counts,
                    "events": [],
                },
                handle,
                indent=2,
            )
        return 0

    info_columns = [c for c in TRANSLOCATION_INFO_COLUMNS if c in header]
    for extra in ("breakpoint", "partner_locus", "partner_gene"):
        if extra not in info_columns:
            info_columns.append(extra)

    workdir = tempfile.mkdtemp(prefix="igv_tx_")
    events = []
    failures = 0

    try:
        for order, row in enumerate(selected, start=1):
            event_id = row.get(TRANSLOCATION_ID) or "event_%d" % order
            stem = safe_name(event_id)
            record = {"event_id": event_id, "order": order}

            for side, chrom_key, pos_key, gene_key in TRANSLOCATION_SIDES:
                chrom = row.get(chrom_key, "")
                position = row.get(pos_key, "")
                if not chrom or not position:
                    continue
                try:
                    start, end = site_bounds(position)
                except ValueError:
                    continue

                other = "b" if side == "a" else "a"
                other_map = {
                    s: (c, p, g) for s, c, p, g in TRANSLOCATION_SIDES
                }
                other_chrom_key, other_pos_key, other_gene_key = other_map[other]

                site = dict(row)
                site["chrom"] = chrom
                site["start"] = start
                site["end"] = end
                site["breakpoint"] = "%s (%s:%s)" % (
                    side.upper(),
                    chrom,
                    position,
                )
                site["partner_locus"] = "%s:%s" % (
                    row.get(other_chrom_key, ""),
                    row.get(other_pos_key, ""),
                )
                site["partner_gene"] = row.get(other_gene_key, "")

                page_name = "%s.%s.html" % (stem, side.upper())
                page_path = os.path.join(outdir, page_name)
                sites_path = os.path.join(
                    workdir, "%s_%s.tsv" % (stem, side)
                )
                write_sites([site], info_columns, sites_path)

                page_title = "%s  %s  %s:%s" % (
                    args.sample,
                    row.get(gene_key, "") or chrom,
                    chrom,
                    position,
                )
                ok = run_create_report(
                    sites_path,
                    info_columns,
                    args.bam,
                    args.fasta,
                    page_path,
                    args.flanking,
                    page_title,
                )
                if not ok:
                    failures += 1
                    continue

                record[side] = {
                    "chrom": chrom,
                    "pos": position,
                    "label": row.get(gene_key, ""),
                    "locus": "%s:%d-%d"
                    % (
                        chrom,
                        max(1, int(float(position)) - args.flanking),
                        int(float(position)) + args.flanking,
                    ),
                    "html": page_name,
                }

            if "a" in record or "b" in record:
                events.append(record)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)

    with open(manifest_path, "w") as handle:
        json.dump(
            {
                "sample": args.sample,
                "flanking": args.flanking,
                "n_events": len(events),
                "selection": counts,
                "events": events,
            },
            handle,
            indent=2,
        )

    # An index page keeps the directory browsable on its own.
    write_index(args.out_html, args.sample, events)

    print(
        "translocations\t%s\t%d events\t%d pages\t%s"
        % (
            args.sample,
            len(events),
            sum(1 for e in events for s in ("a", "b") if s in e),
            outdir,
        )
    )
    if failures:
        eprint("%d breakpoint page(s) failed to render" % failures)
        return 2
    return 0


def write_index(out_html, sample, events):
    """Write a plain index of the generated breakpoint pages."""
    parts = [
        "<!DOCTYPE html><html lang='en'><head><meta charset='utf-8'>",
        "<title>%s translocation breakpoints</title>" % sample,
        "<style>body{font-family:system-ui,-apple-system,'Segoe UI',sans-serif;"
        "margin:2rem;color:#333}table{border-collapse:collapse}"
        "td,th{border:1px solid #ddd;padding:.35rem .6rem;font-size:.9rem}"
        "th{background:#f5f5f5;text-align:left}</style></head><body>",
        "<h2>%s &mdash; translocation breakpoints</h2>" % sample,
        "<table><tr><th>Event</th><th>Breakpoint A</th>"
        "<th>Breakpoint B</th></tr>",
    ]
    for event in events:
        cells = []
        for side in ("a", "b"):
            entry = event.get(side)
            if entry:
                cells.append(
                    "<a href='%s'>%s %s:%s</a>"
                    % (
                        entry["html"],
                        entry.get("label", ""),
                        entry["chrom"],
                        entry["pos"],
                    )
                )
            else:
                cells.append("&mdash;")
        parts.append(
            "<tr><td>%s</td><td>%s</td><td>%s</td></tr>"
            % (event["event_id"], cells[0], cells[1])
        )
    parts.append("</table></body></html>")

    os.makedirs(os.path.dirname(os.path.abspath(out_html)) or ".", exist_ok=True)
    with open(out_html, "w") as handle:
        handle.write("\n".join(parts))


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def build_parser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--mode",
        required=True,
        choices=("somatic", "translocations"),
        help="evidence class to render",
    )
    parser.add_argument("--sample", required=True, help="sequencing identifier")
    parser.add_argument(
        "--sites-tsv",
        required=True,
        help="clinical.tsv (somatic) or mm_annotated.tsv (translocations)",
    )
    parser.add_argument("--bam", required=True, help="alignment for this mode")
    parser.add_argument("--fasta", required=True, help="matching reference")
    parser.add_argument(
        "--out-html",
        required=True,
        help="report page (somatic) or index page (translocations)",
    )
    parser.add_argument(
        "--out-dir",
        default=None,
        help="directory for per-breakpoint pages; defaults to the directory "
        "of --out-html (translocations mode only)",
    )
    parser.add_argument(
        "--flanking",
        type=int,
        default=5000,
        help="bases of context either side of each site (default 5000)",
    )
    parser.add_argument(
        "--sv-types",
        default="TRA",
        help="comma-separated sv_type values to render, case-insensitive "
        "(default TRA). The annotated table also holds DEL, INS, INV and DUP, "
        "which are not rearrangements and are not rendered by default. Pass "
        "an empty string to render every type.",
    )
    parser.add_argument(
        "--all-sv-types",
        dest="sv_types",
        action="store_const",
        const="",
        help="render every sv_type; equivalent to --sv-types ''",
    )
    parser.add_argument(
        "--intrachromosomal",
        dest="interchromosomal_only",
        action="store_false",
        help="also render events whose breakpoints share a chromosome",
    )
    parser.set_defaults(interchromosomal_only=True)
    parser.add_argument(
        "--min-callers",
        type=int,
        default=1,
        help="minimum n_callers to render (default 1). Single-caller "
        "rearrangements at low read support are retained deliberately.",
    )
    parser.add_argument(
        "--max-events",
        type=int,
        default=200,
        help="cap on events rendered, highest supporting-read count first; "
        "0 means no cap (translocations only)",
    )
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)

    if not os.path.isfile(args.sites_tsv):
        fail("input table not found: %s" % args.sites_tsv)
    if not os.path.isfile(args.bam):
        fail("BAM not found: %s" % args.bam)
    if not os.path.isfile(args.fasta):
        fail("reference FASTA not found: %s" % args.fasta)
    if not shutil.which("create_report"):
        fail("create_report is not on PATH (activate the awgs_sv environment)")

    if args.mode == "somatic":
        return run_somatic(args)
    return run_translocations(args)


if __name__ == "__main__":
    sys.exit(main())
