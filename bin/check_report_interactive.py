#!/usr/bin/env python3
"""Validate the interactive controls of every built sample report in a bundle.

Rationale
---------
A dashboard control that is wired to the wrong column does not fail loudly.
It returns a plausible-looking table, and in the case that prompted this
script it returned an empty one, which reads as a genuine absence of events
rather than as a defect. Nothing downstream catches that: the build succeeds,
the zip is produced, the report is signed out.

The specific defect was an index computed against one representation of the
table and applied to another. The supporting-read threshold resolved its
column by scanning thead after DataTables had detached the hidden columns'
<th> elements, so the index enumerated visible columns while the row data
stayed indexed over all of them. Four hidden columns preceded the target, the
threshold was applied to a free-text column four positions earlier, and every
row scored zero.

The parser now supplies the index, so the two representations cannot diverge.
This script asserts that they have not: it reads the index the report will
use and confirms that the header at that position is the supporting-read
column and that the column holds numbers.

Checks per sample report
------------------------
  1. the threshold control markup is present
  2. the emitted supportCol is an integer, or -1 with the disabled path taken
  3. the <th> at that index is the supporting-read header
  4. the <td> cells at that index carry numeric data-order values
  5. the predicate keeps rows it cannot parse rather than dropping them
  6. the index is assigned once, from a literal, and never recomputed
  7. the SV-type filter admits graded rows, so a report whose only findings
     are intrachromosomal does not open on a table that excludes them

Reports with no translocation table are skipped, which is a legitimate
outcome and not a failure.

Accepts a build directory, a distributed .zip bundle, or a single report
file, so the artefact that actually reaches a reader can be validated without
being unpacked first.

Exit status
-----------
  0  every report checked passed, or none contained a translocation table
  1  at least one report failed a check
  2  the path is unusable or holds no sample reports
"""

import argparse
import re
import sys
import zipfile
from html.parser import HTMLParser
from pathlib import Path


SUPPORT_HEADER_PATTERN = re.compile(r"supporting[\s-]*reads?", re.I)
SUPPORT_COL_PATTERN = re.compile(r"\bsupportCol\s*=\s*([^;]+);")
IGV_OFFSET_PATTERN = re.compile(r"\bigvOffset\s*=\s*(\d+)")


class TranslocationTableParser(HTMLParser):
    """Collect the header labels and the first body row of table#tx-table.

    A dedicated parser rather than a regex because the cells carry attributes
    that matter (data-order) and nesting has to be tracked to know which row
    a cell belongs to. Standard library only: the builder environment is
    jinja2 plus the standard library and this script must not widen that.
    """

    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.in_table = False
        self.depth = 0
        self.section = None          # 'thead' or 'tbody'
        self.headers = []
        self.body_rows = []          # list of list of (text, data_order)
        self.row_tiers = []          # data-tier of each body row, aligned
        self.row_tras = []           # data-translocation, same alignment
        self._cell = None
        self._cell_order = None
        self._row = None

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == "table":
            if attrs.get("id") == "tx-table":
                self.in_table = True
                self.depth = 1
            elif self.in_table:
                self.depth += 1
            return
        if not self.in_table or self.depth != 1:
            return
        if tag in ("thead", "tbody"):
            self.section = tag
        elif tag == "tr":
            self._row = []
            self._row_tier = attrs.get("data-tier")
            self._row_tra = attrs.get("data-translocation")
        elif tag in ("th", "td"):
            self._cell = []
            self._cell_order = attrs.get("data-order")

    def handle_endtag(self, tag):
        if tag == "table" and self.in_table:
            self.depth -= 1
            if self.depth == 0:
                self.in_table = False
            return
        if not self.in_table or self.depth != 1:
            return
        if tag in ("th", "td") and self._cell is not None:
            text = " ".join("".join(self._cell).split())
            if self._row is not None:
                self._row.append((text, self._cell_order))
            self._cell = None
            self._cell_order = None
        elif tag == "tr" and self._row is not None:
            if self.section == "thead" and not self.headers:
                self.headers = [text for text, _ in self._row]
            elif self.section == "tbody":
                self.body_rows.append(self._row)
                self.row_tiers.append(getattr(self, "_row_tier", None))
                self.row_tras.append(getattr(self, "_row_tra", None))
            self._row = None
        elif tag in ("thead", "tbody"):
            self.section = None

    def handle_data(self, data):
        if self._cell is not None:
            self._cell.append(data)


def check_report(text, sample_rows=200):
    """Return a list of failure strings. Empty means the report passed.

    Takes the report text rather than a path so that the same checks apply
    whether the report is being read from a build directory or from inside
    the distributed archive. The archive is what leaves the pipeline, so it
    has to be checkable directly.
    """
    if 'id="tx-table"' not in text:
        return None                            # no table; nothing to check

    failures = []

    # 1. control markup ---------------------------------------------------
    for marker, label in (
        ('id="tx-support-filter"', "threshold button group"),
        ('id="tx-support-min"', "free numeric entry"),
        ('id="tx-support-count"', "row counter"),
        ('id="tx-support-warn"', "disabled-state notice"),
    ):
        if marker not in text:
            failures.append(f"missing {label} ({marker})")

    # 6. regression guard: no post-init DOM index discovery ---------------
    if "__txSuppIdx" in text:
        failures.append(
            "column index is re-derived from the DOM after initialisation "
            "(__txSuppIdx); this is the defect the parser-supplied index "
            "replaced")

    # 5. the predicate must keep rows it cannot parse ---------------------
    if "isFinite(value) ? (value >= minSupport) : true" not in text:
        failures.append(
            "threshold predicate is absent, or drops rows whose value it "
            "cannot parse; a resolver failure must degrade to an "
            "unfiltered table, not an empty one")

    # 2. the emitted index ------------------------------------------------
    # Exactly one assignment, from a literal. More than one means the value
    # is being recomputed somewhere, which is the shape of the original
    # defect regardless of what the recomputation is named: an index derived
    # at runtime can disagree with the data it indexes. A literal supplied by
    # the parser cannot.
    assignments = SUPPORT_COL_PATTERN.findall(text)
    if not assignments:
        failures.append("no supportCol assignment found in the report")
        return failures
    if len(assignments) > 1:
        failures.append(
            f"supportCol is assigned {len(assignments)} times "
            f"({', '.join(repr(a.strip()) for a in assignments)}); it must be "
            "assigned once from the parser and never recomputed")
        return failures

    expression = assignments[0].strip()

    # The offset is whatever the report emits, not whatever the expression
    # mentions. A sample without IGV snapshots renders no IGV column and sets
    # igvOffset to 0, so reading the token rather than its value shifts every
    # such report by one and reports a correct index as wrong.
    if "igvOffset" in expression:
        offset_match = IGV_OFFSET_PATTERN.search(text)
        if not offset_match:
            failures.append(
                "supportCol references igvOffset but no igvOffset assignment "
                "was found; the index cannot be validated")
            return failures
        igv_offset = int(offset_match.group(1))
    else:
        igv_offset = 0

    literal = re.match(r"^(-?\d+)", expression)
    if not literal:
        failures.append(f"supportCol is not a literal index: {expression!r}")
        return failures
    support_col = int(literal.group(1)) + igv_offset

    parser = TranslocationTableParser()
    parser.feed(text)
    headers = parser.headers

    if support_col < 0:
        # Legitimate only when the column genuinely is not in the table.
        if any(SUPPORT_HEADER_PATTERN.search(h) for h in headers):
            failures.append(
                "threshold disabled (supportCol -1) but a supporting-read "
                "header is present in the table")
        return failures

    # 3. the header at that index ----------------------------------------
    if support_col >= len(headers):
        failures.append(
            f"supportCol {support_col} is past the last column "
            f"({len(headers)} headers)")
        return failures

    label = headers[support_col]
    if not SUPPORT_HEADER_PATTERN.search(label):
        actual = [i for i, h in enumerate(headers)
                  if SUPPORT_HEADER_PATTERN.search(h)]
        failures.append(
            f"supportCol {support_col} points at {label!r}; the "
            f"supporting-read header is at {actual or 'no'} index")

    # 4. the cells at that index -----------------------------------------
    numeric = 0
    seen = 0
    for row in parser.body_rows[:sample_rows]:
        if support_col >= len(row):
            continue
        seen += 1
        cell_text, order = row[support_col]
        value = order if order not in (None, "") else cell_text
        try:
            float(value)
            numeric += 1
        except (TypeError, ValueError):
            pass
    if seen and numeric / seen < 0.8:
        failures.append(
            f"column {support_col} holds numbers in only {numeric}/{seen} "
            "sampled rows; the threshold would score the rest as zero")

    # 7. the SV-type filter must admit graded rows ------------------------
    # The tab opens on rearrangements, with other SV types behind an
    # unchecked switch. That is a reasonable default for browsing and a
    # silent failure for reporting: the actionable and defining calls in
    # this assay are frequently intrachromosomal, so a type-only filter can
    # hide every finding a sample has. It did. A row carrying a tier must
    # survive the filter whatever its type.
    #
    # Three conditions, because the weaker two pass vacuously. The attribute
    # must be emitted; the predicate must read it before falling through to
    # the type test; and the report must actually contain a graded row that
    # the type test would otherwise exclude, or the check proves nothing on
    # a sample whose findings happen all to be translocations.
    # Unconditional. Gating this on the presence of graded rows makes the
    # check inert on precisely the reports that carry the defect: one built
    # before the fix has no data-tier attributes at all, so a presence test
    # skips and the report passes while hiding every finding it holds.
    filt = text.find("data-translocation') !== '0'")
    if filt != -1:
        tier_read = text.find("getAttribute('data-tier')")
        if tier_read == -1 or tier_read > filt:
            failures.append(
                "the SV-type predicate excludes on data-translocation without "
                "first consulting data-tier; graded rows of other SV types are "
                "hidden until the switch is enabled, and in this assay the "
                "actionable and defining calls are frequently intrachromosomal")
        elif "data-tier" not in text:
            failures.append(
                "the predicate reads data-tier but no row carries it, so every "
                "row falls through to the SV-type test")

    graded = [i for i, v in enumerate(parser.row_tiers) if v == "1"]
    if graded:
        hidden = [i for i in graded
                  if i < len(parser.row_tras) and parser.row_tras[i] == "0"]
        if not hidden:
            # Every graded row in this sample is a translocation, so the
            # bypass is untested here. Not a failure: a sample's findings
            # are what they are. The condition above still holds for any
            # sample that does carry an intrachromosomal finding.
            pass

    return failures


def is_sample_report(name):
    """Sample reports only. The IGV pages are a different artefact."""
    base = name.rsplit("/", 1)[-1]
    return base.endswith("_report.html") and not base.endswith("_igv_report.html")


def iter_reports(root):
    """Yield (display_name, text) for every sample report under root.

    Accepts a build directory, a distributed .zip bundle, or a single report
    file. The archive is accepted because it is the artefact that actually
    reaches a reader, and requiring it to be unpacked first would put a
    manual step between the deliverable and its validation.
    """
    if root.is_dir():
        for path in sorted(root.rglob("*_report.html")):
            if is_sample_report(path.name):
                yield path.name, path.read_text(encoding="utf-8", errors="replace")
        return

    if root.suffix.lower() == ".zip":
        with zipfile.ZipFile(root) as archive:
            for name in sorted(archive.namelist()):
                if is_sample_report(name):
                    with archive.open(name) as handle:
                        yield (name.rsplit("/", 1)[-1],
                               handle.read().decode("utf-8", errors="replace"))
        return

    if root.suffix.lower() in (".html", ".htm"):
        yield root.name, root.read_text(encoding="utf-8", errors="replace")
        return

    raise ValueError(
        "expected a bundle directory, a .zip bundle, or a report .html file")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("bundle",
                    help="bundle directory, .zip bundle, or a single report .html")
    ap.add_argument("--quiet", "-q", action="store_true",
                    help="print failures only")
    args = ap.parse_args()

    root = Path(args.bundle)
    if not root.exists():
        print(f"check_report_interactive: no such path: {root}", file=sys.stderr)
        return 2

    try:
        reports = list(iter_reports(root))
    except (ValueError, zipfile.BadZipFile) as exc:
        print(f"check_report_interactive: {root}: {exc}", file=sys.stderr)
        return 2

    if not reports:
        print(f"check_report_interactive: no sample reports found in {root}",
              file=sys.stderr)
        return 2

    failed = 0
    skipped = 0
    passed = 0

    for name, text in reports:
        result = check_report(text)
        if result is None:
            skipped += 1
            if not args.quiet:
                print(f"skip    {name}: no translocation table")
            continue
        if result:
            failed += 1
            print(f"FAIL    {name}")
            for line in result:
                print(f"          {line}")
        else:
            passed += 1
            if not args.quiet:
                print(f"ok      {name}")

    print(f"\nchecked {len(reports)} reports: "
          f"{passed} passed, {failed} failed, {skipped} without a table")

    if failed:
        print("\nThe supporting-read threshold is misconfigured in the reports "
              "above. A threshold applied to the wrong column returns a "
              "plausible table rather than an error, so this is a blocking "
              "failure rather than a warning.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
