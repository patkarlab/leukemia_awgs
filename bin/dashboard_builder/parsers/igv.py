"""Extract the IGV-reports embedded tableJson and build a chr:pos:ref:alt -> unique_id lookup.

igv-reports renders an HTML page that embeds, as a <script> blob, a JS object:

    const tableJson = {"headers": ["unique_id", "CHROM", "POSITION", "REF", "ALT", ...],
                       "rows": [[0, "chr1", 92478757, "C", "CAGAG", ...], ...]}

The rendered table assigns each row id="row_<unique_id>" and an onclick handler that
loads the corresponding IGV session. To cross-link from the parent dashboard's
clinical-variants table to a specific IGV row, we need that lookup at build time.

This module:
  1. Reads the igv_report.html
  2. Extracts the tableJson literal via BeautifulSoup + a small regex
  3. Parses it as JSON
  4. Builds {f"{chrom}:{pos}:{ref}:{alt}": unique_id}

It also injects a tiny, idempotent hash-router <script> at the bottom of the
report. The router listens for window.location.hash changes inside the iframe
and clicks the matching row -- which lets the parent dashboard select a variant
by setting iframe.src to "...#row_<uid>". Cross-origin navigation (setting src
on a sibling iframe) is allowed even under file://; reading the iframe's DOM is
not. This is the only way to wire variant-card -> IGV selection that works for
both file:// and http:// loads.

Returns None if the file is absent or tableJson cannot be parsed.
"""

import json
import logging
import re
from pathlib import Path

from bs4 import BeautifulSoup


_TABLEJSON_RE = re.compile(r"const\s+tableJson\s*=\s*(\{.*?\})\s*;?\s*$", re.DOTALL | re.MULTILINE)


# Sentinel comments wrap the injected script so the build can find and replace
# its own prior injection (idempotent re-runs) without corrupting the IGV
# report's actual content.
_HASH_ROUTER_BEGIN = "<!-- BEGIN tspipe-dashboard-builder hash-router (idempotent) -->"
_HASH_ROUTER_END   = "<!-- END tspipe-dashboard-builder hash-router -->"

_HASH_ROUTER_SCRIPT = """
<script>
// tspipe-igv-hash-router
// Select a variant in this IGV report by URL hash. The parent dashboard
// navigates iframe.src to "<report>#row_<uid>"; this listener (running INSIDE
// the iframe -- so always same-origin to the report it patches) clicks the
// matching row once IGV.js has finished initializing. Works under file://
// because the parent never has to read this iframe's DOM.
(function () {
  function rowIdFromHash() {
    var h = window.location.hash || "";
    if (h.indexOf("#row_") !== 0) return null;
    return h.substring(1);  // drop the '#', keep "row_<uid>"
  }

  function trySelect() {
    var rowId = rowIdFromHash();
    if (!rowId) return true;  // nothing to do, stop polling
    if (typeof igvBrowser === "undefined") return false;  // IGV not ready
    var row = document.getElementById(rowId);
    if (!row) {
      // Hash points at a row that doesn't exist in this report's table.
      // Give it a few seconds in case the table is still rendering, then
      // log and stop (handled by the polling timeout).
      return false;
    }
    // The original onclick reads event.target.parentElement.id expecting a
    // <td>. Click the first child cell so parentElement is the <tr>.
    var firstCell = row.firstElementChild;
    if (firstCell) {
      firstCell.click();
    } else {
      row.click();
    }
    row.scrollIntoView({behavior: "smooth", block: "center"});
    return true;
  }

  function startPolling() {
    if (!rowIdFromHash()) return;
    if (trySelect()) return;
    var start = Date.now();
    var timer = setInterval(function () {
      if (trySelect()) {
        clearInterval(timer);
      } else if (Date.now() - start > 30000) {
        clearInterval(timer);
        var rowId = rowIdFromHash();
        console.warn("tspipe-igv-hash-router: gave up waiting after 30s. " +
                     "Hash=" + window.location.hash + ", row " +
                     (rowId && document.getElementById(rowId) ? "present" : "missing") +
                     ", igvBrowser " +
                     (typeof igvBrowser !== "undefined" ? "ready" : "not ready") + ".");
      }
    }, 200);
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", startPolling);
  } else {
    startPolling();
  }
  window.addEventListener("hashchange", startPolling);
})();
</script>
"""


def inject_hash_router(igv_report_path):
    """Idempotently inject the hash-router script into an IGV report HTML.

    Strategy: locate (and remove) any prior injection block delimited by the
    sentinel comments, then insert the current block immediately before the
    closing </body> tag. If there is no </body>, append at the end of the
    file (graceful degradation -- browsers still execute the script).

    Returns True if the file was modified, False if no injection was needed
    (e.g. the report was missing) or the file is already current.
    """
    path = Path(igv_report_path)
    if not path.exists():
        return False

    try:
        with open(path, "r", encoding="utf-8") as fh:
            html = fh.read()
    except OSError as exc:
        logging.warning("Could not read IGV report for patching: %s", exc)
        return False

    block = _HASH_ROUTER_BEGIN + _HASH_ROUTER_SCRIPT + _HASH_ROUTER_END

    # Strip any prior injection (re-running the build should refresh it).
    if _HASH_ROUTER_BEGIN in html:
        pattern = re.escape(_HASH_ROUTER_BEGIN) + r".*?" + re.escape(_HASH_ROUTER_END)
        stripped = re.sub(pattern, "", html, flags=re.DOTALL)
        # If stripping yields the same html (post-trim), nothing to do --
        # but we still want to re-inject the current version of the script
        # in case it changed between builder versions.
        html = stripped

    # Insert immediately before </body> if present; else append.
    body_close = re.search(r"</body\s*>", html, re.IGNORECASE)
    if body_close:
        idx = body_close.start()
        new_html = html[:idx] + "\n" + block + "\n" + html[idx:]
    else:
        new_html = html.rstrip() + "\n" + block + "\n"

    try:
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(new_html)
    except OSError as exc:
        logging.warning("Could not write patched IGV report: %s", exc)
        return False

    return True


def extract_lookup(igv_report_path):
    path = Path(igv_report_path)
    if not path.exists():
        return None

    try:
        with open(path, "r", encoding="utf-8") as fh:
            html = fh.read()
    except OSError:
        return None

    soup = BeautifulSoup(html, "html.parser")
    # tableJson lives inside a <script> tag — search all of them.
    table_json = None
    for tag in soup.find_all("script"):
        if not tag.string:
            continue
        text = tag.string
        if "tableJson" not in text:
            continue
        # The regex anchors on 'const tableJson = {...};' on a single line in the source.
        # In practice the entire object is on one line; match a balanced top-level JSON.
        # Find first '{' after 'tableJson' and balance braces ourselves to be safe.
        start = text.find("tableJson")
        brace_start = text.find("{", start)
        if brace_start == -1:
            continue
        depth = 0
        end = -1
        in_str = False
        esc = False
        for i in range(brace_start, len(text)):
            c = text[i]
            if in_str:
                if esc:
                    esc = False
                elif c == "\\":
                    esc = True
                elif c == '"':
                    in_str = False
                continue
            if c == '"':
                in_str = True
            elif c == "{":
                depth += 1
            elif c == "}":
                depth -= 1
                if depth == 0:
                    end = i + 1
                    break
        if end == -1:
            continue
        try:
            table_json = json.loads(text[brace_start:end])
            break
        except json.JSONDecodeError:
            continue

    if not table_json:
        return None

    headers = table_json.get("headers", [])
    rows = table_json.get("rows", [])
    try:
        idx_id    = headers.index("unique_id")
        idx_chrom = headers.index("CHROM")
        idx_pos   = headers.index("POSITION")
        idx_ref   = headers.index("REF")
        idx_alt   = headers.index("ALT")
    except ValueError:
        return None

    lookup = {}
    for row in rows:
        try:
            key = f"{row[idx_chrom]}:{row[idx_pos]}:{row[idx_ref]}:{row[idx_alt]}"
            lookup[key] = row[idx_id]
        except (IndexError, TypeError):
            continue

    return {"lookup": lookup, "n": len(lookup)}
