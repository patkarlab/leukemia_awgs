#!/usr/bin/env python3
"""
embed_report_assets.py

Rewrite generated report HTML so every local dependency is carried inside the
file, leaving a single document that renders identically wherever it is opened.

Why this exists
---------------
The builder emits HTML that references its stylesheets, scripts and figures by
relative path (assets/css/..., baf_loh/figures/....png). That works on the
server where the bundle was built. It stops working the moment a report is
copied out of the bundle, mailed, or extracted to a different depth on a
Windows machine, and the failure is silent: a broken image icon or an unstyled
page, with no indication of what went missing.

Embedding removes the class of failure rather than any single instance of it.

What gets inlined
-----------------
- <link rel="stylesheet" href="..."> becomes a <style> block, with any url()
  references inside the CSS resolved and inlined in turn
- <script src="..."> becomes an inline <script>
- <img src="...">, and CSS url(...), become base64 data URIs
- Absolute URLs (http, https, protocol-relative), existing data: URIs and
  in-page anchors are left untouched

What deliberately does not get inlined
--------------------------------------
<iframe src="..."> is left alone. The IGV breakpoint pages are whole documents
of their own, each already carrying its own embedded alignment slice; folding
them into the parent would multiply the page size by the number of events. They
travel alongside the report instead, and the script reports whether they are
present so a missing IGV tree is noticed here rather than by a reader.

Standard library only.

Exit codes
----------
0  success
2  usage or input error
"""

import argparse
import base64
import mimetypes
import os
import re
import sys


# Attributes carrying a local dependency, by tag.
LINK_PATTERN = re.compile(
    r"""<link\b[^>]*?\brel\s*=\s*["']?stylesheet["']?[^>]*?>""",
    re.IGNORECASE,
)
SCRIPT_PATTERN = re.compile(
    r"""<script\b([^>]*?)\bsrc\s*=\s*["']([^"']+)["']([^>]*?)>\s*</script\s*>""",
    re.IGNORECASE,
)
IMG_PATTERN = re.compile(
    r"""(<img\b[^>]*?\bsrc\s*=\s*["'])([^"']+)(["'])""",
    re.IGNORECASE,
)
IFRAME_PATTERN = re.compile(
    r"""<iframe\b[^>]*?\bsrc\s*=\s*["']([^"']+)["']""",
    re.IGNORECASE,
)
HREF_IN_TAG = re.compile(r"""\bhref\s*=\s*["']([^"']+)["']""", re.IGNORECASE)
CSS_URL = re.compile(r"""url\(\s*(['"]?)([^'")]+)\1\s*\)""", re.IGNORECASE)

SKIP_PREFIXES = ("data:", "http://", "https://", "//", "#", "mailto:", "about:")

# Text assets are inlined as text; everything else becomes a data URI.
TEXT_SUFFIXES = (".css", ".js")


class Stats(object):
    def __init__(self):
        self.css = 0
        self.js = 0
        self.images = 0
        self.css_urls = 0
        self.missing = []
        self.iframes_present = 0
        self.iframes_missing = []
        self.igv_embedded = 0
        self.directories = []


def is_external(url):
    lowered = url.strip().lower()
    return any(lowered.startswith(p) for p in SKIP_PREFIXES)


def note_unresolved(target, stats):
    """Record a reference that did not land on a file.

    A reference pointing at a directory is a malformed src or href in the
    template, not a file the bundle step forgot; separating the two stops the
    missing-asset list sending anyone after a file that never existed.
    """
    if os.path.isdir(target):
        stats.directories.append(target)
    else:
        stats.missing.append(target)


def resolve(base_dir, url):
    """Resolve a relative reference against the document directory."""
    cleaned = url.split("?", 1)[0].split("#", 1)[0]
    cleaned = cleaned.replace("\\", "/")
    return os.path.normpath(os.path.join(base_dir, cleaned))


def read_bytes(path):
    with open(path, "rb") as handle:
        return handle.read()


def data_uri(path):
    """Base64 data URI for a binary asset."""
    mime, _ = mimetypes.guess_type(path)
    if not mime:
        mime = "application/octet-stream"
    encoded = base64.b64encode(read_bytes(path)).decode("ascii")
    return "data:%s;base64,%s" % (mime, encoded)


def inline_css_urls(css_text, css_dir, stats, depth=0):
    """Inline url() references inside a stylesheet.

    Fonts and background images referenced from CSS are the assets most often
    forgotten, because nothing in the HTML mentions them.
    """
    if depth > 3:
        return css_text

    def replace(match):
        quote, url = match.group(1), match.group(2)
        if is_external(url):
            return match.group(0)
        target = resolve(css_dir, url)
        if not os.path.isfile(target):
            stats.missing.append(target)
            return match.group(0)
        stats.css_urls += 1
        return "url(%s%s%s)" % (quote, data_uri(target), quote)

    return CSS_URL.sub(replace, css_text)


def inline_stylesheets(html, base_dir, stats):
    def replace(match):
        tag = match.group(0)
        href_match = HREF_IN_TAG.search(tag)
        if not href_match:
            return tag
        href = href_match.group(1)
        if is_external(href):
            return tag
        target = resolve(base_dir, href)
        if not os.path.isfile(target):
            note_unresolved(target, stats)
            return tag
        try:
            css = read_bytes(target).decode("utf-8", "replace")
        except OSError:
            stats.missing.append(target)
            return tag
        css = inline_css_urls(css, os.path.dirname(target), stats)
        stats.css += 1
        return "<style>\n/* inlined from %s */\n%s\n</style>" % (
            os.path.basename(target),
            css,
        )

    return LINK_PATTERN.sub(replace, html)


def inline_scripts(html, base_dir, stats):
    def replace(match):
        before, src, after = match.group(1), match.group(2), match.group(3)
        if is_external(src):
            return match.group(0)
        target = resolve(base_dir, src)
        if not os.path.isfile(target):
            note_unresolved(target, stats)
            return match.group(0)
        try:
            code = read_bytes(target).decode("utf-8", "replace")
        except OSError:
            stats.missing.append(target)
            return match.group(0)
        # A literal </script> inside the source would terminate the block
        # early; the escaped form is equivalent to the parser.
        code = code.replace("</script", "<\\/script")
        stats.js += 1
        attrs = (before + after).replace("defer", "").replace("async", "").strip()
        opening = "<script %s>" % attrs if attrs else "<script>"
        return "%s\n/* inlined from %s */\n%s\n</script>" % (
            opening,
            os.path.basename(target),
            code,
        )

    return SCRIPT_PATTERN.sub(replace, html)


def inline_images(html, base_dir, stats):
    def replace(match):
        prefix, src, suffix = match.group(1), match.group(2), match.group(3)
        if is_external(src):
            return match.group(0)
        target = resolve(base_dir, src)
        if not os.path.isfile(target):
            note_unresolved(target, stats)
            return match.group(0)
        stats.images += 1
        return prefix + data_uri(target) + suffix

    return IMG_PATTERN.sub(replace, html)


IGV_BUTTON = re.compile(
    r"""(<button\b[^>]*?\bclass\s*=\s*["'][^"']*\btx-igv-btn\b[^"']*["'][^>]*?>)""",
    re.IGNORECASE,
)
IGV_ATTR = re.compile(
    r"""\bdata-html-(a|b)\s*=\s*["']([^"']*)["']""", re.IGNORECASE
)


def embed_igv_pages(html, base_dir, limit, stats):
    """Replace breakpoint page paths on IGV buttons with data: URIs.

    An igv-reports page is a whole self-contained document, so it can be
    carried inside the report as a base64 data URI and handed straight to the
    iframe. The viewer then works from the report file alone, with no sibling
    directory to lose in transit.

    The pages are large, so only the first ``limit`` buttons are inlined. The
    table is already sorted by supporting reads, so those are the
    best-evidenced events. Buttons beyond the limit keep their relative paths
    and still work while the igv/ tree is alongside.
    """
    if limit <= 0:
        return html

    state = {"done": 0}

    def replace_button(match):
        tag = match.group(1)
        if state["done"] >= limit:
            return tag

        replacements = {}
        for attr_match in IGV_ATTR.finditer(tag):
            side, rel = attr_match.group(1), attr_match.group(2)
            if not rel or is_external(rel):
                continue
            target = resolve(base_dir, rel)
            if not os.path.isfile(target):
                stats.missing.append(target)
                continue
            try:
                encoded = base64.b64encode(read_bytes(target)).decode("ascii")
            except OSError:
                stats.missing.append(target)
                continue
            replacements[side] = "data:text/html;base64," + encoded

        if not replacements:
            return tag

        def swap(attr_match):
            side = attr_match.group(1)
            if side in replacements:
                return 'data-html-%s="%s"' % (side, replacements[side])
            return attr_match.group(0)

        state["done"] += 1
        stats.igv_embedded += len(replacements)
        return IGV_ATTR.sub(swap, tag)

    return IGV_BUTTON.sub(replace_button, html)


def check_iframes(html, base_dir, stats):
    """Report iframe targets without inlining them."""
    for match in IFRAME_PATTERN.finditer(html):
        src = match.group(1)
        if is_external(src) or src.strip().lower().startswith("about:"):
            continue
        target = resolve(base_dir, src)
        if os.path.isfile(target):
            stats.iframes_present += 1
        else:
            stats.iframes_missing.append(target)


def process(path, backup, stats, embed_igv=0):
    base_dir = os.path.dirname(os.path.abspath(path))
    with open(path, encoding="utf-8", errors="replace") as handle:
        html = handle.read()
    original_size = len(html)

    html = inline_stylesheets(html, base_dir, stats)
    html = inline_scripts(html, base_dir, stats)
    html = inline_images(html, base_dir, stats)
    if embed_igv:
        html = embed_igv_pages(html, base_dir, embed_igv, stats)
    check_iframes(html, base_dir, stats)

    if backup and not os.path.exists(path + ".preembed"):
        with open(path + ".preembed", "w", encoding="utf-8") as handle:
            handle.write(open(path, encoding="utf-8", errors="replace").read())

    with open(path, "w", encoding="utf-8") as handle:
        handle.write(html)

    return original_size, len(html)


def find_reports(root):
    """Locate the report HTML files the builder produced."""
    found = []
    for dirpath, dirnames, filenames in os.walk(root):
        # The IGV pages are self-contained already and must not be rewritten.
        dirnames[:] = [d for d in dirnames if d != "igv"]
        for name in filenames:
            if not name.endswith(".html"):
                continue
            # <sample>_igv_report.html also ends in _report.html but is an
            # igv-reports page, already self-contained, and patched in place
            # by the builder. Rewriting it would be pointless at best.
            if name.endswith("_igv_report.html"):
                continue
            if name.endswith("_report.html") or name == "cohort_index.html":
                found.append(os.path.join(dirpath, name))
    return sorted(found)


def human(n):
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return "%.1f %s" % (n, unit)
        n /= 1024.0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("bundle", help="built report bundle directory")
    parser.add_argument(
        "--embed-igv",
        type=int,
        default=0,
        metavar="N",
        help="inline the breakpoint pages for the first N IGV buttons as "
        "data: URIs, so the paired viewer works from the report file alone. "
        "Costs roughly the size of 2N igv-reports pages. 0 disables.",
    )
    parser.add_argument(
        "--no-backup",
        action="store_true",
        help="do not write a .preembed copy of each rewritten file",
    )
    args = parser.parse_args(argv)

    if not os.path.isdir(args.bundle):
        print("ERROR: not a directory: %s" % args.bundle, file=sys.stderr)
        return 2

    reports = find_reports(args.bundle)
    if not reports:
        print("ERROR: no *_report.html or cohort_index.html under %s"
              % args.bundle, file=sys.stderr)
        return 2

    total_before = 0
    total_after = 0
    all_missing = []
    all_directories = []
    iframe_missing = []

    for path in reports:
        stats = Stats()
        before, after = process(path, not args.no_backup, stats, args.embed_igv)
        total_before += before
        total_after += after
        rel = os.path.relpath(path, args.bundle)
        print(
            "%-52s %8s -> %8s  css=%d js=%d img=%d cssurl=%d"
            % (rel, human(before), human(after), stats.css, stats.js,
               stats.images, stats.css_urls)
        )
        if stats.igv_embedded:
            print("%-52s igv pages embedded: %d" % ("", stats.igv_embedded))
        if stats.iframes_present or stats.iframes_missing:
            print(
                "%-52s iframes: %d present, %d missing"
                % ("", stats.iframes_present, len(stats.iframes_missing))
            )
        all_missing.extend(stats.missing)
        all_directories.extend(stats.directories)
        iframe_missing.extend(stats.iframes_missing)

    print("")
    print("Reports rewritten : %d" % len(reports))
    print("Total size        : %s -> %s" % (human(total_before), human(total_after)))

    if all_missing:
        unique = sorted(set(all_missing))
        print("")
        print("Assets referenced but not found (%d):" % len(unique))
        for path in unique[:20]:
            print("  " + os.path.relpath(path, args.bundle))
        if len(unique) > 20:
            print("  ... and %d more" % (len(unique) - 20))
        print("")
        print("These stay as relative references and will not render. They are")
        print("usually files the bundle step did not collect.")

    if all_directories:
        unique = sorted(set(all_directories))
        print("")
        print("References that resolve to a directory, not a file (%d):" % len(unique))
        for path in unique[:10]:
            print("  " + os.path.relpath(path, args.bundle))
        print("")
        print("These are malformed src or href values in a template, not files")
        print("the bundle is missing. Harmless, but they belong in the template.")

    if iframe_missing:
        unique = sorted(set(iframe_missing))
        print("")
        print("IGV pages referenced but not present (%d):" % len(unique))
        for path in unique[:10]:
            print("  " + os.path.relpath(path, args.bundle))
        if len(unique) > 10:
            print("  ... and %d more" % (len(unique) - 10))
        print("")
        print("The paired breakpoint viewer will show empty panes for these.")
        print("Rebuild the bundle with the igv/ tree present, or accept that")
        print("this archive is being distributed without IGV.")

    return 0


if __name__ == "__main__":
    sys.exit(main())
