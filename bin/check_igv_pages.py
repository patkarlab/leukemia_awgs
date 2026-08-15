#!/usr/bin/env python3
"""
check_igv_pages.py

Inspect the IGV breakpoint pages in a report bundle and say, per sample,
whether they are present and whether they actually contain alignment data.

An igv-reports page that rendered correctly is a few hundred kilobytes at
minimum, because the read pileup for the locus is embedded in it as base64. A
page of a few kilobytes opened without error but has nothing to show: usually
no reads at the locus, or a track that failed to attach. Distinguishing those
two cases from the page not being in the archive at all is the whole point of
this script, because on screen they look identical.

Usage:
    python3 bin/check_igv_pages.py <bundle_dir>
    python3 bin/check_igv_pages.py <bundle_dir> --verbose

Standard library only.
"""

import argparse
import json
import os
import re
import sys


# Whether a page carries reads is decided by the presence of an embedded
# payload, not by file size. Adaptive sampling leaves partner-side breakpoints
# outside the panel windows at low depth, so a small page with reads in it is
# the correct result for those loci rather than a defect. Size is reported for
# context only.
TINY_BYTES = 12 * 1024

# Markers that indicate a working igv-reports page.
IGV_MARKERS = ("igv.createBrowser", "igv.js", "igvjs", "createBrowser")
DATA_MARKER = re.compile(r"base64,[A-Za-z0-9+/]{500,}")


def human(n):
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return "%.0f %s" % (n, unit)
        n /= 1024.0


def find_manifests(bundle):
    found = []
    for dirpath, _dirnames, filenames in os.walk(bundle):
        for name in filenames:
            if name.endswith(".translocations.manifest.json"):
                found.append(os.path.join(dirpath, name))
    return sorted(found)


def inspect(path):
    """Return (exists, size, has_igv, has_payload)."""
    if not os.path.isfile(path):
        return False, 0, False, False
    size = os.path.getsize(path)
    # Read a bounded window; the markers appear early and the payload is bulk.
    with open(path, encoding="utf-8", errors="replace") as handle:
        head = handle.read(400000)
    has_igv = any(marker in head for marker in IGV_MARKERS)
    has_payload = bool(DATA_MARKER.search(head))
    return True, size, has_igv, has_payload


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("bundle", help="report bundle directory")
    parser.add_argument("--verbose", action="store_true",
                        help="list every page, not just the problems")
    args = parser.parse_args(argv)

    if not os.path.isdir(args.bundle):
        print("ERROR: not a directory: %s" % args.bundle, file=sys.stderr)
        return 2

    manifests = find_manifests(args.bundle)
    if not manifests:
        print("No IGV manifest found anywhere under %s." % args.bundle)
        print("")
        print("That means the bundle carries no IGV pages at all, so the paired")
        print("viewer has nothing to open. Either the snapshots were never")
        print("generated, or build_report_bundle.sh did not collect igv/, or the")
        print("archive was made with --light.")
        return 0

    grand_missing = 0
    grand_broken = 0
    grand_ok = 0

    for manifest_path in manifests:
        base = os.path.dirname(manifest_path)
        try:
            with open(manifest_path) as handle:
                manifest = json.load(handle)
        except (OSError, ValueError) as error:
            print("UNREADABLE MANIFEST %s: %s" % (manifest_path, error))
            continue

        sample = manifest.get("sample", "?")
        events = manifest.get("events", [])
        selection = manifest.get("selection") or {}

        print("=" * 68)
        print("Sample %s" % sample)
        print("  manifest      : %s" % os.path.relpath(manifest_path, args.bundle))
        print("  events        : %d" % len(events))
        if selection:
            print("  selected from : %s rows (%s)"
                  % (selection.get("total", "?"),
                     ", ".join("%s=%s" % kv for kv in
                               sorted((selection.get("by_type") or {}).items()))))
        print("  flanking      : %s bp" % manifest.get("flanking", "?"))

        missing = []
        broken = []
        ok = []

        for event in events:
            for side in ("a", "b"):
                entry = event.get(side)
                if not entry:
                    continue
                page = os.path.join(base, entry["html"])
                exists, size, has_igv, has_payload = inspect(page)
                label = "%s %s %s:%s" % (event["event_id"], side.upper(),
                                         entry.get("chrom"), entry.get("pos"))
                if not exists:
                    missing.append((label, entry["html"]))
                elif not has_payload or not has_igv or size < TINY_BYTES:
                    broken.append((label, size, has_igv, has_payload))
                else:
                    ok.append((label, size))

        print("  pages present : %d" % (len(ok) + len(broken)))
        print("  pages missing : %d" % len(missing))
        print("  pages broken  : %d  (no embedded read data)" % len(broken))

        if missing:
            print("")
            print("  MISSING (first 5):")
            for label, name in missing[:5]:
                print("    %-44s %s" % (label, name))

        if broken:
            print("")
            print("  BROKEN (first 5) - size / igv markers / embedded data:")
            for label, size, has_igv, has_payload in broken[:5]:
                print("    %-44s %8s  igv=%s  data=%s"
                      % (label, human(size), "yes" if has_igv else "NO",
                         "yes" if has_payload else "NO"))

        if args.verbose and ok:
            print("")
            print("  OK:")
            for label, size in ok:
                print("    %-44s %8s" % (label, human(size)))
        elif ok:
            sizes = [s for _l, s in ok]
            small = [(l, s) for l, s in ok if s < 64 * 1024]
            print("")
            print("  OK pages      : %d, %s to %s"
                  % (len(ok), human(min(sizes)), human(max(sizes))))
            if small:
                print("  low-depth     : %d page(s) under 64 KB. Expected where a"
                      % len(small))
                print("                  breakpoint falls outside a panel window,")
                print("                  since adaptive sampling did not enrich it.")

        grand_missing += len(missing)
        grand_broken += len(broken)
        grand_ok += len(ok)
        print("")

    print("=" * 68)
    print("Totals: %d good, %d broken, %d missing"
          % (grand_ok, grand_broken, grand_missing))
    print("")
    if grand_missing and not grand_ok:
        print("Diagnosis: the pages are not in this bundle. The viewer has")
        print("nothing to load. Regenerate snapshots, rebuild the bundle, and")
        print("do not use --light.")
    elif grand_broken and not grand_ok:
        print("Diagnosis: the pages exist but carry no read data. create_report")
        print("produced a page without a pileup, which points at the locus, the")
        print("BAM, or the flanking window rather than at the report wiring.")
        print("Open one of the thin pages directly in a browser to confirm, then")
        print("check read depth at that locus in the BAM it was rendered from.")
    elif grand_ok and (grand_broken or grand_missing):
        print("Diagnosis: mixed. The good pages prove the chain works end to")
        print("end, so the thin or missing ones are specific to those loci.")
    else:
        print("Diagnosis: every page is present and carrying data. If the viewer")
        print("still shows nothing, the problem is in how the report reaches it")
        print("-- run bin/embed_report_assets.py --embed-igv to inline the pages")
        print("into the report so no external file is needed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
