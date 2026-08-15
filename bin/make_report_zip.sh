#!/usr/bin/env bash
#
# make_report_zip.sh
#
# Package a built report bundle into a zip for release to reporting staff.
#
# The archive is written to the directory the command was run from. Pass --out
# to put it somewhere else.
#
# Usage:
#   tools/make_report_zip.sh <bundle_dir> [--light] [--force] [--out <dir>]
#
# Examples:
#   tools/make_report_zip.sh report_v7_20260713_24h --light
#   tools/make_report_zip.sh report_v7_20260713_24h --out /goast/reports
#
# --light omits the igv/ directories. The IGV pages embed alignment slices and
# dominate the archive size; without them every tab still renders and only the
# paired breakpoint viewer is unavailable. Use --light for anything that has to
# be moved off the server, and the full archive for local archival.
#
# The zip is built with the bundle directory as its top-level entry, so it
# unpacks into a single folder rather than scattering files.

set -euo pipefail

BUNDLE=""
LIGHT=0
FORCE=0
# Default: the working directory the command was invoked from.
OUTDIR="$PWD"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --light) LIGHT=1; shift ;;
    --force) FORCE=1; shift ;;
    --out)   OUTDIR="${2:?--out needs a directory}"; shift 2 ;;
    -h|--help)
      sed -n '2,25p' "$0" | sed 's/^# \{0,1\}//'
      exit 0 ;;
    *)
      if [[ -z "$BUNDLE" ]]; then BUNDLE="$1"; shift
      else echo "ERROR: unexpected argument: $1" >&2; exit 1; fi ;;
  esac
done

if [[ -z "$BUNDLE" ]]; then
  echo "Usage: make_report_zip.sh <bundle_dir> [--light] [--out <dir>]" >&2
  exit 1
fi

BUNDLE="${BUNDLE%/}"
if [[ ! -d "$BUNDLE" ]]; then
  echo "ERROR: bundle directory not found: $BUNDLE" >&2
  exit 1
fi
if ! command -v zip > /dev/null; then
  echo "ERROR: zip is not installed" >&2
  exit 1
fi

mkdir -p "$OUTDIR"
# Resolved to an absolute path before the cd below, so a relative --out is
# interpreted against the invocation directory rather than the bundle's parent.
OUTDIR="$(cd "$OUTDIR" && pwd)"

NAME="$(basename "$BUNDLE")"
if [[ "$LIGHT" -eq 1 ]]; then
  ZIP="${OUTDIR}/${NAME}_light.zip"
else
  ZIP="${OUTDIR}/${NAME}.zip"
fi
rm -f "$ZIP"

# Zip from the parent so the archive contains <bundle>/... at its root.
PARENT="$(cd "$(dirname "$BUNDLE")" && pwd)"

echo "Bundle : ${PARENT}/${NAME}"
echo "Archive: ${ZIP}"
if [[ "$LIGHT" -eq 1 ]]; then
  echo "Mode   : light (igv/ excluded)"
else
  echo "Mode   : full"
fi
echo ""

# --light removes the IGV tree, but the reports were built against a bundle
# that had it, so their IGV buttons stay live and open empty panes. Refuse by
# default rather than ship a report with controls that lead nowhere.
if [[ "$LIGHT" -eq 1 && "$FORCE" -eq 0 ]]; then
  if grep -rlq 'tx-igv-btn' "$BUNDLE" --include='*_report.html' 2>/dev/null; then
    echo "ERROR: these reports contain live IGV buttons, and --light removes the" >&2
    echo "       pages they open. Readers would get empty panes with no" >&2
    echo "       explanation. Options:" >&2
    echo "" >&2
    echo "  1. Drop --light. With rearrangement-only selection the IGV tree is" >&2
    echo "     a fraction of its former size and the full archive is portable." >&2
    echo "  2. Rebuild the dashboard against a bundle with no igv/ directory," >&2
    echo "     so the buttons are never rendered, then zip that." >&2
    echo "  3. Pass --force if the recipient has been told IGV is omitted." >&2
    exit 1
  fi
fi

cd "$PARENT"
# Excluded from every archive: pre-embed backups (a stale unstyled duplicate
# of each report), the intermediate tarball, and editor/OS debris.
COMMON_EXCLUDES=( -x "*.preembed" "*.tar.gz" "*.bak" "*.bak_*" "*/.DS_Store" "*/Thumbs.db" )

if [[ "$LIGHT" -eq 1 ]]; then
  zip -r -q "$ZIP" "$NAME" "${COMMON_EXCLUDES[@]}" "${NAME}/*/igv/*" "${NAME}/igv/*"
else
  zip -r -q "$ZIP" "$NAME" "${COMMON_EXCLUDES[@]}"
fi

# Windows Explorer still refuses paths beyond 260 characters when extracting,
# and these bundles nest sample/igv/translocations/<event>.A.html. Warn before
# the archive reaches a machine that cannot unpack it.
LONGEST=$(unzip -l "$ZIP" | awk '{print $4}' | awk '{ print length }' | sort -rn | head -1)
BUDGET=$((260 - 40))
if [[ -n "$LONGEST" && "$LONGEST" -gt "$BUDGET" ]]; then
  echo "WARNING: longest internal path is ${LONGEST} characters." >&2
  echo "         Windows Explorer fails past 260 including the extraction folder." >&2
  echo "         Extract near the drive root, or use --light." >&2
  echo "" >&2
fi

echo "Done."
ls -lh "$ZIP"
echo ""
echo "Entries: $(unzip -l "$ZIP" | tail -1 | awk '{print $2}')"
echo ""
echo "Top level:"
unzip -l "$ZIP" | awk '{print $4}' | grep -E "^${NAME}/[^/]+/?$" | sort -u | head -20

if [[ "$LIGHT" -eq 1 ]]; then
  echo ""
  echo "Note: the paired-breakpoint IGV viewer will report missing pages in"
  echo "      this archive. Every other tab is complete."
fi
