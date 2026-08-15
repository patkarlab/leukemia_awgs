#!/usr/bin/env bash
#
# build_report_bundle.sh
#
# Assemble a per-sample report bundle from a results directory. The bundle is
# what the dashboard builder runs against, and what gets archived or handed
# over, so it holds real files rather than symlinks into the work tree.
#
# Layout produced:
#   <bundle>/<sample>/snv/            clinical and filtered variant tables
#   <bundle>/<sample>/translocations/ merged and annotated SV tables
#   <bundle>/<sample>/cnv/            ichorCNA figure and fit parameters
#   <bundle>/<sample>/qc/             on-target QC plots and tables
#   <bundle>/<sample>/baf_loh/        per-sample BAF/LOH figures
#   <bundle>/<sample>/igv/            IGV snapshot pages and manifest
#   <bundle>/baf_loh/                 cohort BAF screen tables
#   <bundle>/filter_summary.tsv       cohort SNV filter summary
#
# Changes from the previous version: baf_loh/ and igv/ are now collected. Both
# were being produced by the pipeline and left behind in the results tree, so
# the dashboard rendered an empty BAF/LOH tab and had no IGV pages to link to.
#
# Usage:
#   bin/build_report_bundle.sh <results_dir> [bundle_name]
# Example:
#   bin/build_report_bundle.sh results_v7_20260713_24h report_v7_20260713_24h

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Locate the column-alias helper.
#
# Run by hand, it sits beside this script. Run as a Nextflow process, this
# script is staged into a task directory on its own and the sibling is not
# there, but bin/ is on PATH. Both are tried, and if neither resolves the
# failure is loud: the previous version fell back to a plain copy and reported
# success, producing variant tables with no alias columns and a report whose
# fields silently came back empty.
ALIAS_SCRIPT="${ALIAS_SCRIPT:-}"
if [[ -z "$ALIAS_SCRIPT" ]]; then
    if [[ -f "${SCRIPT_DIR}/alias_variant_table.py" ]]; then
        ALIAS_SCRIPT="${SCRIPT_DIR}/alias_variant_table.py"
    elif command -v alias_variant_table.py > /dev/null 2>&1; then
        ALIAS_SCRIPT="$(command -v alias_variant_table.py)"
    fi
fi
if [[ -z "$ALIAS_SCRIPT" ]]; then
    echo "ERROR: alias_variant_table.py not found beside this script or on PATH." >&2
    echo "       Without it the variant browser cannot resolve its columns and" >&2
    echo "       the Variants tabs render empty. Set ALIAS_SCRIPT to override." >&2
    exit 1
fi

RESULTS="${1:?Usage: build_report_bundle.sh <results_dir> [bundle_name]}"
BUNDLE="${2:-report_$(basename "$RESULTS")}"

if [[ ! -d "$RESULTS" ]]; then
  echo "ERROR: results dir not found: $RESULTS" >&2
  exit 1
fi

# Derive the sample list from the annotated TSVs (most reliable anchor).
mapfile -t SAMPLES < <(
  find "$RESULTS" \( -name '*.leukemia_annotated.tsv' -o -name '*.mm_annotated.tsv' \) -printf '%f\n' 2>/dev/null \
    | sed 's/\.mm_annotated\.tsv$//' | sort -u
)
if [[ ${#SAMPLES[@]} -eq 0 ]]; then
  mapfile -t SAMPLES < <(
    find "$RESULTS" -name '*clinical.tsv' -printf '%f\n' 2>/dev/null \
      | sed -E 's/(\.somatic_candidates)?(\.withAD)?\.v6_clinical\.tsv$//' | sort -u
  )
fi
if [[ ${#SAMPLES[@]} -eq 0 ]]; then
  echo "ERROR: no samples found under $RESULTS" >&2
  exit 1
fi

echo "Samples detected: ${SAMPLES[*]}"
rm -rf "$BUNDLE"
mkdir -p "$BUNDLE"

# copy_first <destdir> <destname> <sample> <find-args...>
copy_first() {
  local destdir="$1"; shift
  local destname="$1"; shift
  local sample="$1"; shift
  local src
  src=$(find "$RESULTS" -type f "$@" 2>/dev/null | grep -F "$sample" | head -1 || true)
  if [[ -n "$src" && -f "$src" ]]; then
    mkdir -p "$destdir"
    cp -L "$src" "$destdir/$destname"
    echo "  + $destname"
  else
    echo "  - (missing) $destname" >&2
  fi
}

# copy_all <destdir> <sample> <find-args...>
# Copies every match, preserving the source filename. Used where the count is
# not known ahead of time, such as IGV breakpoint pages.
copy_all() {
  local destdir="$1"; shift
  local sample="$1"; shift
  local count=0
  while IFS= read -r src; do
    [[ -f "$src" ]] || continue
    mkdir -p "$destdir"
    cp -L "$src" "$destdir/$(basename "$src")"
    count=$((count + 1))
  done < <(find "$RESULTS" -type f "$@" 2>/dev/null | grep -F "$sample" || true)
  if [[ $count -gt 0 ]]; then
    echo "  + ${count} file(s) -> $(basename "$destdir")/"
  fi
  return 0
}

# copy_keep_name <destdir> <sample> <find-args...>
# Copies the first match under its original basename, skipping the copy if a
# file of that name is already present.
copy_keep_name() {
  local destdir="$1"; shift
  local sample="$1"; shift
  local src
  src=$(find "$RESULTS" -type f "$@" 2>/dev/null | grep -F "$sample" | head -1 || true)
  if [[ -n "$src" && -f "$src" ]]; then
    mkdir -p "$destdir"
    local base
    base=$(basename "$src")
    if [[ ! -f "$destdir/$base" ]]; then
      cp -L "$src" "$destdir/$base"
      echo "  + $base (original name)"
    fi
  fi
}

for s in "${SAMPLES[@]}"; do
  echo "== $s =="
  d="$BUNDLE/$s"

  # SNV. One file per class, written under the name the dashboard builder
  # discovers, with alias columns added so its variant browser can read them.
  #
  # The browser resolves columns by exact name and was written against a table
  # using capitalised names, while this pipeline's filter emits mostly
  # lowercase ones. Only REF_COUNT, ALT_COUNT and Filter match in both, which
  # is why those three were the only fields that rendered. alias_variant_table
  # appends a capitalised duplicate of each column, preserving the originals,
  # so both readers are satisfied from a single file.
  #
  # Names come from build.py's own discovery patterns. Override with
  # SNV_ALIAS_CLINICAL / SNV_ALIAS_FILTERED.
  alias_clin="${SNV_ALIAS_CLINICAL:-${s}_somaticseq_clinical_final.tsv}"
  alias_filt="${SNV_ALIAS_FILTERED:-${s}_somaticseq_filtered.tsv}"

  src_clin=$(find "$RESULTS" -type f -name '*clinical.tsv' 2>/dev/null | grep -F "$s" | head -1 || true)
  src_filt=$(find "$RESULTS" -type f -name '*filtered.tsv'  2>/dev/null | grep -F "$s" | head -1 || true)

  for pair in "clin:${src_clin}:${alias_clin}" "filt:${src_filt}:${alias_filt}"; do
    kind="${pair%%:*}"; rest="${pair#*:}"
    src="${rest%%:*}"; dest="${rest#*:}"
    if [[ -n "$src" && -f "$src" ]]; then
      mkdir -p "$d/snv"
      # SNV_EXTRA_ALIASES is a space-separated list of NAME=SOURCE pairs,
      # for a consumer that reads a column name the built-in set misses.
      #   SNV_EXTRA_ALIASES="start=pos Locus=chrom" bin/build_report_bundle.sh ...
      extra_args=()
      for a in ${SNV_EXTRA_ALIASES:-}; do extra_args+=( --extra-alias "$a" ); done
      # Not guarded by a fallback on purpose. A plain copy here produces a
      # bundle that looks complete and a report with empty variant fields,
      # which is harder to notice than a failed run.
      python3 "$ALIAS_SCRIPT" "$src" "$d/snv/${dest}" "${extra_args[@]+"${extra_args[@]}"}"
      cp -L "$d/snv/${dest}" "$d/${dest}"
      echo "  + ${dest} (aliased columns, snv/ and sample root)"
    else
      echo "  - (missing) ${kind} variant table" >&2
    fi
  done

  # Focal duplications. The dashboard's FLT3 tab reads a consensus table with
  # its own schema; the AL callers emit one file per gene with a different one.
  # Translated here rather than in the parser, so the dashboard stays byte-for-
  # byte the myeloma version and the AL side owns the mapping.
  dup_files=$(find "$RESULTS" -type f \( -name "*.tandem_dup.tsv" -o -name "*.intragenic_dup.tsv" \) 2>/dev/null | grep -F "$s" || true)
  if [[ -n "$dup_files" ]]; then
    mkdir -p "$d"
    python3 - "$d/${s}_flt3_consensus.tsv" $dup_files <<'PYDUP'
import csv, sys
out, srcs = sys.argv[1], sys.argv[2:]
cols = ["sample","status","n_tools","tools","length_bp","length_range","pos_hg38",
        "vaf_pct_min","vaf_pct_max","vaf_pct_mean","ar_min","ar_max","ar_mean",
        "hgvsc","hgvsp","domain","inserted_seq","raw_calls"]
rows = []
for src in srcs:
    for r in csv.DictReader(open(src), delimiter="\t"):
        conf = (r.get("confidence") or "").lower()
        # status drives the parser's positive count. Only a graded call counts;
        # everything else is recorded as negative so it still shows in the table
        # without inflating the headline.
        status = r.get("label") or r.get("gene") or "duplication"
        if conf and conf not in ("high", "moderate"):
            status = "negative"
        ar = r.get("allelic_ratio") or ""
        vaf = f"{float(ar)*100:.1f}" if ar else ""
        length = r.get("dup_length_bp") or r.get("duplication_span_bp") or ""
        rows.append({
            "sample": r.get("sample",""), "status": status,
            "n_tools": "1", "tools": r.get("callers","call_tandem_dup"),
            "length_bp": length, "length_range": length,
            "pos_hg38": f"{r.get('chrom','')}:{r.get('ref_pos') or r.get('acceptor_pos','')}",
            "vaf_pct_min": vaf, "vaf_pct_max": vaf, "vaf_pct_mean": vaf,
            "ar_min": ar, "ar_max": ar, "ar_mean": ar,
            "hgvsc": "", "hgvsp": "",
            "domain": r.get("hotspot",""), "inserted_seq": "",
            "raw_calls": "; ".join(f"{k}={v}" for k, v in r.items()
                                   if k in ("gene","in_frame","support_reads",
                                            "spanning_reads","tandem_identity",
                                            "confidence","split_read_support",
                                            "n_signals","note") and v),
        })
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
    w.writeheader(); w.writerows(rows)
print(f"  + {out.split('/')[-1]} ({len(rows)} duplication call(s))")
PYDUP
  fi

  # Translocations
  copy_first "$d/translocations" "${s}.leukemia_annotated.tsv"   "$s" \\
      \( -name '*.leukemia_annotated.tsv' -o -name '*.mm_annotated.tsv' \)
  # The AL annotator names the dictionary hit known_pair; the variant browser
  # resolves columns by exact name and looks for known_mm_pair. Append the
  # alias rather than renaming, so both readers are satisfied from one file and
  # the AL column stays authoritative.
  if [[ -f "$d/translocations/${s}.leukemia_annotated.tsv" ]]; then
    python3 - "$d/translocations/${s}.leukemia_annotated.tsv" <<'PYALIAS'
import csv, sys
p = sys.argv[1]
rows = list(csv.DictReader(open(p), delimiter="\t"))
if rows and "known_pair" in rows[0] and "known_mm_pair" not in rows[0]:
    cols = list(rows[0].keys())
    cols.insert(cols.index("known_pair") + 1, "known_mm_pair")
    for r in rows:
        r["known_mm_pair"] = r.get("known_pair", "")
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)
PYALIAS
  fi
  copy_first "$d/translocations" "${s}.translocations.tsv" "$s" -name '*.translocations.tsv'

  # CNV: the all-solutions figure and the fit parameters only.
  copy_first "$d/cnv" "${s}.ichor_all_sols.pdf" "$s" -path '*ichor*' -name '*all_sols*.pdf'
  copy_first "$d/cnv" "${s}.ichor_params.txt" "$s" -path '*ichor*' -name '*params.txt'

  # ichorCNA emits every ploidy/normal-fraction solution into one multi-page
  # PDF. A PDF page cannot be ticked, so the pages are rasterised here and the
  # CNV tab renders them as selectable plot cards. That is what lets a reviewer
  # choose the fit that is actually right for the sample: on this data the
  # top-likelihood solution is often within one log-likelihood unit of two
  # others, so the automatic pick is not authoritative.
  if [[ -f "$d/cnv/${s}.ichor_all_sols.pdf" ]]; then
    if command -v pdftoppm > /dev/null 2>&1; then
      mkdir -p "$d/cnv/all_sols"
      pdftoppm -png -r 110 "$d/cnv/${s}.ichor_all_sols.pdf" \
               "$d/cnv/all_sols/${s}.sol" 2>/dev/null || true
      n=$(find "$d/cnv/all_sols" -name "${s}.sol*.png" | wc -l)
      echo "  + cnv/all_sols/ (${n} solution page(s) as selectable images)"
    else
      echo "  - (skipped) all-solution page images: pdftoppm not on PATH" >&2
    fi
  fi

  # QC: adaptive-sampling plots and per-region coverage.
  copy_first "$d/qc" "${s}.region_coverage.tsv" "$s" -name '*.region_coverage.tsv'
  copy_first "$d/qc" "${s}.region_coverage.png" "$s" -name '*.region_coverage.png'
  copy_first "$d/qc" "${s}.readlen_hist.png"    "$s" -name '*.readlen_hist.png'
  copy_first "$d/qc" "${s}.qscore_hist.png"     "$s" -name '*.qscore_hist.png'
  copy_first "$d/qc" "${s}.readlen_qscore.tsv"  "$s" -name '*.readlen_qscore.tsv'

  # BAF / LOH: per-sample figures. The cohort tables are copied once below.
  copy_all "$d/baf_loh" "$s" -path '*baf_loh*' -name '*.png'

  # IGV: breakpoint pages, the somatic page, and the manifest that maps
  # events to pages. Directory structure is flattened per evidence class so
  # the dashboard's relative hrefs stay short.
  copy_all "$d/igv/translocations" "$s" -path '*igv*translocations*' -name '*.html'
  copy_all "$d/igv/translocations" "$s" -path '*igv*' -name '*.manifest.json'
  copy_all "$d/igv/somatic"        "$s" -path '*igv*somatic*' -name '*.html'
  copy_all "$d/igv/focal_dup"      "$s" -path '*igv*focal_dup*' -name '*.html'

  # The IGV tab loads exactly one page, resolved by build.py as
  #   effective_dir / f"{sample}_igv_report.html"
  # so the somatic page is placed there. The builder then extracts a row
  # lookup from it and injects a hash router in place, which is what gives the
  # clinical variant cards their working IGV links. Override with IGV_ALIAS
  # only if that expression changes.
  som=$(find "$RESULTS" -type f -ipath '*igv*somatic*' -name '*.html' 2>/dev/null | grep -F "$s" | head -1 || true)
  if [[ -n "$som" && -f "$som" ]]; then
    igv_dest="${IGV_ALIAS:-${s}_igv_report.html}"
    cp -L "$som" "$d/${igv_dest}"
    echo "  + ${igv_dest} (IGV tab source)"
  else
    echo "  - (missing) somatic IGV page; IGV tab and variant IGV links stay empty" >&2
  fi
done

# Cohort-level artefacts (single copies, not per-sample).
summary=$(find "$RESULTS" -name 'v6_filter_summary.tsv' 2>/dev/null | head -1 || true)
if [[ -n "$summary" ]]; then
  cp -L "$summary" "$BUNDLE/filter_summary.tsv"
  echo "+ cohort filter_summary.tsv"
fi

# The pipeline publishes BAF/LOH under more than one path, and not all of them
# hold the cohort tables. Pick the directory that actually contains a
# cohort*.tsv rather than whichever find reaches first.
baf_dir=""
while IFS= read -r candidate; do
  if compgen -G "${candidate}/cohort*.tsv" > /dev/null; then
    baf_dir="$candidate"
    break
  fi
done < <(find "$RESULTS" -type d -name 'baf_loh' 2>/dev/null | sort)

if [[ -z "$baf_dir" ]]; then
  # Fall back to any baf_loh directory so figures are still collected.
  baf_dir=$(find "$RESULTS" -type d -name 'baf_loh' 2>/dev/null | head -1 || true)
  [[ -n "$baf_dir" ]] && echo "  note: no cohort*.tsv found; using $baf_dir" >&2
fi

if [[ -n "$baf_dir" ]]; then
  mkdir -p "$BUNDLE/baf_loh"
  find "$baf_dir" -maxdepth 1 -type f -name '*.tsv' -exec cp -L {} "$BUNDLE/baf_loh/" \;

  # Cohort figures. The pipeline module writes them to baf_cn_figures/ while
  # the standalone script writes figures/, so any immediate subdirectory
  # holding PNGs is taken and normalised to figures/ in the bundle, which is
  # where the dashboard looks. Matching only figures/ silently skipped the
  # cohort plots on pipeline runs.
  fig_count=0
  for figdir in "$baf_dir"/*/; do
    [[ -d "$figdir" ]] || continue
    if compgen -G "${figdir}*.png" > /dev/null 2>&1; then
      mkdir -p "$BUNDLE/baf_loh/figures"
      cp -L "${figdir}"*.png "$BUNDLE/baf_loh/figures/" 2>/dev/null || true
      fig_count=$((fig_count + 1))
      echo "  cohort figures from $(basename "$figdir")" >&2
    fi
  done
  if [[ "$fig_count" -eq 0 ]]; then
    echo "  WARNING: no cohort BAF/LOH figures found under ${baf_dir}" >&2
  fi
  n_baf=$(find "$BUNDLE/baf_loh" -type f | wc -l)
  echo "+ cohort baf_loh/ from ${baf_dir} (${n_baf} files)"
  if [[ "$n_baf" -eq 0 ]]; then
    echo "  WARNING: baf_loh directory found but empty; the BAF/LOH tab will render empty" >&2
  fi
else
  echo "- (missing) cohort baf_loh/" >&2
fi

# The tarball is off by default: tools/make_report_zip.sh produces the
# distributable archive, and writing both leaves two multi-hundred-megabyte
# copies of the same tree on disk. Set BUNDLE_TAR=1 if the tar is wanted.
if [[ "${BUNDLE_TAR:-0}" == "1" ]]; then
  tar czf "${BUNDLE}.tar.gz" "$BUNDLE"
  echo ""
  echo "Tarball:     ${BUNDLE}.tar.gz"
  du -sh "${BUNDLE}.tar.gz"
fi
# ---------------------------------------------------------------------------
# Alignments.
#
# The reports name a locus and a read count; reviewing either means opening the
# alignment. Included by default so the bundle is self-sufficient rather than a
# set of assertions the reader has to take on trust.
#
# Sliced, not whole-genome. A 20 h adaptive-sampling run is around 30 GB per
# sample per reference, so three samples across both references is 180 GB and
# the bundle stops being something anyone can receive. Almost all of that is
# rejected off-target reads at about 1x that nothing called from and nobody
# reviews. Sliced to the panel the same alignment is a few hundred megabytes
# and still carries every read behind every call.
#
# Both references travel because the two tracks are not interchangeable:
# fusions and intragenic duplications were called on T2T, variants and tandem
# duplications on hg38. Opening a call against the wrong assembly gives
# coordinates that look plausible and are wrong.
#
# BUNDLE_NO_BAMS=1 to skip. BUNDLE_BAM_REGIONS to narrow further, e.g.
#   BUNDLE_BAM_REGIONS="chr13:28030000-28040000 chr17:44200000-44220000"
if [[ -z "${BUNDLE_NO_BAMS:-}" ]]; then
  if ! command -v samtools > /dev/null 2>&1; then
    echo "WARNING: samtools not on PATH; alignments omitted from the bundle." >&2
  else
    for s in "${SAMPLES[@]}"; do
      for bam in $(find "$RESULTS" -type f -name "${s}.t2t.bam" -o -type f -name "${s}.hg38.bam" 2>/dev/null); do
        base="$(basename "$bam")"
        case "$base" in
          *.t2t.bam)  ref_bed="${PANEL_BED_T2T:-${SCRIPT_DIR}/../assets/AML_panel_t2t_chr.bed}" ;;
          *.hg38.bam) ref_bed="${PANEL_BED_HG38:-${SCRIPT_DIR}/../assets/AML_panel_hg38.bed}" ;;
          *) continue ;;
        esac
        mkdir -p "$BUNDLE/$s/bam"
        dest="$BUNDLE/$s/bam/$base"
        if [[ -n "${BUNDLE_BAM_REGIONS:-}" && "$base" == *.hg38.bam ]]; then
          # Region strings are assembly-specific, so they are applied to hg38
          # only. A T2T alignment falls back to its own BED rather than being
          # sliced at hg38 coordinates.
          # shellcheck disable=SC2086
          samtools view -b -o "$dest" "$bam" ${BUNDLE_BAM_REGIONS}
        elif [[ -f "$ref_bed" ]]; then
          samtools view -b -L "$ref_bed" -o "$dest" "$bam"
        else
          echo "  - (missing) panel BED for $base: $ref_bed" >&2
          continue
        fi
        samtools index "$dest"
        echo "  + $s/bam/$base ($(du -h "$dest" | cut -f1), sliced)"
      done
    done
  fi
fi

echo ""
echo "Bundle tree: $BUNDLE/"
du -sh "$BUNDLE"
echo "Archive it with: bin/make_report_zip.sh $BUNDLE"
echo ""
echo "Per-sample contents:"
find "$BUNDLE" -mindepth 2 -maxdepth 2 -type d | sed "s|^$BUNDLE/||" | sort | uniq -c | sort -rn | head -20
