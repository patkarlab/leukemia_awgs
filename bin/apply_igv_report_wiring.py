#!/usr/bin/env python3
"""
apply_igv_report_wiring.py

Wire IGV_SNAPSHOTS, EMBED_REPORT_ASSETS and REPORT_ZIP into the main workflow.

The workflow already calls REPORT_BUNDLE and DASHBOARD behind an ordering gate,
so this patch is additive: it adds the IGV stage before the gate, mixes IGV
into the gate, and extends the reporting block with the two packaging steps.
Nothing existing is replaced.

Edits
-----
workflows/<pipeline>.nf
    1. three includes, after the DASHBOARD include
    2. an IGV_SNAPSHOTS block, after the HG38_TRACK call
    3. IGV mixed into the `ready` gate, before bundle_name is computed
    4. EMBED_REPORT_ASSETS and REPORT_ZIP, after the DASHBOARD call

nextflow.config
    5. IGV and report parameters, after the hg38_mmi entry

On the join
-----------
The four channels are joined on meta so that a sample's SV table stays paired
with its own alignments. A `combine` here would cross samples, which in a
clinical report means one patient's breakpoints rendered against another's
reads.

`v6_report` is joined with `remainder: true` on purpose. A sample with no
on-panel clinical SNVs never emits it, and a plain join would drop that sample
from the IGV stage entirely, taking its translocation pages with it. With
remainder the sample survives and the clinical table arrives as an empty list,
which the process treats as "no somatic snapshots" rather than an error.

Safety
------
Every anchor is validated before anything is written. A timestamped .bak is
kept per file. Each edit carries a sentinel and is skipped if already present,
so the script is safe to re-run. --dry-run reports the plan and writes nothing.

Usage
-----
    python3 bin/apply_igv_report_wiring.py --dry-run
    python3 bin/apply_igv_report_wiring.py
    python3 bin/apply_igv_report_wiring.py --repo /goast/mm-awgs-nextflow
"""

import argparse
import glob
import os
import shutil
import sys
from datetime import datetime


INCLUDES = """include { IGV_SNAPSHOTS       } from '../modules/local/igv_snapshots.nf'
include { EMBED_REPORT_ASSETS } from '../modules/local/embed_report_assets.nf'
include { REPORT_ZIP          } from '../modules/local/report_zip.nf'
"""

IGV_BLOCK = """
    // 4b. IGV snapshots.
    //
    // Two evidence classes per sample: paired breakpoint pages for each
    // rearrangement against T2T, and one clinical SNV page against hg38. The
    // SNV page is published under the filename the dashboard builder
    // resolves, which also gives the variant cards their IGV links.
    //
    // Joined on meta so each sample's tables stay with its own alignments.
    // v6_report uses remainder because a sample with no on-panel clinical
    // SNVs never emits one; without it that sample would be dropped here and
    // lose its translocation pages too.
    if (!params.skip_igv && !params.skip_t2t_track && !params.skip_hg38_track) {
        igv_input = T2T_TRACK.out.mm_annotated_tsv
            .join(T2T_TRACK.out.t2t_bam_bai)
            .join(HG38_TRACK.out.hg38_bam_bai)
            .join(HG38_TRACK.out.v6_report, remainder: true)
            .filter { it[0] != null && it[1] != null }
            .map { meta, mm, tbam, tbai, hbam, hbai, clin ->
                tuple(meta, mm, clin ?: [], tbam, tbai, hbam, hbai)
            }

        IGV_SNAPSHOTS(
            igv_input,
            file(params.t2t_fasta),
            file(params.t2t_fai),
            file(params.hg38_fasta),
            file(params.hg38_fai)
        )
    }
"""

GATE_BLOCK = """        // IGV pages are collected by the bundle, so the bundle must not
        // start before they are published.
        if (!params.skip_igv && !params.skip_t2t_track && !params.skip_hg38_track) {
            ready = ready.mix(IGV_SNAPSHOTS.out.igv.map { it -> 'ok' })
        }

"""

PACKAGING_BLOCK = """
            // Inline every local dependency, then package. Without the embed
            // step the reports reference their stylesheets and figures by
            // relative path and break as soon as they are moved.
            if (!params.skip_report_package) {
                EMBED_REPORT_ASSETS(DASHBOARD.out.bundle)
                REPORT_ZIP(EMBED_REPORT_ASSETS.out.bundle)
            }
"""

PARAMS_BLOCK = """
    // ===== IGV snapshots =====
    skip_igv              = false
    igv_flanking          = 5000
    // sv_type values rendered. The annotated table covers the whole merged
    // callset, so without this every deletion and insertion gets a page pair.
    igv_sv_types          = 'TRA'
    // Deliberately 1: single-caller rearrangements at low read support are
    // what this panel exists to recover.
    igv_min_callers       = 1
    igv_max_events        = 200

    // ===== Report packaging =====
    skip_report_package   = false
    // Inline breakpoint pages for the first N events as data URIs. 0 keeps
    // reports small; the archive carries igv/ regardless.
    report_embed_igv      = 0
    report_zip_light      = false
    report_zip_force      = false
"""


class PatchError(Exception):
    pass


def read(path):
    with open(path, encoding="utf-8") as handle:
        return handle.read()


def backup_and_write(path, text):
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    shutil.copy2(path, "%s.bak_igv_wiring_%s" % (path, stamp))
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(text)


def after_line(text, needle, description):
    """Offset just past the end of the line containing needle."""
    index = text.find(needle)
    if index == -1:
        raise PatchError("anchor not found: %s (%r)" % (description, needle))
    end = text.find("\n", index)
    return len(text) if end == -1 else end + 1


def after_block(text, needle, description):
    """Offset just past the closing brace of the block containing needle.

    Used for the `if (...) { ... }` wrappers around the track calls, where the
    insertion point is after the whole conditional rather than after a line.
    """
    index = text.find(needle)
    if index == -1:
        raise PatchError("anchor not found: %s (%r)" % (description, needle))
    depth = 0
    seen = False
    position = index
    while position < len(text):
        char = text[position]
        if char == "{":
            depth += 1
            seen = True
        elif char == "}":
            depth -= 1
            if seen and depth <= 0:
                end = position + 1
                if text[end : end + 1] == "\n":
                    end += 1
                return end
        position += 1
    raise PatchError("unbalanced braces after anchor: %s" % description)


def find_workflow(repo):
    candidates = sorted(glob.glob(os.path.join(repo, "workflows", "*.nf")))
    hits = [p for p in candidates if "REPORT_BUNDLE" in read(p)]
    if not hits:
        raise PatchError(
            "no workflow under workflows/ calls REPORT_BUNDLE; nothing to patch"
        )
    if len(hits) > 1:
        raise PatchError(
            "more than one workflow calls REPORT_BUNDLE: %s" % ", ".join(hits)
        )
    return hits[0]


def patch_workflow(path):
    text = read(path)
    actions = []

    if "IGV_SNAPSHOTS" in text and "include { IGV_SNAPSHOTS" in text:
        actions.append("SKIP  workflow includes: already applied")
    else:
        offset = after_line(
            text, "include { DASHBOARD", "DASHBOARD include"
        )
        text = text[:offset] + INCLUDES + text[offset:]
        actions.append("EDIT  workflow: add IGV/embed/zip includes")

    if "IGV_SNAPSHOTS(" in text:
        actions.append("SKIP  workflow IGV stage: already applied")
    else:
        offset = after_block(
            text, "if (!params.skip_hg38_track)", "HG38_TRACK conditional"
        )
        text = text[:offset] + IGV_BLOCK + text[offset:]
        actions.append("EDIT  workflow: add IGV_SNAPSHOTS stage after HG38_TRACK")

    if "IGV_SNAPSHOTS.out.igv.map" in text:
        actions.append("SKIP  workflow ordering gate: already applied")
    else:
        index = text.find("bundle_name = params.report_bundle_name")
        if index == -1:
            raise PatchError("anchor not found: bundle_name assignment")
        line_start = text.rfind("\n", 0, index) + 1
        text = text[:line_start] + GATE_BLOCK + text[line_start:]
        actions.append("EDIT  workflow: mix IGV into the ordering gate")

    if "EMBED_REPORT_ASSETS(" in text:
        actions.append("SKIP  workflow packaging: already applied")
    else:
        offset = after_block(
            text, "if (!params.skip_dashboard)", "DASHBOARD conditional"
        )
        # Insert inside the conditional, immediately before its closing brace.
        closing = text.rfind("}", 0, offset)
        line_start = text.rfind("\n", 0, closing) + 1
        text = text[:line_start] + PACKAGING_BLOCK + text[line_start:]
        actions.append("EDIT  workflow: add EMBED_REPORT_ASSETS and REPORT_ZIP")

    return text, actions


def patch_config(path):
    text = read(path)
    if "skip_igv" in text:
        return text, ["SKIP  nextflow.config: parameters already present"]
    offset = after_line(text, "hg38_mmi", "hg38_mmi parameter")
    text = text[:offset] + PARAMS_BLOCK + text[offset:]
    return text, ["EDIT  nextflow.config: add IGV and report parameters"]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--repo", default=".", help="repository root")
    parser.add_argument("--dry-run", action="store_true",
                        help="report the plan, write nothing")
    args = parser.parse_args(argv)

    config_path = os.path.join(args.repo, "nextflow.config")
    if not os.path.isfile(config_path):
        print("ERROR: nextflow.config not found under %s" % args.repo,
              file=sys.stderr)
        return 2

    required = [
        os.path.join(args.repo, "modules", "local", "igv_snapshots.nf"),
        os.path.join(args.repo, "modules", "local", "embed_report_assets.nf"),
        os.path.join(args.repo, "modules", "local", "report_zip.nf"),
        os.path.join(args.repo, "bin", "igv_snapshots.py"),
        os.path.join(args.repo, "bin", "embed_report_assets.py"),
        os.path.join(args.repo, "tools", "make_report_zip.sh"),
    ]
    missing = [p for p in required if not os.path.isfile(p)]
    if missing:
        print("ERROR: install the modules and scripts first. Missing:",
              file=sys.stderr)
        for path in missing:
            print("  " + os.path.relpath(path, args.repo), file=sys.stderr)
        return 2

    try:
        workflow_path = find_workflow(args.repo)
        workflow_text, workflow_actions = patch_workflow(workflow_path)
        config_text, config_actions = patch_config(config_path)
    except PatchError as error:
        print("ERROR: %s" % error, file=sys.stderr)
        print("Nothing was written.", file=sys.stderr)
        return 2

    print("workflow: %s" % os.path.relpath(workflow_path, args.repo))
    for action in workflow_actions + config_actions:
        print(action)

    edits = [a for a in workflow_actions + config_actions if a.startswith("EDIT")]
    if not edits:
        print("\nNothing to do; all edits already applied.")
        return 0

    if args.dry_run:
        print("\nDry run: %d edit(s) planned, nothing written." % len(edits))
        return 0

    backup_and_write(workflow_path, workflow_text)
    backup_and_write(config_path, config_text)
    print("\nApplied %d edit(s). Backups written alongside each file." % len(edits))
    print("\nVerify the wiring without compute:")
    print("  nextflow run main.nf -profile docker -stub-run \\")
    print("    --sample_sheet <sheet> --outdir results_stub")
    print("  grep -E 'IGV_SNAPSHOTS|EMBED_REPORT_ASSETS|REPORT_ZIP' .nextflow.log | tail")
    return 0


if __name__ == "__main__":
    sys.exit(main())
