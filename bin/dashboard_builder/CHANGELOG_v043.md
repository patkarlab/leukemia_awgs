# dashboard_builder v0.4.3-cnv-visual-report

Release date: 2026-05-22
Previous: v0.4.2-reporting-trim

## Headline change: visual CNV reporting workflow

The CNV tab and the Reporting tab now form a complete workflow for assembling
a print-ready CNV section in a clinical report. The pathologist browses the
plot galleries in the CNV tab, ticks "Include" on the plots that show real
events, writes an interpretation under each one in the Reporting tab, and
prints to PDF.

## What changed

### parsers/cnv.py

**Bug fix that should have its own release**: the genome-wide scatter PNG
and diagram PDF were never being rendered in real builds. The parser was
looking for `<sample>_final-scatter.png` (underscore separator) but cnvkit
emits `<sample>.final-scatter.png` (dot separator). Same for the diagram.
The parser now globs both forms and recognises whichever exists.

**New plot directories surfaced**:

| Directory | What it contains |
|---|---|
| `cnvkit_plots/overview/`        | Whole-panel context plots (genome-wide, gene-summary heatmap) |
| `cnvkit_plots/combined/`        | Per-chromosome multi-panel views (23 PNGs) |
| `cnvkit_plots/per_chromosome/`  | Focused per-chromosome views (already supported, now properly parsed) |
| `cnvkit_plots/per_gene/`        | Per-gene plots, now split into priority + other |

**Priority gene list**, hardcoded in the parser as
`PRIORITY_GENES = ["KMT2A", "IKZF1", "CDKN2A", "CDKN2B", "PTPN11", "ETV6", "TET2", "TP53"]`.
The per-gene gallery presents these at the top in this order; everything
else falls into a collapsed "All other panel genes" `<details>` block. If
the pipeline doesn't emit a plot for one of the priority genes (e.g.
PTPN11 isn't currently produced), it's silently skipped and a one-line
note appears in the gallery header.

### CNV tab (in sample_report.html.j2)

The phase-3-stubbed "Files found: N per-chromosome, M per-gene" alert is
gone. The tab now renders, in order:

1. The existing **clinical CNV calls** DataTable (pre-existing, now
   actually accompanied by the plots underneath).
2. **Genome-wide overview** -- scatter + overview/ PNGs.
3. **Per-chromosome combined views** -- combined/ PNGs (23 per typical
   panel run).
4. **Per-chromosome focused views** -- per_chromosome/ PNGs.
5. **Clinically important genes** -- the 8 priority genes from
   PRIORITY_GENES.
6. **All other panel genes** -- collapsed `<details>`, opens on click.
7. **Genome-wide diagram (PDF)** -- shown as an iframe but explicitly
   not part of the plot-selection mechanism (PDFs don't embed in
   print-to-PDF of the Reporting tab).

Each plot is rendered as a Bootstrap card with the image, the label, an
**"Include in report" checkbox**, and a link to open the full-size PNG in
a new tab. Selecting the checkbox adds the plot to the Reporting tab.

A live count badge at the top right of the CNV tab shows how many plots
are currently selected for the report.

### Reporting tab

A new **"CNV findings"** section appears below the variants table when
any CNV plots are selected. Each selected plot renders as:

- The image at full-tab width, capped at 480px height with object-fit.
- A "Remove" button to deselect from the report.
- An editable **interpretation textarea** below the image. The value
  persists in localStorage under `tspipe-cnv-caption:<sample>:<plot_id>`,
  so it survives page reloads and rebuilds. Designed to be filled with
  a one-sentence clinical interpretation like
  *"Gain at 1q21 spanning RIT1 / BRINP3 consistent with chr1q duplication"*.

A **"Clear CNV selections"** button on the section header clears all CNV
selections + their captions in one step. The existing **"Clear all"**
button now wipes both variants (with their tier values) and CNVs (with
their captions), with a confirmation dialog that mentions both.

The **"Copy as TSV"** export now emits two sections separated by a blank
line:
```
# Variants
Gene\tGenomic variant\t...\tAMP/ASCO/CAP Tier (somatic)
DNMT3A\tchr2:25234373:C>T\t...\tTier II

# CNV findings
Label\tImage path\tInterpretation
KMT2A\tcnvkit_plots/per_gene/<sample>_gene_KMT2A.png\tGain spanning ...
```

### Print-friendly Reporting tab

A new `@media print` CSS block hides the sidebar, tabs, and editing
controls when the user does Browser -> Print -> Save as PDF on the
Reporting tab. CNV cards get `page-break-inside: avoid` so an image and
its caption never get split across a page boundary. Tier inputs and
caption textareas lose their borders and backgrounds when printing so
they read as plain text.

## Files touched

- `parsers/cnv.py` -- rewritten. Backward-compat aliases preserved
  (`per_chrom_pngs`, `per_gene_pngs`) so any code still using the old
  keys keeps working.
- `templates/macros.html.j2` -- new `render_cnv_plot_card` macro.
- `templates/sample_report.html.j2` -- CNV tab + Reporting tab markup
  + initCNVGallery call.
- `assets/js/variant-browser.js` -- new `setCnvCaption`/`getCnvCaption`
  on the reporting store, new `initCNVGallery()`, rewritten
  `initReportingTab()` that handles both variants and CNVs.
- `assets/css/dashboard.css` -- CNV card styles + print rules.
- `build.py` -- version bump only.

## Smoke tests

8 new headless checks in `smoketest_v043_cnv.js` against a rendered
example with real CNV data:

1. CNV tab no longer shows the phase-3 stub alert.
2. ~166 CNV plot cards rendered (scatter + overview + combined +
   per-chrom + priority + other).
3. Equal number of "Include in report" checkboxes.
4. All 7 available priority gene cards present (PTPN11 expected
   absent because the pipeline doesn't produce a plot for it).
5. Genome-wide scatter card uses the correct `.final-scatter.png`
   path (regression test for the dot/underscore parser bug).
6. CNV clinical table source HTML has all 118 data rows.
7. Reporting tab has the new CNV findings section.
8. Reporting tab variants section structure intact (no regression).

All 6 v0.4.1 CancerVar checks continue to pass against the same
rebuilt example, confirming no regression to the variants pipeline.

## Build invocation (unchanged)

```bash
python tools/dashboard_builder/build.py /path/to/run_dir \
    --verbose \
    --annotate-genebe --genebe-user EMAIL --genebe-key KEY \
    --annotate-oncokb --oncokb-token TOKEN \
    --annotate-cancervar
```

No new CLI flags. The CNV behaviour is unconditionally enabled because
the artifacts are already present in every run.

## Pending from prior handoffs (unchanged)

- Real-browser end-to-end test of OncoKB with a live token.
- Filtered variants externalisation (768 KB inline JSON -> sidecar).
- Exon column from VEP (separately, if the pipeline ever emits it).
- End-to-end test under `python3 -m http.server`.

## Specifically for the next session

- **Are PTPN11 per-gene plots intentionally omitted by the pipeline?**
  If the pipeline should be producing one, that's a `cnv_plots.py`
  fix on the pipeline side. The dashboard handles its absence
  gracefully.
- **Is the priority gene list canonical?** I hardcoded the 8 you
  named. If the panel evolves (e.g. adds RUNX1, NF1, ATM as priority
  CNV targets), the list moves in `parsers/cnv.py`.
- **Real-browser test of the full workflow**: open the CNV tab, tick
  a few plot cards across overview / combined / per-chrom / per-gene,
  switch to Reporting tab, write captions, then Browser -> Print ->
  Save as PDF and confirm the output looks like a clean clinical
  report.
