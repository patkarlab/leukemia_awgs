# dashboard_builder v0.4.5 — `--subdir` option for nested per-sample layouts

## What changed

A new `--subdir SUBDIR` flag on `build.py` causes the builder to treat
`<sample>/<SUBDIR>/` (instead of `<sample>/` directly) as each sample's
working directory. With `--subdir clinical`, the builder:

- Reads every per-sample file (`<sample>_hsmetrics.txt`,
  `<sample>_somaticseq_clinical_final.tsv`, `<sample>_igv_report.html`,
  `cnvkit_plots/...`, `cnv_consensus/...`, etc.) from `<sample>/clinical/`.
- Writes the per-sample report HTML to
  `<sample>/clinical/<sample>_report.html`.
- Writes annotation cache JSONs (when `--annotate-*` flags are used) into
  `<sample>/clinical/` alongside their source files.
- Applies the idempotent IGV hash-router patch to the IGV report at its
  actual location inside `<sample>/clinical/`.

Default (`--subdir ""`) is unchanged from v0.4.4: the builder operates
directly on `<sample>/` with no intermediate path component.

This was motivated by the nf-core port of `targeted-seq-pipeline`, which
publishes each sample's user-facing outputs under
`<run_dir>/<sample>/clinical/`. The flag lets the same builder run
unchanged against both layouts.

## Templates touched

- `templates/cohort_index.html.j2` — the two `<a href>` link targets now
  include `{{ ctx.report_subdir }}` between the sample name and the
  report filename. When `--subdir` is empty, `report_subdir` is `""` and
  the links are byte-identical to v0.4.4.

- `templates/sample_report.html.j2` — the "Back to cohort" link now uses
  `{{ cohort_link | default('../cohort_index.html') }}`. The default
  reproduces the v0.4.4 behavior; build.py supplies the correct
  `../../cohort_index.html` (or deeper) when a subdir is set.

## build.py touched

- `collect_sample_context()` gained a `subdir=""` keyword argument and
  uses a local `effective_dir = sample_dir / subdir` for every file
  lookup (HsMetrics, variants, FLT3, coverage, CNV, IGV, fastp,
  existing-dashboard, listing, and all annotation caches).
- `build()` gained `subdir=""`, validates it (no `..`, no absolute,
  no leading/trailing `/`), computes the relative `assets_prefix` and
  `cohort_link` for the report depth automatically, and threads the
  subdir-aware report URL into the cohort template via `report_subdir`.
- CLI gained `--subdir` (default `""`).
- `BUILDER_VERSION` bumped to `0.4.5-subdir`.

## Parsers, macros, assets

No changes. All eleven parsers continue to take `sample_dir` as input;
build.py passes them `effective_dir` instead when `--subdir` is set.

## Backward compatibility

Every existing invocation without `--subdir` produces byte-identical
output to v0.4.4 (verified with a small synthetic two-sample fixture in
both flat and `clinical/`-nested modes).

## Validation done

Two-sample synthetic fixture, both layouts:

| Layout                   | Cohort URLs                                    | sample assets_prefix | sample → cohort link        |
|--------------------------|------------------------------------------------|----------------------|------------------------------|
| `<sample>/`              | `./SAMPLE_A/SAMPLE_A_report.html`              | `../assets`          | `../cohort_index.html`       |
| `<sample>/clinical/`     | `./SAMPLE_A/clinical/SAMPLE_A_report.html`     | `../../assets`       | `../../cohort_index.html`    |

Both produce a valid cohort_index.html and per-sample report at the
expected location with the expected asset and back-link depths.

## nf-core integration

Two new processes, both in `modules/local/dashboard.nf`:

- `DASHBOARD_STAGE` (per-sample, label `process_single`): rolls a sample's
  clinical-bound files into a single `<sample>_clinical/` directory.
  Symlinks every input except `*_igv_report.html`, which is copied so
  the in-place hash-router patch lands on a module-owned file.

- `DASHBOARD_BUILD` (once per run, label `process_low`): collects all
  staged directories, renames each into the `<sample>/clinical/`
  layout, runs `build.py . --subdir clinical`, and publishes
  `cohort_index.html`, `assets/`, and each `<sample>_report.html` (plus
  the patched `<sample>_igv_report.html`) into
  `${params.outdir}/${meta.run_id}/`.

Container: `quay.io/biocontainers/multiqc:1.25--pyhdfd78af_0` for
DASHBOARD_BUILD (ships Python 3 + jinja2 + pandas). Pin to whatever tag
your MULTIQC process already uses to avoid an additional image pull.
