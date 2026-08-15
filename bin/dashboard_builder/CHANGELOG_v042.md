# dashboard_builder v0.4.2-reporting-trim

Release date: 2026-05-22
Previous: v0.4.1-cancervar

## What's changed

### Reporting tab: Exon column removed

The Exon column has been dropped from the Reporting tab. It was reserving
space for a field that the pipeline does not currently emit in the clinical
TSV -- the cells were rendering blank for every selected variant, taking up
horizontal space without carrying information.

The Reporting tab column order is now:

  Gene | Genomic variant | Protein variant | Transcript Variant
       | COSMIC database reference | VAF (%) | AMP/ASCO/CAP Tier (somatic)

The "AMP/ASCO/CAP Tier (somatic)" column remains the only editable column.
When CancerVar annotation is enabled at build time, the tier is pre-filled
on first selection with CancerVar's tier classification (Tier I/II/III/IV),
and can be edited or cleared like any manual value (behavior unchanged
from v0.4.1).

The Exon snapshot field is still captured at toggle-time (in the localStorage
selection record) for forward compatibility, in case the pipeline starts
emitting an EXON column upstream -- but it is no longer rendered or copied
to TSV.

### TSV export format change

The "Copy as TSV" output drops the Exon column accordingly. If you have
downstream consumers of the pasted TSV that expected an 8-column layout
(Gene, Genomic, Protein, Transcript, **Exon**, COSMIC, VAF, Tier), they
will need to expect 7 columns now.

## No other changes

Everything else from v0.4.1 carries forward unchanged:

- CancerVar annotation pass (`--annotate-cancervar`) and detail block
- OncoKB annotation pass (`--annotate-oncokb --oncokb-token TOKEN`)
- GeneBe annotation pass and detail block
- Tier pre-fill from CancerVar on first variant selection
- IGV hash-router injection
- Lazy iframe loading
- Filtered tab deferred initialization

## Files touched

Only `templates/sample_report.html.j2` (four lines: header `<th>`, COLUMNS
array, rowToTsv, row render `<td>`) and the version string in `build.py`.
No parser or JS-asset changes.
