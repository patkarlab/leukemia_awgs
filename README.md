# leukemia_awgs

**Adaptive whole-genome sequencing (Oxford Nanopore) for acute leukaemia —
a dual-reference, multi-caller pipeline covering AML and B/T-ALL, ported from
[patkarlab/mm-awgs-nextflow](https://github.com/patkarlab/mm-awgs-nextflow).**

Status: `v0.1.0-dev`. See [Scope](#v01-scope) and [Known limitations](#known-limitations).

## Why dual reference

| Reference | Used for | Why |
|---|---|---|
| **T2T-CHM13v2.0** | SV calling, fusion annotation, intragenic duplication, on-target QC | Fusion partners cluster in regions hg38 represents poorly. KMT2A sits in a segmentally duplicated stretch of 11q23 where hg38 alt-contigs pull away the supplementary alignments that intragenic duplication calling depends on. The Ig and TR loci are fully assembled in T2T and fragmented in hg38. |
| **hg38** | SNV, indel, large CNV, phased germline, VEP annotation, tandem duplication | Mature tooling. ClairS-TO, ichorCNA, Clair3 and VEP all ship hg38-coordinate reference panels. Clinical exon numbering and the ELN 2022 risk tables assume hg38. |

Both tracks consume the same MinKNOW BAM and run in parallel.

## One pipeline, two panels

Lineage is assigned before this reflex assay fires, so an AML and an ALL panel
never share a flow cell. They do share a pipeline. Panel selection is per
sample, via an optional `panel` column in the sample sheet falling back to
`--panel`, so a batch spanning several flow cells runs in a single invocation
and each sample resolves its own BEDs.

| | AML | ALL |
|---|---|---|
| Design targets | 120 | 133 |
| T2T regions / bases | 117 / 29.94 Mb | 126 / 40.74 Mb |
| hg38 regions / bases | 117 / 29.43 Mb | 126 / 38.00 Mb |
| Fraction of T2T-CHM13v2.0 | 0.960% | 1.307% |
| Shared region labels | 54 | 54 |

Adding a third panel means adding a block to `conf/panels.config` and dropping
its BEDs in `assets/`. No workflow code changes.

## What this pipeline does

```
sample sheet (sample_id, minknow_bam, panel, ...)
    │
    ├── T2T track ──┬─► Sniffles2 + CuteSV + Severus ─► SURVIVOR merge
    │               │        └─► AL fusion annotation ─► support enrichment ─► junction merge
    │               ├─► on-target QC (per-region coverage, read length, Q)
    │               └─► intragenic duplication (split reads + insertions + exon coverage)
    │
    └── hg38 track ─┬─► ClairS-TO (somatic SNV/indel, panel-restricted)
                    ├─► ichorCNA (large CNV from off-target reads)
                    ├─► Clair3 phased germline ─► VEP ─► AL panel filter
                    └─► tandem duplication (CIGAR insertions + tandem verification)
```

## Nothing biological is hardcoded

Every clinical prior, gene list, coordinate and synonym lives in a TSV or BED
under `assets/`. No gene symbol, genomic coordinate or known mutation appears
in any `.py`, `.nf` or `.config` file.

| Asset | What it drives | Extend by |
|---|---|---|
| `al_fusion_dictionary.tsv` | 200 named recurrent pairs | adding a row |
| `al_fusion_anchors.tsv` | 54 promiscuous anchors reported partner-agnostically | adding a row |
| `focal_duplication_targets.tsv` | which genes are scanned for ITD/PTD, on which reference, with which thresholds | adding a row |
| `gene_aliases.tsv` | synonyms VEP emits that no BED label carries | adding a row |
| `excluded_loci.tsv` | loci excluded from the SNV report (Ig, TR) | adding a row |
| `<PANEL>_panel_*.bed` | panel membership, calling regions, gene intervals | rebuilding with `bin/build_panel.py` |

The focal-duplication callers take a gene symbol and resolve its interval from
the panel BED at runtime, so they cannot drift out of step with the panel.
`build_panel.py` builds the hg38 BEDs from a target table plus RefSeq, so the
panel itself is reproducible from its design table.

## Two priors acute leukaemia needs that myeloma did not

**Promiscuous anchors.** KMT2A has over a hundred described partners, NUP98
more than thirty, and any RARA junction is APL until proven otherwise. A
named-pair dictionary cannot cover these, and missing a KMT2A junction because
its partner was novel would be a clinically material false negative. Anchors
are reported whenever they appear on either side of a junction, whatever the
partner is, including when the partner is off-panel.

**Single-sided cytoband matching.** Two defining B-ALL entities have an
off-panel partner by design: PBX1 in TCF3::PBX1, and DUX4 in IGH::DUX4, which
is unmappable anyway inside the D4Z4 macrosatellite. When one side is on-panel
and the other resolves only to a cytoband, a dictionary row carrying the
expected band for its partner still names the event. Band comparison is
tiered — exact sub-band, then major band, then arm — because T2T-CHM13 band
boundaries differ from GRCh38's and exact-string matching would lose a
defining entity to a one-band shift. Arm-level agreement is recorded as
`partial_arm` and should be read as a lead, not a call.

## Quick start

```bash
git clone https://github.com/patkarlab/leukemia_awgs
cd leukemia_awgs

# Confirm all four descriptions of each panel agree before anything else.
# Exits non-zero on a missing or empty BED; --strict also fails on warnings.
python3 bin/check_panel_consistency.py --panel AML --dictionary
python3 bin/check_panel_consistency.py --panel ALL --dictionary

# Smoke test, no real data needed
nextflow run main.nf -profile stub -stub-run --sample_sheet tests/sample_sheet_stub.csv

# Real run
cp assets/sample_sheet_template.csv my_samples.csv
$EDITOR my_samples.csv
nextflow run main.nf -profile conda,docker,gandalf \
    --sample_sheet my_samples.csv --outdir results/
```

Reference paths default to the gandalf locations; override any with
`--<name> <value>`. See `--help`.

## v0.1 scope

**In:** dual-reference realignment; three-caller SV ensemble with SURVIVOR
merge; AL fusion annotation with promiscuous anchors and cytoband fallback;
ClairS-TO, ichorCNA, Clair3 phased germline and VEP; AL panel variant filter;
tandem-duplication and intragenic-duplication callers; on-target QC; per-sample
panel resolution.

**Deferred to v0.2:** report bundle and HTML dashboard (the myeloma versions
have MM-specific parsers); IGV snapshot generation; BAF/LOH screen; ELN 2022
risk assignment; cohort summary assembly.

**Deferred to v0.3:** the SPARSH methylation classifier hook; iAMP21 detection
from on-target depth; ploidy classification.

## Known limitations

1. **The ALL T2T BEDs carry liftover drift the hg38 BEDs do not.** The panel
   is designed in hg38 (`bin/build_panel.py` reads RefSeq hg38), and the T2T
   BEDs are downstream of it. Three ALL regions are materially larger on T2T,
   totalling 2.33 Mb, which is 5.7% of the T2T panel:

   | Region | hg38 | T2T | Ratio |
   |---|---|---|---|
   | `RANBP2` | 166,328 bp | 1,223,906 bp | 7.4x |
   | `IGK_locus` | 1,450,000 bp | 2,324,577 bp | 1.6x |
   | `IGL_locus` | 910,000 bp | 1,306,194 bp | 1.4x |

   The Ig loci are defensible: they are genuinely larger and fully assembled
   in T2T, where the hg38 `NAMED_REGIONS` windows are fixed intervals that may
   under-cover. `RANBP2` is not. It sits at 2q13 inside a segmental
   duplication block, and a 7.4x expansion of a gene-body-plus-flank window is
   the signature of a chain mapping across an SD rather than a real coverage
   requirement.

   This costs sequencing, not just tidiness: the NC_-named T2T BED is the
   MinKNOW adaptive-sampling configuration, so 1.06 Mb of flowcell capacity
   (2.6% of the ALL panel) is being spent on that expansion. `HOXB_cluster` on
   the AML panel shows the same pattern at smaller scale (110 kb to 310 kb,
   0.7% of that panel). Worth re-deriving the T2T BEDs from the current hg38
   design and checking `RANBP2` specifically.

2. **The ALL T2T chrY coverage does not match the design.** The design
   specifies `PAR1_CRLF2_P2RY8_Y` as a fixed 500 kb window at chrY:1,150,000
   to 1,650,000, and the hg38 BED has it. The T2T BEDs instead carry a 175 kb
   `P2RY8` gene-body region at chrY:1,275,950 to 1,451,516. P2RY8::CRLF2
   arises from an interstitial PAR1 deletion whose breakpoints are distributed
   across PAR1, so a gene-body window will miss breakpoints falling outside
   P2RY8 itself.

   Two further consequences follow from PAR1 being real sequence on both X and
   Y in T2T-CHM13v2.0 rather than N-masked as in the hg38 analysis set. PAR1
   reads become multi-mapping, so anything filtering on MAPQ sees depleted
   coverage there, which includes on-target QC and the intragenic-duplication
   caller at its default `min_mapq = 20`. Sniffles and CuteSV run at
   `mapq = 0` and are unaffected. Worth confirming how your T2T reference
   handles the Y PARs before reading PAR1 coverage numbers.

3. **QC and caller thresholds are inherited, not re-derived.**
   `qc_depth_threshold = 15` and the SV callers' `min_support = 2` come from
   the myeloma pipeline, where the panel was 0.776% of the genome. These
   panels are 0.96% and 1.31%, so adaptive-sampling enrichment per region is
   lower and per-region depth should be roughly proportionally lower at equal
   flowcell output. Lower depth at `min_support = 2` pushes toward false
   negatives, so it is sensitivity that needs rechecking against a real AL run,
   not specificity.

4. **The duplication callers are untested against real data.** The
   tandem-duplication caller recovers a synthetic 45 bp duplication at the
   correct position, length, frame and allelic ratio, and correctly rejects a
   gene absent from the panel BED. The intragenic-duplication caller has no
   test beyond stub execution. Both report evidence, not diagnoses; confirm
   with an orthogonal assay before clinical reporting.

5. **Allelic ratio is a read-level ratio.** The clinical ELN mutant:wild-type
   ratio from fragment-analysis PCR is a different measurement. They are not
   interchangeable without local validation.

6. **Four `emerging`-tier dictionary rows carry `panel:` references** rather
   than citations, meaning they reflect the panel design's own designation and
   have not been independently verified. Curate them before clinical use.

## Panel and dictionary documentation

- [`assets/PANEL_README.md`](assets/PANEL_README.md) — panel composition,
  enrichment expectations, dictionary coverage, and gaps worth confirming
- [`assets/al_fusion_dictionary.README`](assets/al_fusion_dictionary.README) —
  dictionary and anchor schema, matching semantics, curation status

## Environment

| Step | Mechanism |
|---|---|
| SV callers, realignment, annotation, QC, duplication callers | conda `awgs_sv` (needs pysam) |
| ClairS-TO | docker `hkubal/clairs-to:latest` |
| ichorCNA | conda env `ichorCNA`, absolute paths, no activation |
| Clair3 | docker `hkubal/clair3:latest` |
| VEP | docker `ensemblorg/ensembl-vep:release_113.0` |

Docker-based steps invoke `docker run` inline rather than via Nextflow's
`container` directive, which would nest docker in docker. See
`conf/envs.config`.

## Porting notes

Four patterns in the myeloma repository do not compile on Nextflow 24.10 and
later and were rewritten here in a form that parses on 23.10 through 26.04: a
`def check_max()` function inside `conf/base.config`; top-level statements in
`main.nf`; an unescaped `${CONDA_PREFIX}` in a process script block; and a
bare-string `publishDir` referencing `meta`. A fifth is version-independent:
inside a workflow body the name `workflow` is shadowed, so completion handlers
must capture the metadata reference at registration time.

## Roadmap

- [x] v0.1 — dual-reference tracks, AL fusion annotation, duplication callers, panel switch
- [ ] v0.2 — report bundle, dashboard, IGV snapshots, BAF/LOH, ELN 2022 risk
- [ ] v0.3 — SPARSH methylation hook, iAMP21, ploidy classification
- [ ] v1.0 — validation cohort, peer-reviewed publication

## License

MIT. See [LICENSE](LICENSE).
