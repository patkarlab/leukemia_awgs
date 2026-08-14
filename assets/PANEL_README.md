# Panel design notes

Two adaptive-sampling panels, selected per sample via `meta.panel`.

| | AML | ALL |
|---|---|---|
| Design targets | 120 | 136 |
| T2T regions / bases | 117 / 29.94 Mb | 128 / 40.66 Mb |
| hg38 regions / bases | 117 / 29.43 Mb | 128 / 38.70 Mb |
| Fraction of T2T-CHM13v2.0 | 0.960% | 1.304% |

Both panels report full consistency across all four of their descriptions
(`bin/check_panel_consistency.py --panel <PANEL> --dictionary`).

The ALL panel was revised after its T2T BEDs were first built: sixteen
MPN/eosinophilia partner genes were dropped, keeping the five kinases
themselves (ABL1, BCR, FGFR1, PDGFRA, PDGFRB), and `PAR1_CRLF2_P2RY8_Y` was
added for the chrY pseudoautosomal copy. The T2T BEDs here have been trimmed
to match; see the known limitations in the top-level README for what still
needs re-deriving from the current hg38 design.

Files, per panel: `<PANEL>_panel_t2t_chr.bed` (chr-named, used by the
pipeline), `<PANEL>_panel_t2t_NC.bed` (NC_-named, for the MinKNOW adaptive
sampling configuration), `<PANEL>_panel_hg38.bed` (SNV calling, currently a
placeholder), and `<PANEL>_panel_targets.tsv` (the design table, with a
`GROUP` column giving each region's rationale).

## Enrichment expectation differs from myeloma

The MM v7 PCN panel covers 0.776% of the genome. These panels cover 0.96% and
1.31%. Adaptive sampling rejects off-target reads, so per-region depth falls
roughly in proportion as the target fraction rises: at equal flowcell output,
the ALL panel should yield roughly 60% of the per-region depth the MM panel
does.

Two consequences. First, `qc_depth_threshold` is inherited from MM at 15 and
has not been re-derived for these panels; treat it as a placeholder until you
have on-target depth from a real AL run. Second, the SV callers are configured
at `min_support = 2` with `mapq = 0`, tuning that was validated at MM depths.
Lower depth pushes toward false negatives at that setting rather than false
positives, so it is the sensitivity that needs rechecking, not the specificity.

## Dictionary coverage

Of 200 rows in `al_fusion_dictionary.tsv`, 186 have both partners present in
the corresponding panel BED. The 14 that do not fall into two groups.

**Off-panel by design.** PBX1 (TCF3::PBX1), DUX4 (IGH::DUX4), CUX1 and BRD9
(NUTM1 partners), ARID1B (BCL11B partner), LAIR1 and SMU1 (Ph-like JAK-class
partners), ZFAND3, STAG2. These carry a `partner_b_band` so the annotator can
still name the event from the on-panel side plus the partner's cytoband. For
TCF3::PBX1 and IGH::DUX4 this is the only route to a call, so the mechanism is
load-bearing, not a convenience.

DUX4 additionally sits in the D4Z4 macrosatellite and is not reliably
mappable even when targeted. The intended detection route is an IGH breakpoint
at 4q35 together with the ERG intragenic deletion that accompanies the
subtype; ERG is on the ALL panel.

**Present on one panel but not the other.** LYN and NTRK3 were AML-only
despite each participating in a described Ph-like ALL fusion whose partner was
already on the ALL panel: LYN::NCOR1 is ABL-class and dasatinib-sensitive,
ETV6::NTRK3 is other-kinase and TRK-inhibitor sensitive under a tumour-agnostic
approval. Both have been added to the ALL panel, reusing the AML intervals
verbatim since the build parameters and reference are identical. Together they
add 731,461 bp, 1.80% of the ALL panel — less than the RANBP2 over-expansion
described in the top-level README, so fixing that funds both with room over.

ABL2 and EWSR1 remain single-panel. ETV6::ABL2 is a kinase fusion of both
lineages and ZNF384::EWSR1 is a recognised ZNF384-rearranged B-ALL, but in each
case the on-panel side carries the subtype signal and the off-panel partner has
a cytoband fallback in the dictionary, so the entity is still nameable. Worth
confirming these two omissions are intentional.

## Gaps worth confirming

These are absences relative to what a comprehensive panel for each disease
would carry. If SNVs are covered by a separate assay and these panels are
deliberately fusion-focused, several of them are not gaps at all — but they
are worth an explicit decision rather than an implicit one.

**AML.** DNMT3A and TET2 are the two most frequently mutated genes in adult
AML and neither is present. NRAS, KRAS and PTPN11 cover the RAS pathway and
are absent. WT1 is absent despite modifying risk in NPM1-mutant and
NUP98::NSD1 disease. DDX41 is absent despite now being a recognised germline
predisposition entity in both WHO-5 and ICC. RAD21, SMC1A, SMC3 and PHF6
complete the cohesin and chromatin comutation set that ASXL1, BCOR, EZH2 and
STAG2 have started.

**ALL.** CDKN2A/CDKN2B are absent. The 9p21 deletion is among the commonest
lesions in both B-ALL and T-ALL and is a standard part of the diagnostic
workup. NOTCH1 and FBXW7 are absent, and they are the two most frequently
mutated genes in T-ALL, with NOTCH1/FBXW7 status entering several risk
algorithms. PTEN is absent despite marking glucocorticoid resistance in T-ALL.
PAX5 is present, which covers PAX5 P80R, and IKZF1 is present, which covers
both the deletion and N159Y.

**Both.** Adaptive sampling captures the whole targeted interval including
introns and regulatory sequence, so a gene on the panel is fully available for
SNV, indel and breakpoint detection. That makes the absent genes above a real
coverage limit rather than a resolution limit.

## Detection routes that need something other than a junction call

Several defining lesions produce no gene-to-gene junction and will not appear
in the fusion table however well the panel performs.

`FLT3-ITD` and `KMT2A-PTD` are handled by dedicated callers in this pipeline
(`bin/call_flt3_itd.py`, `bin/call_kmt2a_ptd.py`). `UBTF-TD` is an exon-13
tandem duplication and `MYB` duplication is intragenic; both are marked
`partner_policy=self` in the anchor table but have no dedicated caller yet,
and would need one modelled on the KMT2A-PTD approach.

High hyperdiploidy, near-haploidy and low hypodiploidy in B-ALL are ploidy
states. The ichorCNA ploidy grid has been widened from the myeloma `c(2,3)` to
`c(1,2,3,4)` and `maxCN` from 5 to 6 to accommodate them, but ichorCNA was
built for tumour-fraction estimation from low-pass off-target reads and is not
a validated ploidy classifier. Treat its output as a screen.

iAMP21 needs finer resolution than the 1 Mb ichorCNA bins can give. The
amplified region is typically a few megabases of 21q22 and the diagnostic
criterion is copy number of the RUNX1 region specifically. RUNX1 is on both
panels, so per-region on-target depth from `QC_ONTARGET` is a better starting
point than the genome-wide CNV track. This is not yet automated.

## Change log

### CDKN2A and CDKN2B added to the ALL panel

9p21 deletion occurs in roughly 70% of T-ALL and 30% of B-ALL, making it one
of the commonest lesions in the disease, and it is the fourth term of the
IKZF1plus classifier alongside IKZF1, PAX5 and PAR1, all three of which were
already on the panel.

The two genes lie about 8 kb apart and merge at the default flank into one
region of roughly 141.6 kb: `chr9:21,917,751-22,059,313` in hg38 and
`chr9:21,932,051-22,073,690` in T2T. The ~14 kb offset between them is the
genuine assembly difference at 9p21. The default flank is deliberate; it
supports deletion detection by depth, which is what the clinical call needs,
rather than breakpoint mapping, which would not resolve at this panel's
on-target depth.

### RANBP2 removed from the ALL panel

T2T RefSeq annotates RANBP2 as a single 1,123,905 bp gene model spanning the
2q13 RGPD paralogue cluster, against roughly 66 kb in hg38. This was verified
as a real annotation rather than a liftover artefact: the GFF holds one gene
entry at that span, so re-deriving returns the same result.

Keeping it cost 1.22 Mb of adaptive-sampling target, 2.9% of the ALL panel,
across segmental duplication where reads are multi-mapping and discarded by
anything filtering on MAPQ. The capacity was spent on reads most of the
pipeline throws away.

`RANBP2::ABL1` detection is unaffected. ABL1 is on the panel and is a
promiscuous anchor, so any ABL1 junction is reportable whatever the partner,
and the dictionary row carries `partner_b_band=2q13` so the entity is still
named from the ABL1 side. What is lost is breakpoint resolution on the RANBP2
side.

This is a general hazard rather than a RANBP2 one. Any gene sitting in a
T2T-assembled segmental duplication may carry a similarly inflated model. The
`--hg38-bed` span comparison in `build_panel_t2t.py` is what surfaces them;
`ZNF362` at 1.91x is the remaining candidate on the ALL panel.

### chrY PAR1 window corrected on the ALL T2T panel

The design specifies `PAR1_CRLF2_P2RY8_Y` as a 500 kb interval and the hg38
BED carried it, but the T2T BEDs held a 175 kb P2RY8 gene-body region instead.
Because adaptive sampling runs against T2T, that 325 kb was never enriched:
P2RY8::CRLF2 arises from an interstitial PAR1 deletion whose breakpoints are
distributed across PAR1, so breakpoints outside the P2RY8 gene body were
absent from the data rather than merely unreported. Samples sequenced against
the previous NC_ BED cannot be recovered by reanalysis.

Replaced with `chrY:1,000,430-1,451,516`, 451,086 bp, spanning CRLF2 (start
1,050,430) to P2RY8 (end 1,401,516) in T2T coordinates plus 50 kb margin,
supplied through `assets/extra_regions_ALL_t2t.bed`.

## Known limitations

**`build_panel_t2t.py` matches labels by token, not exactly.** The X-side
`PAR1_CRLF2_P2RY8` and the Y-side `PAR1_CRLF2_P2RY8_Y` share every token but
`Y`, so they collide and the emitted chrY label needed manual correction. A
rebuild of the ALL panel will reproduce that and needs the same repair until
named intervals are matched exactly. Token matching is correct for compound
labels like `TAL1/STIL` and wrong for X/Y paralogous intervals.

**Two spans still exceed their hg38 counterparts.** `ZNF362` at 1.91x is
unexplained and worth checking against the T2T annotation the way RANBP2 was.
`IGK_locus` at 1.60x and the X-side PAR1 interval at 1.80x are expected: the
Ig loci are genuinely larger and fully assembled in T2T, where the hg38
`NAMED_REGIONS` entries are fixed windows that may under-cover.

**AML 17p / TP53 window deferred.** TP53 currently receives a centred window
of plus or minus 500 kb. Widening it was considered and deferred pending
review of data from the batch currently sequencing. Note when revisiting that
`CENTERED` in `build_panel.py` is shared across both panels and TP53 is
`SCOPE=BOTH`, so changing its half-width widens TP53 on the ALL panel too
unless `CENTERED` is made panel-aware first.

## Regenerating the derived records

`BASES` in the target tables and `panel_summary.tsv` are computed, not
authored. After any panel change:

    python3 bin/update_panel_records.py \
        --assets  assets \
        --refgene <path>/refGene.txt.gz \
        --sizes   <path>/hg38.chrom.sizes

`BASES` is the per-target span before merging, matching what `build_panel.py`
accumulates. It is not the sum of the emitted BED regions: overlapping targets
merge into one region whose span is smaller than the individual spans added
together.
