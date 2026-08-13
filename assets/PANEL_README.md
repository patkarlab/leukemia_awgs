# Panel design notes

Two adaptive-sampling panels, selected per sample via `meta.panel`.

| | AML | ALL |
|---|---|---|
| Design targets | 120 | 133 |
| T2T regions / bases | 117 / 29.94 Mb | 126 / 40.74 Mb |
| hg38 regions / bases | 117 / 29.43 Mb | 126 / 38.00 Mb |
| Fraction of T2T-CHM13v2.0 | 0.960% | 1.307% |
| Shared region labels | 54 | 54 |

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

**Present on one panel but not the other.** LYN, NTRK3, ABL2 and EWSR1 are on
the AML panel only, yet each participates in a described ALL fusion:
LYN::NCOR1 and ETV6::NTRK3 are Ph-like ABL-class and other-kinase fusions
respectively, ETV6::ABL2 is a kinase fusion of both lineages, and
ZNF384::EWSR1 is a recognised ZNF384-rearranged B-ALL. Each is TKI-relevant.
Worth confirming these omissions are intentional.

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
