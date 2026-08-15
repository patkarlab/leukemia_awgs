# dashboard_builder v0.4.1-cancervar

Release date: 2026-05-22
Previous: v0.4.0-phase1+reporting

## What's new

### CancerVar annotation (AMP/ASCO/CAP 2017 somatic tiers)

A new opt-in build-time annotation step that adds, per clinical variant:

- A **four-level tier classification** (Tier I / II / III / IV) under the
  AMP/ASCO/CAP 2017 guideline -- the somatic equivalent of ACMG/AMP and the
  framework most molecular tumor boards report against.
- The **12 CBP criteria scores** (each 0/1/2) showing the rule-based reasoning
  behind the tier call (therapeutic association, hotspot membership,
  computational predictions, etc.).
- The **OPAI** ("Oncogenic Pathogenicity AI") score in [0, 1] from CancerVar's
  semi-supervised deep-learning model -- useful as a tiebreaker when the
  rule-based engine returns Tier III (e.g. the test sample's TET2 variant
  comes back Tier III but with OPAI 0.90, flagging it as likely oncogenic
  despite insufficient guideline-level evidence).

The endpoint is anonymous -- no token or API key required, no registration.

### Reporting tab tier pre-fill

When CancerVar annotation is enabled at build time, ticking a variant in the
Clinical tab pre-fills the editable tier column in the Reporting tab with the
CancerVar tier ("Tier I", "Tier II", etc.). The input remains fully editable
and persists to localStorage -- the pathologist can confirm, override, or
clear the value like any manual entry. Pre-fill only fires on the first
selection-on transition for a given variant; if the input has any existing
value it is left untouched.

### Reporting tab column rename

The editable tier column header changed from "ACMG/AMP Tier" to
"AMP/ASCO/CAP Tier (somatic)". The TSV export header matches. ACMG/AMP is
strictly the germline framework; AMP/ASCO/CAP (Li et al, JMD 2017) is the
correct citation for somatic variant tiering on a leukemia panel.

The input's `maxlength` was widened from 6 to 16 to accommodate "Tier III" /
"Tier IV" prefills and any short pathologist annotation.

## Files added

- `parsers/cancervar.py` -- new parser, mirrors the OncoKB module shape.
- `<sample>/<sample>_cancervar_cache.json` -- written by the builder when
  `--annotate-cancervar` is set, carries both positive entries (with tier,
  OPAI, 12 CBP rows) and negative entries (variants the API couldn't
  classify, cached so they aren't re-queried on rebuild).

## Files changed

- `build.py` -- added `--annotate-cancervar` flag (no token arg), wired
  through `build()` and `collect_sample_context()`,
  `_cancervar_cache.json` added to `FILE_DESCRIPTIONS`. Version bumped to
  `0.4.1-cancervar`.
- `assets/js/variant-browser.js`:
  - new `renderCancervarBlock()` between GeneBe and OncoKB, with tier pill
    (reuses ACMG color palette: I red / II orange / III gray / IV green),
    color-coded OPAI badge, and 12 CBP criteria in a collapsed `<details>`
    element.
  - extended the selection snapshot to carry `cancervar_tier`, with
    pre-fill of the editable tier on first selection-on if no manual
    value already exists.
- `templates/sample_report.html.j2`:
  - new `cancervarAnnotations` embedded global.
  - passed through to the clinical browser (filtered browser gets `{}` for
    the same payload-size reasons as the existing OncoKB pass-through).
  - column header rename, prose update, wider tier input.

## API endpoint

```
GET https://cancervar.wglab.org/api_new.php
    ?queryType=position&build=hg38&chr=<C>&pos=<P>&ref=<R>&alt=<A>
```

No authentication required. Politeness sleep of 0.5 s between requests --
the server is on older Apache/PHP 7.2 infrastructure, so don't hammer it.

## Behavior on the test sample (26CGH825-MYCNV)

12 clinical variants queried, 4 classified, 8 not classified, 0 errors:

| Variant            | Gene   | Tier    | OPAI |
|--------------------|--------|---------|------|
| chr2:25234373:C:T  | DNMT3A | Tier II | 0.54 |
| chr4:105272774:C:T | TET2   | Tier III| 0.90 |
| chr11:32435102:C:T | WT1    | Tier III| 0.01 |
| chr11:118505854:C:T| KMT2A  | Tier IV | 0.00 |

The 8 misses are all complex indels (FLT3-ITD-style insertions, NPM1-style
insertions, internal tandem duplications). CancerVar's curation is heavily
SNV-focused, so this is expected -- the negative cache prevents re-querying
on subsequent builds.

## Build invocation

```bash
python tools/dashboard_builder/build.py /path/to/run_dir \
    --verbose \
    --annotate-genebe --genebe-user EMAIL --genebe-key KEY \
    --annotate-oncokb --oncokb-token TOKEN \
    --annotate-cancervar
```

All four flags compose; they each gate an independent network pass with its
own per-sample cache.

## Known limitations and design notes

- CancerVar's API has no documented rate limit; the 0.5 s politeness sleep is
  conservative. On a 12-variant sample the annotation pass takes ~6 s.
- Negative cache entries are durable on disk. If CancerVar's database is
  updated upstream, a previously-unclassified variant will not be re-queried
  unless the cache file is deleted. This is intentional (matches OncoKB
  behavior); if upstream changes are suspected, delete the cache and rebuild.
- The tier pre-fill only fires on the first selection-on transition for a
  variant. If a user clears the input, then unchecks and re-checks, the
  CancerVar tier is offered again (because the cleared input means there is
  no existing value to preserve). This is a deliberate trade-off in favor of
  visibility over strict idempotence.

## Pending from v0.4.0 handoff (unchanged)

- Exon column in the TSV / from VEP.
- CNV selection for Reporting tab.
- Real-browser end-to-end test of OncoKB with a live token.
- Filtered variants externalization (768 KB inline JSON moves to a sidecar).
- MobiDetails audit cache visibility.
- End-to-end test under `python3 -m http.server`.

## Reference

Li Q, Ren Z, Cao K, Li MM, Zhou Y, Wang K. CancerVar: an Artificial
Intelligence empowered platform for clinical interpretation of somatic
mutations in cancer. Sci Adv. 2022.
https://www.science.org/doi/10.1126/sciadv.abj1624
