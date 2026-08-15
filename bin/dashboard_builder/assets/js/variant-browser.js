/* variant-browser.js — master-detail UI for somaticseq variant TSVs.
 *
 * Used by both the "Variants - Clinical" and "Variants - All Filtered" tabs.
 * Initialised once per tab via initVariantBrowser(config).
 *
 * config = {
 *   containerId:     string         // id of the empty <div> to render into
 *   variants:        Array<Object>  // rows from the TSV
 *   filterIds:       Array<string>  // which filters to show (keys of FILTER_DEFS)
 *   igvLookup:       Object|null    // {chr:pos:ref:alt -> unique_id} or null
 *   igvFrameId:      string|null    // id of the IGV iframe
 *   igvTabSelector:  string|null    // selector of the IGV tab button (for tab switching)
 *   genebeAnnotations: Object|null  // {chr:pos:ref:alt -> annotation dict} or null
 * }
 */

(function () {
  "use strict";

  // ---------------------------------------------------------------------------
  // Reporting selection store (window-level so the Reporting tab can read it)
  //
  // We track a per-sample list of selected items in localStorage. Each item is
  //   { kind: "variant"|"cnv"|..., id: string, snapshot: {...} }
  // The snapshot carries everything the Reporting table needs to render the row
  // without re-walking the variant TSV. This survives page reloads and is
  // independent of the variant TSV being re-embedded across builds.
  //
  // A separate per-sample, per-variant key holds the editable ACMG/AMP Tier
  // value the pathologist types into the Reporting table (also localStorage,
  // so it persists across reloads).
  // ---------------------------------------------------------------------------
  window.tspipeReporting = window.tspipeReporting || (function () {
    const SELECTION_KEY_PREFIX = "tspipe-report:";
    const TIER_KEY_PREFIX = "tspipe-tier:";
    let sampleKey = null;
    const listeners = [];

    function storageAvail() {
      try {
        const t = "__tspipe_probe__";
        window.localStorage.setItem(t, t);
        window.localStorage.removeItem(t);
        return true;
      } catch (e) {
        return false;
      }
    }
    const HAS_STORAGE = storageAvail();

    function setSample(name) { sampleKey = name; }

    function load() {
      if (!sampleKey || !HAS_STORAGE) return [];
      try {
        const raw = window.localStorage.getItem(SELECTION_KEY_PREFIX + sampleKey);
        return raw ? JSON.parse(raw) : [];
      } catch (e) { return []; }
    }

    function save(items) {
      if (!sampleKey || !HAS_STORAGE) return;
      try {
        window.localStorage.setItem(SELECTION_KEY_PREFIX + sampleKey, JSON.stringify(items));
      } catch (e) { /* quota or disabled */ }
      listeners.forEach(function (cb) {
        try { cb(items); } catch (e) { /* one listener's bug shouldn't break others */ }
      });
    }

    function isSelected(kind, id) {
      // "selected" historically meant "in the report". With three-state triage,
      // an item can be present but marked "exclude"; those are NOT selected.
      // Legacy items written before triage have no decision -> treated include.
      return load().some(function (it) {
        return it.kind === kind && it.id === id && (it.decision || "include") === "include";
      });
    }

    // Three-state triage: "" (neutral, not in store) | "include" | "exclude".
    function getDecision(kind, id) {
      const it = load().find(function (x) { return x.kind === kind && x.id === id; });
      if (!it) return "";
      return it.decision || "include";
    }

    function setDecision(kind, id, decision, snapshot) {
      const items = load();
      const idx = items.findIndex(function (it) { return it.kind === kind && it.id === id; });
      if (!decision) {
        // Neutral: drop the item entirely and clear its exclude reason.
        if (idx >= 0) items.splice(idx, 1);
        setExcludeReason(id, "");
        save(items);
        return;
      }
      if (idx >= 0) {
        items[idx].decision = decision;
        if (snapshot) items[idx].snapshot = snapshot;
      } else {
        items.push({ kind: kind, id: id, decision: decision, snapshot: snapshot || {} });
      }
      if (decision !== "exclude") setExcludeReason(id, "");
      save(items);
    }

    // Exclude reason -- the short note on why a reviewed variant is not reported.
    // Same persistence model as tier values.
    const EXCLUDE_REASON_KEY_PREFIX = "tspipe-exclude-reason:";

    function getExcludeReason(id) {
      if (!sampleKey || !HAS_STORAGE) return "";
      try {
        return window.localStorage.getItem(EXCLUDE_REASON_KEY_PREFIX + sampleKey + ":" + id) || "";
      } catch (e) { return ""; }
    }

    function setExcludeReason(id, text) {
      if (!sampleKey || !HAS_STORAGE) return;
      try {
        if (text) {
          window.localStorage.setItem(EXCLUDE_REASON_KEY_PREFIX + sampleKey + ":" + id, text);
        } else {
          window.localStorage.removeItem(EXCLUDE_REASON_KEY_PREFIX + sampleKey + ":" + id);
        }
      } catch (e) {}
    }

    function toggle(kind, id, snapshot) {
      const items = load();
      const idx = items.findIndex(function (it) { return it.kind === kind && it.id === id; });
      if (idx >= 0) {
        items.splice(idx, 1);
      } else {
        items.push({ kind: kind, id: id, snapshot: snapshot || {} });
      }
      save(items);
    }

    function clearAll() { save([]); }

    function getTier(id) {
      if (!sampleKey || !HAS_STORAGE) return "";
      try {
        return window.localStorage.getItem(TIER_KEY_PREFIX + sampleKey + ":" + id) || "";
      } catch (e) { return ""; }
    }

    function setTier(id, tier) {
      if (!sampleKey || !HAS_STORAGE) return;
      try {
        if (tier) {
          window.localStorage.setItem(TIER_KEY_PREFIX + sampleKey + ":" + id, tier);
        } else {
          window.localStorage.removeItem(TIER_KEY_PREFIX + sampleKey + ":" + id);
        }
      } catch (e) {}
    }

    // CNV interpretation captions -- the editable text the pathologist writes
    // below each selected CNV plot. Same persistence model as tier values.
    const CNV_CAPTION_KEY_PREFIX = "tspipe-cnv-caption:";

    function getCnvCaption(id) {
      if (!sampleKey || !HAS_STORAGE) return "";
      try {
        return window.localStorage.getItem(CNV_CAPTION_KEY_PREFIX + sampleKey + ":" + id) || "";
      } catch (e) { return ""; }
    }

    function setCnvCaption(id, text) {
      if (!sampleKey || !HAS_STORAGE) return;
      try {
        if (text) {
          window.localStorage.setItem(CNV_CAPTION_KEY_PREFIX + sampleKey + ":" + id, text);
        } else {
          window.localStorage.removeItem(CNV_CAPTION_KEY_PREFIX + sampleKey + ":" + id);
        }
      } catch (e) {}
    }

    function onChange(cb) { listeners.push(cb); }

    // Extract just the COSMIC identifiers (COSV*, COSM*, COSN*) from VEP's
    // ampersand-delimited Existing_variation column. Returns ID=...;ID=...
    // (semicolon-joined) or "" if none are present.
    function extractCosmicIds(existingVariation) {
      if (!existingVariation || existingVariation === "-1") return "";
      return String(existingVariation).split("&")
        .filter(function (s) { return /^COS[VMN]/.test(s); })
        .map(function (s) { return "ID=" + s; })
        .join(";");
    }

    return {
      setSample: setSample,
      load: load,
      isSelected: isSelected,
      getDecision: getDecision,
      setDecision: setDecision,
      getExcludeReason: getExcludeReason,
      setExcludeReason: setExcludeReason,
      toggle: toggle,
      clearAll: clearAll,
      getTier: getTier,
      setTier: setTier,
      getCnvCaption: getCnvCaption,
      setCnvCaption: setCnvCaption,
      onChange: onChange,
      extractCosmicIds: extractCosmicIds,
      hasStorage: function () { return HAS_STORAGE; },
    };
  })();

  // ---------------------------------------------------------------------------
  // Filter definitions
  // ---------------------------------------------------------------------------
  const FILTER_DEFS = {
    somaticseq_verdict: {
      label: "SomaticSeq Verdict",
      field: "SomaticSeq_Verdict",
      type: "multi-checkbox",
    },
    filter_status: {
      label: "Filter Status",
      field: "Filter",
      type: "multi-checkbox",
    },
    impact: {
      label: "IMPACT",
      field: "IMPACT",
      type: "multi-checkbox",
      preferredOrder: ["HIGH", "MODERATE", "LOW", "MODIFIER"],
    },
    consequence: {
      label: "Consequence",
      field: "Consequence",
      type: "multi-checkbox",
      splitOn: "&",
      scroll: true,
    },
    oncovi: {
      label: "OncoVI Classification",
      field: "OncoVI_Classification",
      type: "multi-checkbox",
    },
    gene: {
      label: "Gene contains",
      field: "Gene",
      type: "text",
    },
    vaf_min: {
      label: "VAF % \u2265",
      field: "VAF_pct",
      type: "min-range",
      min: 0, max: 100, step: 0.5, default: 0,
    },
    // Button-group filters: clinical convention uses discrete thresholds rather than
    // a continuous slider. Each option is {label, predicate(rowValueAsNumber)}.
    alt_count_buttons: {
      label: "ALT count",
      field: "ALT_COUNT",
      type: "button-group",
      options: [
        { id: "any", label: "Any",  test: function () { return true; } },
        { id: "gt1", label: ">1",  test: function (v) { return v !== null && v > 1; } },
        { id: "gt2", label: ">2",  test: function (v) { return v !== null && v > 2; } },
        { id: "gt5", label: ">5",  test: function (v) { return v !== null && v > 5; } },
      ],
      default: "any",
    },
    pop_af_max: {
      label: "Max population AF \u2264",
      field: "Max_AF",
      type: "max-range-novel-aware",
      min: 0, max: 1, step: 0.0001, default: 1.0,
      novelValues: ["", "-1", "-1.0"],
    },
  };

  const SORT_OPTIONS = [
    { id: "original",    label: "Original order" },
    { id: "gene",        label: "Gene (A\u2192Z)" },
    { id: "vaf_desc",    label: "VAF (high\u2192low)" },
    { id: "vaf_asc",     label: "VAF (low\u2192high)" },
    { id: "alt_desc",    label: "ALT count (high\u2192low)" },
  ];

  // Detail-view field groupings. Fields not present in the row are skipped.
  const DETAIL_GROUPS = [
    ["Identity", [
      ["Sample", "Sample"],
      ["Gene", "Gene"],
      ["Chr", "Chr"],
      ["Start", "Position (start)"],
      ["End", "Position (end)"],
      ["Ref", "Ref"],
      ["Alt", "Alt"],
      ["COSMIC_ID", "COSMIC ID"],
      ["rsID", "rsID"],
      ["Existing_variation", "Existing variation"],
    ]],
    ["Annotation", [
      ["Consequence", "Consequence"],
      ["IMPACT", "IMPACT"],
      ["HGVSc", "HGVSc"],
      ["HGVSp", "HGVSp"],
      ["HGVSg", "HGVSg"],
      ["MANE_SELECT", "MANE Select"],
      ["Canonical", "Canonical"],
    ]],
    ["VariantValidator", [
      ["VV_HGVSc", "VV HGVSc"],
      ["VV_HGVSp", "VV HGVSp"],
      ["VV_HGVSg", "VV HGVSg"],
      ["VV_Transcript", "VV transcript"],
      ["VV_Valid", "VV valid"],
      ["VV_Warnings", "VV warnings"],
    ]],
    ["Calls", [
      ["VariantCaller_Count", "Caller count"],
      ["Callers", "Callers"],
      ["REF_COUNT", "REF count"],
      ["ALT_COUNT", "ALT count"],
      ["VAF_pct", "VAF %"],
      ["SomaticSeq_Verdict", "SomaticSeq verdict"],
      ["Filter", "Filter status"],
    ]],
    ["Pathogenicity", [
      ["SIFT", "SIFT"],
      ["PolyPhen", "PolyPhen"],
      ["ClinVar", "ClinVar"],
      ["OncoVI_Score", "OncoVI score"],
      ["OncoVI_Classification", "OncoVI classification"],
      ["OncoVI_Criteria", "OncoVI criteria"],
    ]],
    ["Population frequency", [
      ["gnomAD_exome_AF", "gnomAD exome AF"],
      ["gnomAD_genome_AF", "gnomAD genome AF"],
      ["AF_1KG", "1000 Genomes AF"],
      ["Max_AF", "Max AF"],
    ]],
    ["Pipeline notes", [
      ["Dedup_Note", "Dedup note"],
      ["Blacklist_Reason", "Blacklist reason"],
      ["Blacklist_Date", "Blacklist date"],
      ["HGVS_ITD", "HGVS-ITD"],
      ["Confirmed_by_FLT3_ITD_ensemble", "Confirmed by FLT3 ITD ensemble"],
    ]],
  ];

  // ---------------------------------------------------------------------------
  // Helpers
  // ---------------------------------------------------------------------------
  function asNum(v) {
    if (v === undefined || v === null || v === "" || v === "-1" || v === "-1.0") return null;
    const n = parseFloat(v);
    return isNaN(n) ? null : n;
  }

  function emptyVal(v) {
    if (v === undefined || v === null || v === "") return '<span class="text-muted">\u2014</span>';
    if (v === "-1" || v === "-1.0") return '<span class="text-muted">not reported</span>';
    return escapeHtml(String(v));
  }

  function escapeHtml(s) {
    return String(s).replace(/[&<>"']/g, function (c) {
      return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c];
    });
  }

  function bestHGVSp(row) {
    const vv = row.VV_HGVSp;
    if (vv && vv !== "-1") return vv;
    const hg = row.HGVSp;
    if (hg && hg !== "-1") return hg;
    return "";
  }

  function bestHGVSc(row) {
    const vv = row.VV_HGVSc;
    if (vv && vv !== "-1") return vv;
    const hg = row.HGVSc;
    if (hg && hg !== "-1") return hg;
    return "";
  }

  function verdictBadgeClass(v) {
    if (v === "PASS") return "bg-success";
    if (v === "REJECT") return "bg-secondary";
    if (v === "LowQual") return "bg-warning text-dark";
    return "bg-light text-dark border";
  }

  // ACMG classification CSS class lookup (GeneBe returns e.g. "Pathogenic", "Likely Pathogenic")
  function acmgPillClass(classification) {
    if (!classification) return "vb-acmg-uncertain";
    const norm = classification.toLowerCase().replace(/[_\s]+/g, "-");
    if (norm.indexOf("likely-pathogenic") !== -1) return "vb-acmg-likely-pathogenic";
    if (norm.indexOf("pathogenic") !== -1) return "vb-acmg-pathogenic";
    if (norm.indexOf("likely-benign") !== -1) return "vb-acmg-likely-benign";
    if (norm.indexOf("benign") !== -1) return "vb-acmg-benign";
    return "vb-acmg-uncertain";
  }

  function uniqueValues(rows, field, splitOn) {
    const counts = new Map();
    for (const r of rows) {
      let raw = r[field];
      if (raw === undefined || raw === null || raw === "") continue;
      const tokens = splitOn ? String(raw).split(splitOn) : [String(raw)];
      for (const t of tokens) {
        const k = t.trim();
        if (!k) continue;
        counts.set(k, (counts.get(k) || 0) + 1);
      }
    }
    return Array.from(counts.entries()).sort(function (a, b) { return b[1] - a[1]; });
  }

  // Strip chr prefix for GeneBe URL slug.
  function genebeUrl(chr, pos, ref, alt) {
    const c = String(chr || "").replace(/^chr/, "");
    return "https://genebe.net/variant/hg38/chr" + c + "-" + pos + "-" + ref + "-" + alt;
  }

  // True iff a TSV cell carries a real value (not blank and not the
  // "-1" sentinel that somaticseq writes for missing fields).
  function hasUsefulValue(v) {
    return v !== undefined && v !== null && v !== "" && v !== "-1" && v !== "-1.0";
  }

  // ---------------------------------------------------------------------------
  // Main entry point
  // ---------------------------------------------------------------------------
  window.initVariantBrowser = function (config) {
    const root = document.getElementById(config.containerId);
    if (!root) return;

    const variants = config.variants || [];
    const filterIds = config.filterIds || [];
    const igvLookup = config.igvLookup || {};
    const hasIGV = !!config.igvFrameId;
    const genebeAnnotations = config.genebeAnnotations || {};
    const oncokbAnnotations = config.oncokbAnnotations || {};
    const cancervarAnnotations = config.cancervarAnnotations || {};
    const enableReportSelect = !!config.enableReportSelect;

    variants.forEach(function (r, i) {
      r._idx = i;
      r._igvKey = (r.Chr || "") + ":" + (r.Start || "") + ":" + (r.Ref || "") + ":" + (r.Alt || "");
    });

    // Keep the triage buttons in sync when the selection store changes from
    // elsewhere (e.g. the Reporting tab's "Clear all"). We don't re-render the
    // whole list -- just flip the .active class on the affected buttons.
    if (enableReportSelect && window.tspipeReporting) {
      window.tspipeReporting.onChange(function (items) {
        const decisionByKey = {};
        items.filter(function (it) { return it.kind === "variant"; })
             .forEach(function (it) { decisionByKey[it.id] = it.decision || "include"; });
        const btns = listEl.querySelectorAll(".vb-triage-btn");
        btns.forEach(function (btn) {
          const k = btn.getAttribute("data-vb-key");
          const want = btn.getAttribute("data-decision");
          const shouldBe = decisionByKey[k] === want;
          btn.classList.toggle("active", shouldBe);
        });
      });
    }

    const state = {
      checkboxes: {},    // empty Set = "all checked", no constraint
      gene: "",
      sliders: {},
      buttonGroups: {},  // filterId -> active option id
      sortBy: "original",
      expandedIdx: null,
    };

    for (const id of filterIds) {
      const def = FILTER_DEFS[id];
      if (!def) continue;
      if (def.type === "multi-checkbox") state.checkboxes[id] = new Set();
      else if (def.type === "min-range" || def.type === "max-range-novel-aware") state.sliders[id] = def.default;
      else if (def.type === "button-group") state.buttonGroups[id] = def.default;
    }

    root.innerHTML =
      '<div class="row g-3">' +
        '<div class="col-lg-3 col-md-4">' +
          '<div class="card vb-filter-panel">' +
            '<div class="card-header d-flex justify-content-between align-items-center">' +
              '<span class="fw-semibold">Filters</span>' +
              '<button type="button" class="btn btn-sm btn-outline-secondary vb-reset">Reset</button>' +
            '</div>' +
            '<div class="card-body vb-filter-body"></div>' +
          '</div>' +
        '</div>' +
        '<div class="col-lg-9 col-md-8">' +
          '<div class="d-flex justify-content-between align-items-center mb-2 flex-wrap gap-2">' +
            '<div class="vb-count text-muted small"></div>' +
            '<div class="d-flex align-items-center gap-2">' +
              '<label class="text-muted small mb-0">Sort:</label>' +
              '<select class="form-select form-select-sm vb-sort" style="width:auto;"></select>' +
            '</div>' +
          '</div>' +
          '<div class="vb-list"></div>' +
        '</div>' +
      '</div>';

    const filterBody = root.querySelector(".vb-filter-body");
    const listEl     = root.querySelector(".vb-list");
    const countEl    = root.querySelector(".vb-count");
    const sortEl     = root.querySelector(".vb-sort");
    const resetBtn   = root.querySelector(".vb-reset");

    sortEl.innerHTML = SORT_OPTIONS.map(function (o) {
      return '<option value="' + o.id + '">' + escapeHtml(o.label) + "</option>";
    }).join("");

    for (const id of filterIds) {
      const def = FILTER_DEFS[id];
      if (!def) continue;
      filterBody.appendChild(renderFilterWidget(id, def, variants, state));
    }

    render();

    // ----- event delegation -----
    filterBody.addEventListener("change", function (ev) {
      const t = ev.target;
      const fid = t.getAttribute("data-filter-id");
      if (!fid) return;
      const def = FILTER_DEFS[fid];
      if (!def) return;

      if (def.type === "multi-checkbox") {
        if (!t.checked) state.checkboxes[fid].add(t.value);
        else            state.checkboxes[fid].delete(t.value);
      } else if (def.type === "min-range" || def.type === "max-range-novel-aware") {
        state.sliders[fid] = parseFloat(t.value);
        const lbl = t.parentElement.querySelector(".vb-slider-value");
        if (lbl) lbl.textContent = formatSliderValue(def, state.sliders[fid]);
      }
      state.expandedIdx = null;
      render();
    });

    filterBody.addEventListener("input", function (ev) {
      const t = ev.target;
      if (t.tagName !== "INPUT" || t.type !== "text") return;
      const fid = t.getAttribute("data-filter-id");
      if (!fid) return;
      state.gene = t.value.trim();
      state.expandedIdx = null;
      render();
    });

    filterBody.addEventListener("click", function (ev) {
      const btn = ev.target.closest(".vb-btn");
      if (!btn) return;
      const fid = btn.getAttribute("data-filter-id");
      const optId = btn.getAttribute("data-opt-id");
      if (!fid || !optId) return;
      state.buttonGroups[fid] = optId;
      // toggle active class on siblings
      const sibs = btn.parentElement.querySelectorAll(".vb-btn");
      sibs.forEach(function (b) { b.classList.remove("active"); });
      btn.classList.add("active");
      state.expandedIdx = null;
      render();
    });

    sortEl.addEventListener("change", function () {
      state.sortBy = sortEl.value;
      state.expandedIdx = null;
      render();
    });

    resetBtn.addEventListener("click", function () {
      for (const id of filterIds) {
        const def = FILTER_DEFS[id];
        if (!def) continue;
        if (def.type === "multi-checkbox") state.checkboxes[id] = new Set();
        else if (def.type === "min-range" || def.type === "max-range-novel-aware") state.sliders[id] = def.default;
        else if (def.type === "button-group") state.buttonGroups[id] = def.default;
      }
      state.gene = "";
      state.sortBy = "original";
      state.expandedIdx = null;
      filterBody.innerHTML = "";
      for (const id of filterIds) {
        const def = FILTER_DEFS[id];
        if (!def) continue;
        filterBody.appendChild(renderFilterWidget(id, def, variants, state));
      }
      sortEl.value = "original";
      render();
    });

    // List interactions: card expand, IGV chip click, IGV jump button, GeneBe link, MobiDetails copy+open
    listEl.addEventListener("click", function (ev) {
      // Report-triage buttons (Include / Exclude) -- clinical browser only.
      const triageBtn = ev.target.closest(".vb-triage-btn");
      if (triageBtn) {
        ev.stopPropagation(); // don't expand the card
        const igvKey = triageBtn.getAttribute("data-vb-key");
        const clicked = triageBtn.getAttribute("data-decision"); // "include" | "exclude"
        if (igvKey && window.tspipeReporting) {
          const row = variants.find(function (r) { return r._igvKey === igvKey; });
          if (row) {
            // CancerVar tier suggestion, when --annotate-cancervar was used at
            // build time. We carry it through the snapshot so the Reporting
            // table can render even if the variant TSV embed changes.
            const cvAnn = cancervarAnnotations[igvKey];
            const cvTier = (cvAnn && cvAnn.tier_label) ? cvAnn.tier_label : "";

            // Build the Reporting snapshot once at toggle-time so the table can
            // render even if the variant TSV embed changes in a future build.
            const snapshot = {
              gene:     row.Gene || "",
              hgvsG:    row.VV_HGVSg || row.HGVSg || "",
              hgvsP:    row.VV_HGVSp || row.HGVSp || "",
              hgvsC:    row.VV_HGVSc || row.HGVSc || "",
              exon:     row.EXON_REPORT || row.EXON || "",
              cosmic:   window.tspipeReporting.extractCosmicIds(row.Existing_variation),
              vaf:      row.VAF_pct || "",
              chr:      row.Chr || "",
              pos:      row.Start || "",
              ref:      row.Ref || "",
              alt:      row.Alt || "",
              cancervar_tier: cvTier,   // for audit; the live value lives in localStorage
            };

            // Clicking the already-active state clears to neutral; otherwise
            // switch to the clicked state. Include/Exclude are mutually exclusive.
            const current = window.tspipeReporting.getDecision("variant", igvKey);
            const next = (current === clicked) ? "" : clicked;
            window.tspipeReporting.setDecision("variant", igvKey, next, snapshot);

            // On a fresh Include with no manual tier yet, prefill the CancerVar
            // suggestion so the Reporting table shows it immediately.
            if (next === "include" && cvTier && !window.tspipeReporting.getTier(igvKey)) {
              window.tspipeReporting.setTier(igvKey, cvTier);
            }
          }
        }
        return;
      }
      // IGV chip on compact card
      const igvChip = ev.target.closest(".vb-igv-chip-ok");
      if (igvChip) {
        ev.stopPropagation();
        const uid = igvChip.getAttribute("data-igv-uid");
        if (uid !== null) jumpToIGV(uid);
        return;
      }
      // IGV jump button in detail view (same behaviour as chip)
      const igvBtn = ev.target.closest(".vb-igv-jump");
      if (igvBtn) {
        ev.stopPropagation();
        if (igvBtn.disabled) return;
        const uid = igvBtn.getAttribute("data-igv-uid");
        if (uid !== null) jumpToIGV(uid);
        return;
      }
      // Copy HGVS dropdown items
      const copyBtn = ev.target.closest(".vb-copy-hgvs");
      if (copyBtn) {
        // Intentionally don't stopPropagation — Bootstrap needs the click to
        // bubble so the dropdown auto-closes. The card-click branch below
        // already returns early for any button click.
        const text = copyBtn.getAttribute("data-hgvs");
        copyToClipboard(text, copyBtn);
        return;
      }
      // Ordinary card click → expand/collapse
      const card = ev.target.closest(".vb-card");
      if (!card) return;
      // Ignore clicks on links/buttons inside the detail view
      if (ev.target.closest("a") || ev.target.closest("button")) return;
      const idx = parseInt(card.getAttribute("data-idx"), 10);
      state.expandedIdx = (state.expandedIdx === idx) ? null : idx;
      render();
    });

    // -------------------------------------------------------------------------
    function render() {
      const filtered = variants.filter(function (r) { return passesAllFilters(r); });
      sortRows(filtered, state.sortBy);

      countEl.innerHTML =
        '<strong>' + filtered.length + '</strong> of ' + variants.length + ' variant' +
        (variants.length === 1 ? "" : "s") + " shown";

      if (filtered.length === 0) {
        listEl.innerHTML = '<div class="tspipe-empty">No variants match the current filters.</div>';
        return;
      }
      const html = [];
      for (const r of filtered) html.push(renderCard(r));
      listEl.innerHTML = html.join("");
    }

    function passesAllFilters(r) {
      // Checkbox filters
      for (const id of Object.keys(state.checkboxes)) {
        const def = FILTER_DEFS[id];
        const unchecked = state.checkboxes[id];
        if (unchecked.size === 0) continue;
        const cell = r[def.field];
        if (cell === undefined || cell === null || cell === "") continue;
        const tokens = def.splitOn ? String(cell).split(def.splitOn).map(function (s) { return s.trim(); })
                                   : [String(cell)];
        let any = false;
        for (const t of tokens) {
          if (!unchecked.has(t)) { any = true; break; }
        }
        if (!any) return false;
      }

      // Gene name search
      if (state.gene) {
        const g = (r.Gene || "").toUpperCase();
        if (g.indexOf(state.gene.toUpperCase()) === -1) return false;
      }

      // Slider filters
      for (const id of Object.keys(state.sliders)) {
        const def = FILTER_DEFS[id];
        const thr = state.sliders[id];
        const v = asNum(r[def.field]);
        if (def.type === "min-range") {
          if (v !== null && v < thr) return false;
        } else if (def.type === "max-range-novel-aware") {
          const raw = r[def.field];
          const isNovel = def.novelValues && def.novelValues.indexOf(String(raw || "")) !== -1;
          if (isNovel) continue;
          if (v !== null && v > thr) return false;
        }
      }

      // Button-group filters
      for (const id of Object.keys(state.buttonGroups)) {
        const def = FILTER_DEFS[id];
        const optId = state.buttonGroups[id];
        if (!optId || optId === "any") continue;
        const opt = def.options.find(function (o) { return o.id === optId; });
        if (!opt) continue;
        const v = asNum(r[def.field]);
        if (!opt.test(v)) return false;
      }

      return true;
    }

    function sortRows(rows, key) {
      const num = function (r, field) {
        const v = asNum(r[field]); return v === null ? -Infinity : v;
      };
      switch (key) {
        case "gene":         rows.sort(function (a, b) { return (a.Gene||"").localeCompare(b.Gene||""); }); break;
        case "vaf_desc":     rows.sort(function (a, b) { return num(b, "VAF_pct") - num(a, "VAF_pct"); }); break;
        case "vaf_asc":      rows.sort(function (a, b) { return num(a, "VAF_pct") - num(b, "VAF_pct"); }); break;
        case "alt_desc":     rows.sort(function (a, b) { return num(b, "ALT_COUNT") - num(a, "ALT_COUNT"); }); break;
        case "callers_desc": rows.sort(function (a, b) { return num(b, "VariantCaller_Count") - num(a, "VariantCaller_Count"); }); break;
        default:             rows.sort(function (a, b) { return a._idx - b._idx; });
      }
    }

    function renderCard(r) {
      const idx = r._idx;
      const isOpen = state.expandedIdx === idx;
      const verdict = r.SomaticSeq_Verdict || "";
      const verdictBadge = verdict
        ? '<span class="badge ' + verdictBadgeClass(verdict) + '">' + escapeHtml(verdict) + "</span>"
        : "";
      const hgvsp = bestHGVSp(r);
      const conseq = r.Consequence || "";
      const callerCount = r.VariantCaller_Count || "";
      const refCount = r.REF_COUNT || "0";
      const altCount = r.ALT_COUNT || "0";
      const vaf = r.VAF_pct || "";

      const filterLabel = (r.Filter && r.Filter !== "PASS")
        ? ' <span class="badge bg-light text-dark border ms-1">' + escapeHtml(r.Filter) + "</span>"
        : "";

      // IGV chip on the compact view — visible without expanding the card.
      let igvChip = "";
      if (hasIGV) {
        const uid = igvLookup[r._igvKey];
        if (uid !== undefined) {
          igvChip = '<span class="vb-igv-chip vb-igv-chip-ok" data-igv-uid="' + uid +
                    '" title="Open this variant in IGV">IGV \u2192</span>';
        } else {
          igvChip = '<span class="vb-igv-chip vb-igv-chip-missing" title="Not in IGV report">no IGV</span>';
        }
      }

      // Report-triage control -- only rendered on the clinical browser. Three
      // states: neutral, Include, Exclude (Include/Exclude mutually exclusive).
      // State is read fresh from the store on every render so it stays in sync
      // when changed from another tab (e.g. Clear all on the Reporting tab).
      let selectCheckbox = "";
      if (enableReportSelect) {
        const decision = window.tspipeReporting ? window.tspipeReporting.getDecision("variant", r._igvKey) : "";
        const key = escapeHtml(r._igvKey);
        selectCheckbox =
          '<div class="btn-group-vertical btn-group-sm vb-triage pt-1 me-1" role="group" ' +
               'aria-label="Report triage">' +
            '<button type="button" class="btn btn-outline-success vb-triage-btn' +
                   (decision === "include" ? " active" : "") + '" ' +
                   'data-vb-key="' + key + '" data-decision="include" ' +
                   'title="Include in report">Incl</button>' +
            '<button type="button" class="btn btn-outline-danger vb-triage-btn' +
                   (decision === "exclude" ? " active" : "") + '" ' +
                   'data-vb-key="' + key + '" data-decision="exclude" ' +
                   'title="Exclude (reviewed, not reported)">Excl</button>' +
          "</div>";
      }

      const compact =
        '<div class="vb-card-compact d-flex justify-content-between gap-2">' +
          selectCheckbox +
          '<div class="flex-grow-1">' +
            '<div class="d-flex align-items-baseline gap-2 flex-wrap">' +
              '<span class="vb-gene">' + escapeHtml(r.Gene || "?") + "</span>" +
              '<span class="text-muted small">' + escapeHtml(conseq) + "</span>" +
              (r.IMPACT ? '<span class="badge bg-light text-dark border">' + escapeHtml(r.IMPACT) + "</span>" : "") +
              igvChip +
            "</div>" +
            '<div class="small font-monospace text-muted mt-1">' +
              escapeHtml(hgvsp || (r.HGVSc || "")) +
              ' \u00b7 ' +
              escapeHtml(r.Chr || "") + ":" + escapeHtml(r.Start || "") + " " +
              escapeHtml(r.Ref || "") + "&gt;" + escapeHtml(r.Alt || "") +
            "</div>" +
          "</div>" +
          '<div class="text-end small">' +
            verdictBadge + filterLabel +
            '<div class="mt-1">VAF ' + escapeHtml(vaf) + "%</div>" +
            '<div class="text-muted">' + escapeHtml(altCount) + " ALT / " + escapeHtml(refCount) + " REF</div>" +
            '<div class="text-muted">' + escapeHtml(callerCount) + " caller" + (callerCount === "1" ? "" : "s") + "</div>" +
          "</div>" +
        "</div>";

      const detail = isOpen ? renderDetail(r) : "";

      return '<div class="vb-card list-group-item list-group-item-action' +
             (verdict === "REJECT" ? " vb-card-reject" : "") +
             (isOpen ? " vb-card-open" : "") +
             '" data-idx="' + idx + '" role="button" tabindex="0">' +
             compact + detail + "</div>";
    }

    function renderDetail(r) {
      const groups = [];
      for (const [heading, fields] of DETAIL_GROUPS) {
        const rowsHtml = [];
        for (const [field, label] of fields) {
          if (!(field in r)) continue;
          const v = r[field];
          rowsHtml.push(
            '<dt class="col-sm-5 text-muted fw-normal small">' + escapeHtml(label) + "</dt>" +
            '<dd class="col-sm-7 small mb-1' +
              (field === "HGVSc" || field === "HGVSg" || field === "VV_HGVSc" || field === "VV_HGVSg" ||
               field === "HGVSp" || field === "VV_HGVSp" ? ' font-monospace' : '') +
            '">' + emptyVal(v) + "</dd>"
          );
        }
        if (rowsHtml.length === 0) continue;
        groups.push(
          '<div class="col-md-6">' +
            '<h6 class="text-uppercase text-muted small mt-3 mb-2">' + escapeHtml(heading) + "</h6>" +
            '<dl class="row mb-0">' + rowsHtml.join("") + "</dl>" +
          "</div>"
        );
      }

      // External-link buttons (always visible)
      const chrClean = String(r.Chr || "").replace(/^chr/, "");
      const gbUrl = genebeUrl(chrClean, r.Start, r.Ref, r.Alt);

      let extButtons =
        '<a href="' + escapeHtml(gbUrl) + '" target="_blank" rel="noopener" class="btn btn-sm btn-outline-primary">' +
          "Open in GeneBe \u2197</a>";

      // Copy VV_HGVS dropdown -- c / p / g. Items for empty/-1 fields are
      // rendered disabled so the user can see which projections are available.
      // We intentionally do NOT render an "Open MobiDetails" button: MD has no
      // reliable anonymous deep link (the /api/variant/{id}/browser/ form
      // shipped in v0.3.1 triggered a "No API key provided" modal client-side,
      // and the /create_g?caller=browser fallback just redirects to MD home
      // without auth). The Copy VV_HGVS dropdown gives the user the HGVS
      // string in any of three projections so they can paste it into MD's
      // search box themselves.
      const hgvsItems = [
        { code: "C", field: "VV_HGVSc", label: "Copy VV_HGVS_C" },
        { code: "P", field: "VV_HGVSp", label: "Copy VV_HGVS_P" },
        { code: "G", field: "VV_HGVSg", label: "Copy VV_HGVS_G" },
      ];
      // Effective exon for copy: VV_Exon -> panel-BED -> VEP EXON, precomputed
      // at build time into EXON_REPORT (falls back to raw EXON for older builds).
      const exonVal = String(r.EXON_REPORT || r.EXON || "").trim();
      const hasExon = hasUsefulValue(exonVal);
      const anyHgvs = hgvsItems.some(function (it) { return hasUsefulValue(r[it.field]); });
      if (anyHgvs || hasExon) {
        let itemsHtml = hgvsItems.map(function (it) {
          const v = r[it.field];
          if (hasUsefulValue(v)) {
            return '<li><button type="button" class="dropdown-item vb-copy-hgvs" ' +
                     'data-hgvs="' + escapeHtml(String(v)) + '" ' +
                     'title="' + escapeHtml(String(v)) + '">' +
                     escapeHtml(it.label) + "</button></li>";
          }
          return '<li><button type="button" class="dropdown-item disabled" disabled ' +
                   'title="' + escapeHtml(it.field) + ' not available">' +
                   escapeHtml(it.label) + " (n/a)</button></li>";
        }).join("");
        // Exon item (copies the bare exon token).
        itemsHtml += '<li><hr class="dropdown-divider"></li>';
        if (hasExon) {
          itemsHtml += '<li><button type="button" class="dropdown-item vb-copy-hgvs" ' +
                         'data-hgvs="' + escapeHtml(exonVal) + '" ' +
                         'title="Exon ' + escapeHtml(exonVal) + '">Copy Exon</button></li>';
        } else {
          itemsHtml += '<li><button type="button" class="dropdown-item disabled" disabled ' +
                         'title="Exon not available">Copy Exon (n/a)</button></li>';
        }
        extButtons +=
          '<div class="dropdown d-inline-block">' +
            '<button type="button" class="btn btn-sm btn-outline-primary dropdown-toggle" ' +
              'data-bs-toggle="dropdown" aria-expanded="false" ' +
              'title="Copy VariantValidator HGVS to clipboard">Copy VV_HGVS</button>' +
            '<ul class="dropdown-menu dropdown-menu-end">' + itemsHtml + "</ul>" +
          "</div>";
      }

      // IGV button (clinical tab only)
      let igvButton = "";
      if (hasIGV) {
        const uid = igvLookup[r._igvKey];
        if (uid !== undefined) {
          igvButton =
            '<button type="button" class="btn btn-sm btn-primary vb-igv-jump" data-igv-uid="' + uid + '">' +
              "View in IGV \u2192</button>";
        } else {
          igvButton =
            '<button type="button" class="btn btn-sm btn-outline-secondary vb-igv-jump" disabled ' +
              'title="This variant is not in the IGV report">Not in IGV report</button>';
        }
      }

      // GeneBe annotation block (only present when --annotate-genebe was used at build time)
      let genebeBlock = "";
      const ann = genebeAnnotations[r._igvKey];
      if (ann) {
        genebeBlock = renderGenebeBlock(ann);
      }

      // CancerVar annotation block (only present when --annotate-cancervar was used at build time)
      let cancervarBlock = "";
      const cvAnn = cancervarAnnotations[r._igvKey];
      if (cvAnn) {
        cancervarBlock = renderCancervarBlock(cvAnn);
      }

      // OncoKB annotation block (only present when --annotate-oncokb was used at build time)
      let oncokbBlock = "";
      const oncoAnn = oncokbAnnotations[r._igvKey];
      if (oncoAnn) {
        oncokbBlock = renderOncokbBlock(oncoAnn);
      }

      return '<div class="vb-card-detail mt-3 pt-3 border-top">' +
               '<div class="row g-3">' + groups.join("") + "</div>" +
               genebeBlock +
               cancervarBlock +
               oncokbBlock +
               '<div class="vb-ext-links mt-3 justify-content-end">' +
                  extButtons + igvButton +
                  '<span class="vb-copied-toast">Copied!</span>' +
               '</div>' +
             "</div>";
    }

    function renderGenebeBlock(ann) {
      const acmgPill = ann.acmg_classification
        ? '<span class="vb-acmg-pill ' + acmgPillClass(ann.acmg_classification) + '">' +
            escapeHtml(ann.acmg_classification) + "</span>"
        : '<span class="text-muted small">not provided</span>';

      const clinvarHtml = ann.clinvar_classification
        ? escapeHtml(ann.clinvar_classification) +
          (ann.clinvar_review_status ? ' <span class="text-muted small">(' + escapeHtml(ann.clinvar_review_status) + ')</span>' : '')
        : '<span class="text-muted">\u2014</span>';

      const clinvarDiseaseHtml = ann.clinvar_disease
        ? escapeHtml(ann.clinvar_disease)
        : '<span class="text-muted">\u2014</span>';

      const acmgCritHtml = ann.acmg_criteria
        ? escapeHtml(ann.acmg_criteria) +
          (ann.acmg_score !== undefined && ann.acmg_score !== null
              ? ' <span class="text-muted small">(score ' + escapeHtml(String(ann.acmg_score)) + ')</span>'
              : '')
        : '<span class="text-muted">\u2014</span>';

      const fmtNum = function (v, digits) {
        if (v === undefined || v === null || v === "") return '<span class="text-muted">\u2014</span>';
        const n = parseFloat(v);
        return isNaN(n) ? escapeHtml(String(v)) : n.toFixed(digits);
      };

      const rows = [
        ["ACMG classification",  acmgPill],
        ["ACMG criteria",        acmgCritHtml],
        ["ClinVar",              clinvarHtml],
        ["ClinVar disease",      clinvarDiseaseHtml],
        ["gnomAD exome AF",      fmtNum(ann.gnomad_exome_af, 6)],
        ["gnomAD genome AF",     fmtNum(ann.gnomad_genome_af, 6)],
        ["REVEL",                ann.revel_score !== undefined && ann.revel_score !== null
                                    ? fmtNum(ann.revel_score, 3)
                                    : '<span class="text-muted">\u2014</span>'],
        ["AlphaMissense",        ann.alphamissense_prediction
                                    ? escapeHtml(String(ann.alphamissense_prediction))
                                    : '<span class="text-muted">\u2014</span>'],
        ["SpliceAI max",         ann.spliceai_max_score !== undefined && ann.spliceai_max_score !== null
                                    ? fmtNum(ann.spliceai_max_score, 3)
                                    : '<span class="text-muted">\u2014</span>'],
      ];

      const rowsHtml = rows.map(function (kv) {
        return '<div class="row g-2 mb-1">' +
                 '<div class="col-sm-4 label">' + kv[0] + '</div>' +
                 '<div class="col-sm-8">' + kv[1] + '</div>' +
               '</div>';
      }).join("");

      const fetchedAt = ann._fetched_at
        ? '<div class="text-muted small mt-2">Fetched ' + escapeHtml(ann._fetched_at) + ' from GeneBe</div>'
        : "";

      return '<div class="vb-genebe-block mt-3">' +
               '<h6 class="text-uppercase text-muted small mb-2">External annotation (GeneBe, build-time)</h6>' +
               rowsHtml + fetchedAt +
             "</div>";
    }

    function renderCancervarBlock(ann) {
      // Tier pill -- map AMP/ASCO/CAP tier to the same ACMG pill colour palette
      // so the visual encoding stays consistent across the three external blocks.
      //   Tier I  -> strong evidence -> "pathogenic" red
      //   Tier II -> potential       -> "likely pathogenic" orange
      //   Tier III -> uncertain      -> "uncertain" gray
      //   Tier IV -> benign          -> "benign" green
      // Match on the leading Roman-numeral root so capitalisation of the third
      // token in the raw slug (e.g. "Tier_III_Uncertain" vs "Tier_III_unknown")
      // doesn't flip the colour.
      function tierPillClass(slug) {
        if (!slug) return "vb-acmg-uncertain";
        if (/^Tier_I(_|$)/.test(slug))   return "vb-acmg-pathogenic";
        if (/^Tier_II(_|$)/.test(slug))  return "vb-acmg-likely-pathogenic";
        if (/^Tier_III(_|$)/.test(slug)) return "vb-acmg-uncertain";
        if (/^Tier_IV(_|$)/.test(slug))  return "vb-acmg-benign";
        return "vb-acmg-uncertain";
      }

      const tierLabel = ann.tier_label || "\u2014";
      const tierPill = '<span class="vb-acmg-pill ' + tierPillClass(ann.tier_slug) + '">' +
                          escapeHtml(tierLabel) + "</span>";
      const tierDesc = ann.tier_description
        ? ' <span class="text-muted small ms-1">' + escapeHtml(ann.tier_description) + "</span>"
        : "";

      // OPAI score -- 0..1, higher = more likely oncogenic. We show it as a
      // formatted decimal alongside a colour cue.
      let opaiHtml = '<span class="text-muted">\u2014</span>';
      if (ann.opai !== null && ann.opai !== undefined) {
        const opaiNum = Number(ann.opai);
        let opaiBadge = "bg-secondary";
        if (!isNaN(opaiNum)) {
          if (opaiNum >= 0.8)      opaiBadge = "bg-danger";
          else if (opaiNum >= 0.5) opaiBadge = "bg-warning text-dark";
          else                     opaiBadge = "bg-success";
        }
        opaiHtml = '<span class="badge ' + opaiBadge + '">' +
                     (isNaN(opaiNum) ? escapeHtml(String(ann.opai)) : opaiNum.toFixed(2)) +
                   "</span>" +
                   ' <span class="text-muted small ms-1">' +
                     "(Oncogenic Pathogenicity AI, 0\u20131)" +
                   "</span>";
      }

      // CBP criteria table -- 12 rows, score colour-coded:
      //   0 = not met (muted), 1 = supporting (blue), 2 = strong (orange).
      // Hidden behind a <details> element to keep the panel compact by default.
      function cbpBadge(score) {
        if (score === 2) return '<span class="badge bg-warning text-dark">2</span>';
        if (score === 1) return '<span class="badge bg-info text-dark">1</span>';
        return '<span class="badge bg-light text-muted border">0</span>';
      }

      let cbpHtml = "";
      const cbp = ann.cbp || [];
      if (cbp.length) {
        const cbpRows = cbp.map(function (c) {
          return "<tr>" +
                   '<td class="font-monospace small">' + escapeHtml(c.id) + "</td>" +
                   '<td class="text-center">' + cbpBadge(c.score) + "</td>" +
                   "<td class=\"small\">" + escapeHtml(c.description || "") + "</td>" +
                 "</tr>";
        }).join("");
        cbpHtml = '<details class="mt-2">' +
                    '<summary class="text-muted small" style="cursor:pointer">' +
                      "Show 12 CBP criteria scores" +
                    "</summary>" +
                    '<table class="table table-sm table-borderless mt-2 mb-0">' +
                      "<thead><tr>" +
                        '<th scope="col" class="small">Criterion</th>' +
                        '<th scope="col" class="small text-center" style="width:60px">Score</th>' +
                        '<th scope="col" class="small">Description</th>' +
                      "</tr></thead>" +
                      "<tbody>" + cbpRows + "</tbody>" +
                    "</table>" +
                  "</details>";
      }

      const rows = [
        ["AMP/ASCO/CAP Tier", tierPill + tierDesc],
        ["OPAI",              opaiHtml],
        ["Gene",              ann.gene ? escapeHtml(ann.gene) : '<span class="text-muted">\u2014</span>'],
      ];

      const rowsHtml = rows.map(function (kv) {
        return '<div class="row g-2 mb-1">' +
                 '<div class="col-sm-4 label">' + kv[0] + '</div>' +
                 '<div class="col-sm-8">' + kv[1] + '</div>' +
               '</div>';
      }).join("");

      const fetchedAt = ann._fetched_at
        ? '<div class="text-muted small mt-2">Fetched ' + escapeHtml(ann._fetched_at) + ' from CancerVar' +
            (ann.build ? " (" + escapeHtml(ann.build) + ")" : "") + "</div>"
        : "";

      return '<div class="vb-genebe-block mt-3">' +
               '<h6 class="text-uppercase text-muted small mb-2">External annotation (CancerVar, build-time)</h6>' +
               rowsHtml + cbpHtml + fetchedAt +
             "</div>";
    }

    function renderOncokbBlock(ann) {
      // Oncogenic status pill -- reuse the ACMG pill styles by approximate
      // semantic match (Oncogenic ~ Pathogenic, Likely Oncogenic ~ Likely Path).
      function oncogenicPillClass(status) {
        if (!status) return "vb-acmg-uncertain";
        const s = status.toLowerCase();
        if (s.indexOf("likely oncogenic") !== -1) return "vb-acmg-likely-pathogenic";
        if (s.indexOf("oncogenic") !== -1) return "vb-acmg-pathogenic";
        if (s.indexOf("likely neutral") !== -1) return "vb-acmg-likely-benign";
        if (s.indexOf("neutral") !== -1) return "vb-acmg-benign";
        return "vb-acmg-uncertain";
      }

      const oncoPill = ann.oncogenic
        ? '<span class="vb-acmg-pill ' + oncogenicPillClass(ann.oncogenic) + '">' +
            escapeHtml(ann.oncogenic) + "</span>"
        : '<span class="text-muted small">not provided</span>';

      const hotspotHtml = ann.hotspot
        ? '<span class="badge bg-warning text-dark">Hotspot</span>'
        : (ann.hotspot === false ? '<span class="text-muted">not a hotspot</span>' : '<span class="text-muted">\u2014</span>');

      const mutEffect = ann.mutation_effect_known
        ? escapeHtml(ann.mutation_effect_known)
        : '<span class="text-muted">\u2014</span>';

      // Highest evidence levels -- show whichever are present.
      const levels = [];
      if (ann.highest_sensitive_level)  levels.push("Sensitivity: " + escapeHtml(ann.highest_sensitive_level));
      if (ann.highest_resistance_level) levels.push("Resistance: " + escapeHtml(ann.highest_resistance_level));
      if (ann.highest_diagnostic_implication_level) levels.push("Diagnostic: " + escapeHtml(ann.highest_diagnostic_implication_level));
      if (ann.highest_prognostic_implication_level) levels.push("Prognostic: " + escapeHtml(ann.highest_prognostic_implication_level));
      const levelsHtml = levels.length
        ? levels.join(" &middot; ")
        : '<span class="text-muted">no FDA/CLIA level</span>';

      // Treatments table -- compact, capped to top 5 to avoid the detail
      // panel ballooning.
      let treatmentsHtml = '<span class="text-muted">\u2014</span>';
      const trts = ann.treatments || [];
      if (trts.length) {
        const shown = trts.slice(0, 5);
        const rows = shown.map(function (t) {
          const drugs = (t.drugs || []).join(" + ");
          const lvl = t.level ? '<span class="badge bg-light text-dark border ms-1">' + escapeHtml(t.level) + "</span>" : "";
          const ct = t.cancer_type ? ' <span class="text-muted small">(' + escapeHtml(t.cancer_type) + ")</span>" : "";
          return "<li>" + escapeHtml(drugs) + lvl + ct + "</li>";
        }).join("");
        const more = trts.length > shown.length
          ? '<li class="text-muted small">\u2026 ' + (trts.length - shown.length) + " more not shown</li>"
          : "";
        treatmentsHtml = '<ul class="mb-0 ps-3">' + rows + more + "</ul>";
      }

      // Optional prose summaries.
      const variantSummary = ann.variant_summary
        ? '<div class="text-muted small mt-2">' + escapeHtml(ann.variant_summary) + "</div>"
        : "";

      const rows = [
        ["Oncogenic",       oncoPill],
        ["Mutation effect", mutEffect],
        ["Hotspot",         hotspotHtml],
        ["Evidence level",  levelsHtml],
        ["Treatments",      treatmentsHtml],
      ];

      const rowsHtml = rows.map(function (kv) {
        return '<div class="row g-2 mb-1">' +
                 '<div class="col-sm-4 label">' + kv[0] + '</div>' +
                 '<div class="col-sm-8">' + kv[1] + '</div>' +
               '</div>';
      }).join("");

      const fetchedAt = ann._fetched_at
        ? '<div class="text-muted small mt-2">Fetched ' + escapeHtml(ann._fetched_at) + ' from OncoKB</div>'
        : "";

      return '<div class="vb-genebe-block mt-3">' +
               '<h6 class="text-uppercase text-muted small mb-2">External annotation (OncoKB, build-time)</h6>' +
               rowsHtml + variantSummary + fetchedAt +
             "</div>";
    }

    // -------------------------------------------------------------------------
    function jumpToIGV(uniqueId) {
      const frame = document.getElementById(config.igvFrameId);
      if (!frame) return;

      // Activate the IGV tab first so the iframe is visible (and so the
      // lazy-load handler's data-src -> src swap has a chance to fire).
      if (config.igvTabSelector) {
        const trigger = document.querySelector(config.igvTabSelector);
        if (trigger && window.bootstrap && window.bootstrap.Tab) {
          window.bootstrap.Tab.getOrCreateInstance(trigger).show();
        }
      }

      // Determine the iframe's base URL (path without fragment) -- prefer the
      // data-src value if the iframe hasn't been lazy-loaded yet, else strip
      // any existing #fragment from src.
      const baseSrc = frame.getAttribute("data-src") ||
                      (frame.src || "").replace(/#.*$/, "");
      if (!baseSrc) return;

      const newSrc = baseSrc + "#row_" + uniqueId;

      // Navigate the iframe to the new fragment. This works cross-origin
      // (including under file://) because the parent is only writing src --
      // never reading the iframe's DOM. The hash-router script injected into
      // the IGV report (by build.py via parsers/igv.py) sees the new hash
      // and clicks the matching row from inside the iframe.
      //
      // Browser behaviour:
      //   - First time: data-src is still set, iframe was never loaded;
      //     setting src triggers the initial load with the hash already in
      //     place, so the hash-router fires once IGV.js is ready.
      //   - Subsequent clicks: if only the hash changed, the iframe doesn't
      //     reload -- the hashchange event fires inside, and the router
      //     re-runs trySelect() to click the new row.
      //   - Re-clicking the same row: setting src to its current value is a
      //     no-op in some browsers (no event). The router handles this case
      //     because the row is already selected from the previous click;
      //     nothing more is needed.
      if (frame.hasAttribute("data-src")) {
        frame.src = newSrc;
        frame.removeAttribute("data-src");
      } else {
        frame.src = newSrc;
      }
    }

    function copyToClipboard(text, btn) {
      // Toast lives on .vb-ext-links — climb out of any wrapper (e.g. dropdown <li>).
      const container = btn.closest(".vb-ext-links") || btn.parentElement;
      const toast = container ? container.querySelector(".vb-copied-toast") : null;
      const showToast = function () {
        if (!toast) return;
        toast.classList.add("visible");
        setTimeout(function () { toast.classList.remove("visible"); }, 1500);
      };
      if (navigator.clipboard && window.isSecureContext) {
        navigator.clipboard.writeText(text).then(showToast).catch(function () {
          fallbackCopy(text); showToast();
        });
      } else {
        fallbackCopy(text); showToast();
      }
    }

    function fallbackCopy(text) {
      const ta = document.createElement("textarea");
      ta.value = text;
      ta.style.position = "fixed";
      ta.style.opacity = "0";
      document.body.appendChild(ta);
      ta.select();
      try { document.execCommand("copy"); } catch (e) { /* ignore */ }
      document.body.removeChild(ta);
    }

  }; // end initVariantBrowser

  // ---------------------------------------------------------------------------
  // Widget builders
  // ---------------------------------------------------------------------------
  function renderFilterWidget(id, def, variants, state) {
    const wrapper = document.createElement("div");
    wrapper.className = "vb-filter mb-3";

    const label = document.createElement("label");
    label.className = "form-label small fw-semibold mb-1";
    label.textContent = def.label;
    wrapper.appendChild(label);

    if (def.type === "multi-checkbox") {
      let values = uniqueValues(variants, def.field, def.splitOn);
      if (def.preferredOrder) {
        const byKey = new Map(values);
        const ordered = [];
        for (const k of def.preferredOrder) if (byKey.has(k)) ordered.push([k, byKey.get(k)]);
        for (const [k, v] of values) if (def.preferredOrder.indexOf(k) === -1) ordered.push([k, v]);
        values = ordered;
      }
      const box = document.createElement("div");
      box.className = "vb-checkbox-group" + (def.scroll ? " vb-checkbox-scroll" : "");
      if (values.length === 0) {
        box.innerHTML = '<div class="text-muted small">(no values)</div>';
      } else {
        for (const [val, count] of values) {
          const item = document.createElement("div");
          item.className = "form-check small";
          const inputId = "cb_" + id + "_" + Math.random().toString(36).slice(2, 8);
          item.innerHTML =
            '<input class="form-check-input" type="checkbox" checked ' +
              'id="' + inputId + '" data-filter-id="' + id + '" value="' + escapeHtml(val) + '">' +
            '<label class="form-check-label" for="' + inputId + '">' +
              escapeHtml(val) + ' <span class="text-muted">(' + count + ")</span>" +
            "</label>";
          box.appendChild(item);
        }
      }
      wrapper.appendChild(box);
    } else if (def.type === "text") {
      const input = document.createElement("input");
      input.type = "text";
      input.className = "form-control form-control-sm";
      input.setAttribute("data-filter-id", id);
      input.placeholder = "Type to filter\u2026";
      input.value = state.gene || "";
      wrapper.appendChild(input);
    } else if (def.type === "min-range" || def.type === "max-range-novel-aware") {
      const cur = state.sliders[id];
      const range = document.createElement("input");
      range.type = "range";
      range.className = "form-range";
      range.min = def.min;
      range.max = def.max;
      range.step = def.step;
      range.value = cur;
      range.setAttribute("data-filter-id", id);
      const valLabel = document.createElement("div");
      valLabel.className = "vb-slider-value small text-muted";
      valLabel.textContent = formatSliderValue(def, cur);
      wrapper.appendChild(range);
      wrapper.appendChild(valLabel);
    } else if (def.type === "button-group") {
      const group = document.createElement("div");
      group.className = "vb-btn-group";
      const active = state.buttonGroups[id];
      for (const opt of def.options) {
        const b = document.createElement("button");
        b.type = "button";
        b.className = "vb-btn" + (opt.id === active ? " active" : "");
        b.setAttribute("data-filter-id", id);
        b.setAttribute("data-opt-id", opt.id);
        b.textContent = opt.label;
        group.appendChild(b);
      }
      wrapper.appendChild(group);
    }
    return wrapper;
  }

  function formatSliderValue(def, v) {
    if (def.field === "Max_AF") {
      if (v >= 1) return "1.0 (off)";
      return Number(v).toFixed(4);
    }
    if (def.field === "VAF_pct") return Number(v).toFixed(1) + " %";
    return String(v);
  }
})();
