"""
Build the unified acute leukemia / myeloid neoplasm adaptive-sampling BED (hg38).

    python build_panel.py --targets panel_targets.tsv --out AL_unified_hg38.bed

Two panels, each built from its own target table: AML_panel_targets.tsv
(myeloid, including MPN / M-LN-Eo / MDS-MPN) and ALL_panel_targets.tsv (B-ALL
and T-ALL). Lineage is assigned by SPARSH before this reflex assay fires, so
the two never share a flow cell. The SCOPE column records whether a target is
specific to that panel or shared with the other; rows with KEEP != Y are
skipped. Within a panel a gene is covered once however many groups nominated it. Emits BED4 (chrom, start, end, name), 0-based
half-open, coordinate sorted. Overlapping regions merge and their names join
with "/", matching the AS_240827 panel convention.

Coverage rules
--------------
Default        : whole gene body plus FLANK either side.
CENTERED genes : gene-body centre plus the stated half-width, for loci whose
                 breakpoints lie outside the gene body. CENTERED_BY_PANEL
                 overrides the half-width for one panel, so a SCOPE=BOTH gene
                 can be widened on one panel while the other is mid-run.
NAMED_REGIONS  : fixed intervals that are not a single RefSeq gene.
"""

import argparse
import gzip
import sys
from collections import defaultdict

PRIMARY = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}

FLANK = 50000

# Gene-body centre plus half-width. Breakpoints for these sit outside the gene.
CENTERED = {
    "TP53":   500000,   # del(17p) breakpoint resolution
    "MECOM": 1000000,   # inv(3)/t(3;3) breakpoints scatter over 1-2 Mb at 3q26
    "BCL11B": 500000,   # 14q32 enhancer hijacking, ~700 kb from the gene
    "MYC":   1000000,   # 8q24 enhancer hijacking; IGH::MYC in B-ALL and
                        # MYC x TR loci in T-ALL both break outside the gene.
                        # A gene-body window misses them: on an MM sample with
                        # 61 chr8 breakpoints between 126.1 and 131.9 Mb, none
                        # fell inside the 107 kb gene-body interval. The MM v7
                        # panel uses plus or minus 2.5 Mb; 1 Mb is the AL
                        # compromise, matching MECOM, pending AL data.
}

# Per-panel overrides. A gene present here uses the panel's value instead of
# the shared one. A value of None means "not centred": fall back to the normal
# gene-body-plus-FLANK rule.
#
# This exists because several genes are SCOPE=BOTH, so a change to CENTERED
# hits both panels at once. When one panel is mid-sequencing, its BED must not
# move even though the other panel's should.
CENTERED_BY_PANEL = {
    # "AML": {"MYC": None},   # hold AML at the default flank while the
    #                         # current batch is still on the instrument
}


def centered_halfwidth(name, panel):
    """Half-width for a centred target, or None if it is not centred."""
    per_panel = CENTERED_BY_PANEL.get((panel or "").strip().upper(), {})
    if name in per_panel:
        return per_panel[name]          # may be None, meaning not centred
    return CENTERED.get(name)

# Fixed intervals that are not single RefSeq genes.
NAMED_REGIONS = {
    "IGH_locus":              ("chr14", 105550000, 106900000),
    "IGK_locus":              ("chr2",   88830000,  90280000),
    "IGL_locus":              ("chr22",  22020000,  22930000),
    "TRA_TRD_locus":          ("chr14",  21620000,  22560000),
    "TRB_locus":              ("chr7",  142290000, 142820000),
    "TRG_locus":              ("chr7",   38230000,  38400000),
    "PAR1_CRLF2_P2RY8":       ("chrX",    1150000,   1650000),
    "PAR1_CRLF2_P2RY8_Y":     ("chrY",    1150000,   1650000),
    "HOXA_cluster":           ("chr7",   27090000,  27300000),
    "HOXB_cluster":           ("chr17",  48520000,  48630000),
}

# 17p CN-LOH tiling was evaluated and removed: at the on-target depth this
# panel supports, per-tile BAF could not be resolved above binomial noise.
# Copy-neutral LOH at 17p is therefore not addressed; the TP53 centred window
# resolves del(17p) breakpoints only. Under the CAM 2026 algorithm, cases
# without locus assessment fall to Tier B rather than Tier A.
TILE_SPEC = {}


def load_refgene(path):
    """Longest transcript per gene symbol, primary chromosomes only."""
    best = {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            chrom = f[2]
            if chrom not in PRIMARY:
                continue
            gene, tx_s, tx_e = f[12], int(f[4]), int(f[5])
            if gene not in best or (tx_e - tx_s) > (best[gene][2] - best[gene][1]):
                best[gene] = (chrom, tx_s, tx_e)
    return best


def load_sizes(path):
    with open(path) as fh:
        return {l.split()[0]: int(l.split()[1]) for l in fh}


def clip(s, e, chrom, sizes):
    return max(0, s), min(e, sizes.get(chrom, e))


def merge(records):
    by_chrom = defaultdict(list)
    for chrom, s, e, name in records:
        by_chrom[chrom].append((s, e, name))

    def chrom_key(c):
        c2 = c.replace("chr", "")
        return (0, int(c2)) if c2.isdigit() else (1, c2)

    out = []
    for chrom in sorted(by_chrom, key=chrom_key):
        iv = sorted(by_chrom[chrom])
        cs, ce, names = iv[0][0], iv[0][1], [iv[0][2]]
        for s, e, name in iv[1:]:
            if s <= ce:
                ce = max(ce, e)
                if name not in names:
                    names.append(name)
            else:
                out.append((chrom, cs, ce, "/".join(names)))
                cs, ce, names = s, e, [name]
        out.append((chrom, cs, ce, "/".join(names)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", required=True,
                    help="AML_panel_targets.tsv or ALL_panel_targets.tsv")
    ap.add_argument("--refgene", default="refGene.txt.gz")
    ap.add_argument("--sizes", default="hg38.chrom.sizes")
    ap.add_argument("--out", default="AL_unified_hg38.bed")

    args = ap.parse_args()

    rg = load_refgene(args.refgene)
    sizes = load_sizes(args.sizes)

    # Read the approved target list; a gene nominated under several groups is
    # emitted once, keeping the group that first nominated it.
    targets, order = {}, []
    with open(args.targets) as fh:
        header = fh.readline()
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            keep, panel, _scope, group, name = line.split("\t")[:5]
            if keep.strip().upper() != "Y":
                continue
            if name not in targets:
                targets[name] = (group, panel)
                order.append(name)

    records, missing, counts = [], [], defaultdict(int)

    for name in order:
        group, panel = targets[name]

        if name in TILE_SPEC:
            chrom, start, ntile, width, spacing = TILE_SPEC[name]
            for i in range(ntile):
                pos = start + i * spacing
                s, e = clip(pos, pos + width, chrom, sizes)
                records.append((chrom, s, e, f"17p_LOH_t{i:02d}"))
                counts[group] += e - s
            continue

        if name in NAMED_REGIONS:
            chrom, s, e = NAMED_REGIONS[name]
            s, e = clip(s, e, chrom, sizes)
            records.append((chrom, s, e, name))
            counts[group] += e - s
            continue

        if name not in rg:
            missing.append(name)
            continue

        chrom, tx_s, tx_e = rg[name]
        hw = centered_halfwidth(name, panel)
        if hw is not None:
            mid = (tx_s + tx_e) // 2
            s, e = clip(mid - hw, mid + hw, chrom, sizes)
        else:
            s, e = clip(tx_s - FLANK, tx_e + FLANK, chrom, sizes)
        records.append((chrom, s, e, name))
        counts[group] += e - s

    final = merge(records)

    with open(args.out, "w") as fh:
        for chrom, s, e, name in final:
            fh.write(f"{chrom}\t{s}\t{e}\t{name}\n")

    total = sum(e - s for _, s, e, _ in final)
    genome = sum(sizes[c] for c in PRIMARY if c in sizes)
    sys.stderr.write(
        f"targets in list : {len(order)}\n"
        f"regions emitted : {len(final)}\n"
        f"total bases     : {total:,} ({100 * total / genome:.3f}% of hg38 primary)\n"
    )
    if missing:
        sys.stderr.write(f"unmatched       : {','.join(missing)}\n")


if __name__ == "__main__":
    main()
