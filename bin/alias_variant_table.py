#!/usr/bin/env python3
"""
alias_variant_table.py

Write a copy of a clinical variant table carrying additional column names, so
that a variant browser expecting a different capitalisation or naming
convention can read it.

Why this is needed
------------------
The dashboard's variant browser resolves columns by exact name. This
pipeline's filter emits mostly lowercase names (gene, consequence, impact,
tumor_af_pct) while the browser was written against a table using capitalised
ones. The three columns that happen to match in both (REF_COUNT, ALT_COUNT,
Filter) are precisely the three that render; everything else comes back empty,
which looks like a missing file but is a naming mismatch on a file that loaded
fine.

Renaming the columns outright would break this pipeline's own parsers, which
read the lowercase names. So the alias copy carries both: every original
column is preserved, and an aliased duplicate is appended for each name the
browser might be looking for. Extra columns are inert to a name-based reader.

No value is transformed, invented or reordered. An alias column is a verbatim
copy of its source column, and an alias whose source is absent is not written.

Usage:
    python3 bin/alias_variant_table.py <source.tsv> <destination.tsv>
    python3 bin/alias_variant_table.py <source.tsv> <destination.tsv> --report

Standard library only.
"""

import argparse
import csv
import os
import sys


# alias name -> source column in this pipeline's filter output.
# Several aliases may point at the same source; that is intentional, since the
# exact convention the browser uses is not knowable from the outside and a
# superset costs nothing.
ALIASES = [
    # Locus
    ("CHROM", "chrom"),
    ("Chrom", "chrom"),
    ("Chr", "chrom"),
    ("POS", "pos"),
    ("Pos", "pos"),
    ("Position", "pos"),
    ("POSITION", "pos"),
    ("position", "pos"),
    ("pos_hg38", "pos"),
    ("POS_HG38", "pos"),
    ("hg38_pos", "pos"),
    ("g_pos", "pos"),
    ("genomic_pos", "pos"),
    ("Start_Position", "pos"),
    # variant-browser.js reads r.Start when it builds the GeneBe link
    # (genebeUrl(chrClean, r.Start, r.Ref, r.Alt)). The position is 1-based
    # VCF coordinate in both, so this is a straight copy, not a conversion.
    ("Start", "pos"),
    ("REF", "ref"),
    ("Ref", "ref"),
    ("ALT", "alt"),
    ("Alt", "alt"),
    # Gene and annotation
    ("GENE", "gene"),
    ("Gene", "gene"),
    ("SYMBOL", "gene"),
    ("Symbol", "gene"),
    ("CONSEQUENCE", "consequence"),
    ("Consequence", "consequence"),
    ("IMPACT", "impact"),
    ("Impact", "impact"),
    ("TRANSCRIPT", "transcript"),
    ("Transcript", "transcript"),
    ("Feature", "transcript"),
    ("HGVSc", "hgvsc"),
    ("HGVSC", "hgvsc"),
    ("HGVSp", "hgvsp"),
    ("HGVSP", "hgvsp"),
    ("EXON", "exon_rank"),
    ("Exon", "exon_rank"),
    ("BIOTYPE", "biotype"),
    ("CANONICAL", "canonical"),
    # Identifiers and population data
    ("rsID", "rs_id"),
    ("RSID", "rs_id"),
    ("dbSNP", "rs_id"),
    ("Existing_variation", "rs_id"),
    ("gnomAD_AF", "pop_af_max"),
    ("POP_AF", "pop_af_max"),
    ("CLINVAR", "clinvar_sig"),
    ("ClinVar", "clinvar_sig"),
    ("CLNSIG", "clinvar_sig"),
    # Allele fraction and depth
    ("VAF", "tumor_af"),
    ("vaf", "tumor_af"),
    ("AF", "tumor_af"),
    ("VAF%", "tumor_af_pct"),
    ("VAF_pct", "tumor_af_pct"),
    ("VAF_PCT", "tumor_af_pct"),
    ("VAFpct", "tumor_af_pct"),
    ("VAF_percent", "tumor_af_pct"),
    ("DEPTH", "DP"),
    ("Depth", "DP"),
    ("TotalDepth", "DP"),
    ("RefCount", "REF_COUNT"),
    ("AltCount", "ALT_COUNT"),
    ("REF_DP", "REF_COUNT"),
    ("ALT_DP", "ALT_COUNT"),
    # Filter status
    ("FILTER", "Filter"),
    ("filter", "Filter"),
    ("FilterStatus", "Filter"),
    # Variant class
    ("TYPE", "variant_type"),
    ("Type", "variant_type"),
    ("VariantType", "variant_type"),
    ("QUAL", "qual"),
    ("Qual", "qual"),
]


# Columns the browser reads that have no single source column, so they are
# composed rather than copied. Each entry is (name, builder, description);
# the builder receives a dict of the source row and returns a string, or ""
# when its inputs are absent. A derived column is never written if it would
# be empty for every row.
def _exon(row):
    """VEP-style exon designation, "rank/total"."""
    rank = (row.get("exon_rank") or "").strip()
    total = (row.get("exon_total") or "").strip()
    if rank and total:
        return "%s/%s" % (rank, total)
    return rank


DERIVED = [
    ("EXON", _exon, "exon_rank/exon_total"),
]


def read_tsv(path):
    with open(path, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = [h.strip() for h in next(reader)]
        except StopIteration:
            return [], []
        rows = [
            list(row) + [""] * (len(header) - len(row))
            for row in reader
            if row and any(cell.strip() for cell in row)
        ]
    return header, rows


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("source", help="clinical or filtered TSV from the filter")
    parser.add_argument("destination", help="aliased copy to write")
    parser.add_argument(
        "--caller-count",
        default=None,
        metavar="N",
        help="write VariantCaller_Count with this constant value. The browser "
        "shows a caller count and offers >2/>3/>4 filters; this pipeline calls "
        "somatic SNVs with a single caller, so the column is absent by default "
        "rather than asserted. Pass 1 to populate it.",
    )
    parser.add_argument(
        "--verdict-from",
        default=None,
        metavar="COLUMN",
        help="populate SomaticSeq_Verdict from this column, e.g. --verdict-from "
        "Filter. Off by default: the label names a caller this pipeline does "
        "not run, and an empty filter is more honest than one that presents "
        "output under another tool's name.",
    )
    parser.add_argument(
        "--extra-alias",
        action="append",
        default=[],
        metavar="NAME=SOURCE",
        help="add an alias not in the built-in list, e.g. --extra-alias "
        "start=pos. Repeatable. Use this when a consumer reads a column name "
        "the built-in set does not cover.",
    )
    parser.add_argument(
        "--report",
        action="store_true",
        help="list the alias columns that were added",
    )
    args = parser.parse_args(argv)

    if not os.path.isfile(args.source):
        print("ERROR: source not found: %s" % args.source, file=sys.stderr)
        return 2

    header, rows = read_tsv(args.source)
    if not header:
        # An empty input is a legitimate outcome; copy the emptiness through
        # rather than failing the bundle.
        with open(args.destination, "w") as handle:
            handle.write("")
        if args.report:
            print("source is empty; wrote empty destination")
        return 0

    index = {name: i for i, name in enumerate(header)}

    aliases = list(ALIASES)
    for spec in args.extra_alias:
        if "=" not in spec:
            print("ERROR: --extra-alias needs NAME=SOURCE, got %r" % spec,
                  file=sys.stderr)
            return 2
        name, source = spec.split("=", 1)
        aliases.append((name.strip(), source.strip()))

    added = []
    for alias, source in aliases:
        if alias in index:
            continue          # already present under that name
        if source not in index:
            continue          # nothing to copy from
        added.append((alias, source, index[source]))

    # Derived columns are composed from more than one source, so they are
    # evaluated per row rather than copied by index.
    dict_rows = [dict(zip(header, row)) for row in rows]

    derived = []
    for name, builder, description in DERIVED:
        if name in index:
            continue
        values = [builder(d) for d in dict_rows]
        if any(v for v in values):
            derived.append((name, description, values))

    constants = []
    if args.caller_count is not None:
        if "VariantCaller_Count" not in index:
            constants.append(("VariantCaller_Count", str(args.caller_count)))
    if args.verdict_from is not None:
        if args.verdict_from not in index:
            print("ERROR: --verdict-from column not present: %s"
                  % args.verdict_from, file=sys.stderr)
            return 2
        if "SomaticSeq_Verdict" not in index:
            derived.append((
                "SomaticSeq_Verdict",
                "copy of %s" % args.verdict_from,
                [d.get(args.verdict_from, "") for d in dict_rows],
            ))

    out_header = (
        list(header)
        + [alias for alias, _s, _i in added]
        + [name for name, _d, _v in derived]
        + [name for name, _v in constants]
    )

    os.makedirs(os.path.dirname(os.path.abspath(args.destination)) or ".",
                exist_ok=True)
    with open(args.destination, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(out_header)
        for position, row in enumerate(rows):
            writer.writerow(
                row
                + [row[i] for _a, _s, i in added]
                + [values[position] for _n, _d, values in derived]
                + [value for _n, value in constants]
            )

    if args.report:
        print("%s -> %s" % (args.source, args.destination))
        print("  rows            : %d" % len(rows))
        print("  source columns  : %d" % len(header))
        print("  alias columns   : %d" % len(added))
        for alias, source, _i in added:
            print("    %-22s <- %s" % (alias, source))
        if derived:
            print("  derived columns : %d" % len(derived))
            for name, description, _v in derived:
                print("    %-22s <- %s" % (name, description))
        if constants:
            print("  constant columns: %d" % len(constants))
            for name, value in constants:
                print("    %-22s =  %s" % (name, value))
    return 0


if __name__ == "__main__":
    sys.exit(main())
