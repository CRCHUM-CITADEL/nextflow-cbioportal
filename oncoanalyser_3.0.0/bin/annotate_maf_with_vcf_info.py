#!/usr/bin/env python3
"""Carry SAGE/PAVE INFO annotations from the source VCF into a mafsmith MAF.

mafsmith's `--retain-ann` only reads CSQ subfields, so every plain INFO field the
caller wrote -- including PAVE's gnomAD frequency (GND_FREQ) and SAGE's TIER -- is
dropped on the floor. This joins them back on the MAF's own coordinate convention.

The join key is (Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2)
with the VCF side normalized exactly the way mafsmith does it
(src/vcf/normalization.rs): strip the common REF/ALT prefix advancing POS, render an
emptied allele as "-", then Start = pos-1 for an insertion and pos otherwise. Keying
on raw VCF POS/REF/ALT instead would silently miss every indel.

Columns are always emitted, even when the VCF declares none of the fields, so every
per-subject MAF has identical width -- the group-level `collectFile` merge would
otherwise produce a ragged file.
"""

import argparse
import gzip
import sys

# INFO id -> output MAF column. gnomAD_AF is spelled the conventional way; the rest
# are namespaced so cBioPortal can surface them via `namespaces:` in the meta file.
INFO_FIELDS = [
    ("GND_FREQ", "gnomAD_AF"),
    ("TIER", "TIER"),
    ("CLNSIG", "CLNSIG"),
    ("CLNSIGCONF", "CLNSIGCONF"),
    ("PON_COUNT", "PON_COUNT"),
    ("MAPPABILITY", "MAPPABILITY"),
    ("MSG", "MSG"),
    ("TNC", "TNC"),
    ("REP_C", "REP_C"),
    ("MH", "MH"),
]


def normalize(pos, ref, alt):
    """Port of mafsmith's normalize() + maf_positions()."""
    prefix = 0
    for r, a in zip(ref, alt):
        if r != a:
            break
        prefix += 1
    if prefix == len(ref) and prefix == len(alt):
        prefix = max(0, prefix - 1)
    pos += prefix
    ref_n = ref[prefix:] or "-"
    alt_n = alt[prefix:] or "-"
    start = pos - 1 if ref_n == "-" else pos
    return start, ref_n, alt_n


def key(chrom, start, ref, alt):
    return (chrom[3:] if chrom.startswith("chr") else chrom, str(start), ref, alt)


def parse_info(field):
    out = {}
    for item in field.split(";"):
        if not item:
            continue
        k, _, v = item.partition("=")
        out[k] = v if _ else "1"  # a bare Flag reads as present
    return out


def read_vcf(path):
    opener = gzip.open if path.endswith(".gz") else open
    index = {}
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            chrom, pos, ref, alts, info = f[0], int(f[1]), f[3], f[4], f[7]
            values = parse_info(info)
            wanted = {i: values.get(i, "") for i, _ in INFO_FIELDS}
            # Index every ALT: mafsmith picks one "effective" alt from a multi-allelic
            # record and we cannot know which, so make them all resolvable.
            for alt in alts.split(","):
                index[key(chrom, *normalize(pos, ref, alt))] = wanted
    return index


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--vcf", required=True)
    p.add_argument("--maf", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--namespace", default="HMF")
    p.add_argument(
        "--min-match-rate",
        type=float,
        default=0.0,
        help="Fail when fewer than this fraction of MAF rows find a VCF record",
    )
    args = p.parse_args()

    index = read_vcf(args.vcf)
    print(f"Indexed {len(index)} VCF allele(s) from {args.vcf}", file=sys.stderr)

    new_cols = [
        col if col == "gnomAD_AF" else f"{args.namespace}.{col}"
        for _, col in INFO_FIELDS
    ]
    blanks = [""] * len(new_cols)

    total = matched = 0
    with open(args.maf) as fh, open(args.output, "w") as out:
        version = fh.readline().rstrip("\n")
        header = fh.readline().rstrip("\n").split("\t")
        out.write(version + "\n")
        out.write("\t".join(header + new_cols) + "\n")

        idx = {c: i for i, c in enumerate(header)}
        for required in ("Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"):
            if required not in idx:
                sys.exit(f"ERROR: MAF is missing the {required} column")

        for line in fh:
            row = line.rstrip("\n").split("\t")
            total += 1
            k = key(
                row[idx["Chromosome"]],
                row[idx["Start_Position"]],
                row[idx["Reference_Allele"]],
                row[idx["Tumor_Seq_Allele2"]],
            )
            hit = index.get(k)
            if hit is None:
                out.write("\t".join(row + blanks) + "\n")
            else:
                matched += 1
                out.write("\t".join(row + [hit[i] for i, _ in INFO_FIELDS]) + "\n")

    rate = (matched / total) if total else 1.0
    print(f"Annotated {matched}/{total} MAF rows ({rate:.1%}) from VCF INFO", file=sys.stderr)
    if total and index and matched == 0:
        sys.exit("ERROR: no MAF row matched any VCF record -- the join key is wrong")
    if rate < args.min_match_rate:
        sys.exit(f"ERROR: match rate {rate:.1%} below --min-match-rate {args.min_match_rate:.1%}")


if __name__ == "__main__":
    main()
