#!/usr/bin/env python3
"""
Fill MAF annotation columns that mafsmith/fastVEP leave empty, from the source VCF.

fastVEP has no dbSNP/ClinVar/COSMIC/1000G sources, so mafsmith always writes dbSNP_RS,
Existing_variation, CLIN_SIG and AF empty. Ampliseq Pisces VCFs carry that information already:
rsIDs in the ID column, and Nirvana INFO fields clinvar / cosmic (GenotypeIndex|value) and
AF1000G (Number=A). Each MAF row is matched to its VCF record by chromosome (chr prefix
ignored), MAF-style start position and alt allele; only empty cells are filled.

Usage: maf_add_vcf_annotations.py <in.maf> <source.vcf[.gz]> <out.maf>
"""
import csv
import gzip
import sys

FILLABLE = ("dbSNP_RS", "Existing_variation", "CLIN_SIG", "AF")
EMPTY = ("", ".", "-", "NA")


def nochr(c):
    return c[3:] if c.startswith("chr") else c


def maf_key(chrom, pos, ref, alt):
    """(chrom, Start_Position, Tumor_Seq_Allele2) as vcf2maf/mafsmith write them: the shared
    leading base(s) are trimmed; an insertion starts at the base before the inserted sequence."""
    p = 0
    while p < min(len(ref), len(alt)) and ref[p] == alt[p] and (len(ref) > 1 or len(alt) > 1):
        p += 1
    ref_t, alt_t = ref[p:] or "-", alt[p:] or "-"
    start = pos + p - 1 if ref_t == "-" else pos + p
    return nochr(chrom), str(start), alt_t.upper()


def per_allele(value, allele_index):
    """Values of a Nirvana 'GenotypeIndex|value,...' field that belong to allele_index (1-based)."""
    out = []
    for item in value.split(","):
        idx, _, val = item.partition("|")
        if val and idx == str(allele_index) and val not in out:
            out.append(val)
    return out


def load_vcf(path):
    opener = gzip.open if path.endswith(".gz") else open
    ann = {}
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            chrom, pos, vid, ref, alts, info = f[0], int(f[1]), f[2], f[3].upper(), f[4], f[7]
            kv = dict(x.split("=", 1) if "=" in x else (x, "") for x in info.split(";") if x)
            ids = [i for i in vid.split(";") if i not in EMPTY]
            for n, alt in enumerate(alts.split(","), start=1):
                if alt in (".", "*") or alt.startswith("<") or "[" in alt or "]" in alt:
                    continue
                af = kv.get("AF1000G", "").split(",")
                cosmic = per_allele(kv.get("cosmic", ""), n)
                ann[maf_key(chrom, pos, ref, alt.upper())] = {
                    "dbSNP_RS": ",".join(i for i in ids if i.startswith("rs")),
                    "Existing_variation": ",".join(ids + [c for c in cosmic if c not in ids]),
                    "CLIN_SIG": ",".join(per_allele(kv.get("clinvar", ""), n)),
                    "AF": af[n - 1] if len(af) >= n and af[n - 1] not in EMPTY else "",
                }
    return ann


def main():
    if len(sys.argv) != 4:
        print(__doc__, file=sys.stderr)
        sys.exit(1)
    maf_in, vcf, maf_out = sys.argv[1:]
    ann = load_vcf(vcf)

    with open(maf_in, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        rows = list(reader)
    col = {name: i for i, name in enumerate(header)}
    missing = [c for c in ("Chromosome", "Start_Position", "Tumor_Seq_Allele2") if c not in col]
    if missing:
        sys.exit(f"ERROR: MAF lacks {', '.join(missing)}")
    for name in FILLABLE:            # add any fillable column mafsmith did not emit
        if name not in col:
            col[name] = len(header)
            header.append(name)

    matched = filled = 0
    for row in rows:
        row.extend([""] * (len(header) - len(row)))
        key = (nochr(row[col["Chromosome"]]), row[col["Start_Position"]], row[col["Tumor_Seq_Allele2"]].upper())
        hit = ann.get(key)
        if not hit:
            continue
        matched += 1
        for name in FILLABLE:
            if row[col[name]] in EMPTY and hit[name]:
                row[col[name]] = hit[name]
                filled += 1

    with open(maf_out, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
    print(f"maf_add_vcf_annotations: {matched}/{len(rows)} MAF rows matched to VCF records, "
          f"{filled} cells filled", file=sys.stderr)


if __name__ == "__main__":
    main()
