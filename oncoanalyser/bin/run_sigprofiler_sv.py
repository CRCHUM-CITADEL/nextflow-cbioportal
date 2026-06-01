#!/usr/bin/env python3
"""
Extract structural variant (SV) type counts from an ESVEE VCF, fit against
COSMIC v3.6 SV-32 reference using scipy NNLS, and produce cBioPortal
GENERIC_ASSAY contribution and counts files.

Inputs:
  --vcf           : ESVEE somatic VCF (.vcf.gz or .vcf) with BND records
  --signatures_db : COSMIC_Human_SV-32_GRCh38_v3.6.csv (extracted from zip)
  --metadata      : cosmic_sv_metadata.tsv (sv_id, category, etiology, main_effect)
  --sample        : sample ID for output column header
  --output_contrib: output path for contribution file
  --output_counts : output path for counts file
"""

import argparse
import csv
import gzip
import os
import re
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.optimize import nnls


# Clustering threshold: breakpoints within this distance (bp) are considered clustered
CLUSTERING_THRESHOLD_BP = 25000

# Size bins for intra-chromosomal SVs (label, min_bp, max_bp)
SIZE_BINS = [
    ("1-10Kb", 1_000, 10_000),
    ("10-100Kb", 10_000, 100_000),
    ("100Kb-1Mb", 100_000, 1_000_000),
    ("1Mb-10Mb", 1_000_000, 10_000_000),
    (">10Mb", 10_000_000, float("inf")),
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--vcf", required=True)
    p.add_argument("--signatures_db", required=True)
    p.add_argument("--metadata", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--output_contrib", required=True)
    p.add_argument("--output_counts", required=True)
    return p.parse_args()


def open_vcf(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def get_bnd_strand(alt):
    """Determine local strand from BND ALT field.

    REF before bracket (e.g., N[chr:pos[, N]chr:pos]) → +
    REF after bracket (e.g., [chr:pos[N, ]chr:pos]N) → -
    """
    if re.match(r"^[A-Za-z.]+[\[\]]", alt):
        return "+"
    return "-"


def get_mate_info(alt):
    """Parse BND ALT to extract mate chrom and position."""
    match = re.search(r"[\[\]](.*?):(\d+)[\[\]]", alt)
    if not match:
        return None, None
    return match.group(1), int(match.group(2))


def classify_sv_type(chrom1, chrom2, strand1, strand2):
    """Classify SV type from strand orientation.

    Same-chromosome: (+,-) DEL, (-,+) TDS, (+,+)/(-,-) INV
    Different chromosome: TRANS
    """
    if chrom1 != chrom2:
        return "trans"
    if strand1 == "+" and strand2 == "-":
        return "del"
    if strand1 == "-" and strand2 == "+":
        return "tds"
    return "inv"


def get_size_bin(size_bp):
    """Return the COSMIC size bin label for a given SV size."""
    for label, min_bp, max_bp in SIZE_BINS:
        if min_bp <= size_bp < max_bp:
            return label
    return ">10Mb"


def parse_bnd_records(vcf_path):
    """Parse ESVEE VCF and return list of SV events as dicts.

    Each SV event: {chrom1, pos1, strand1, chrom2, pos2, strand2, sv_type, size}
    """
    bnds = {}  # id -> (chrom, pos, alt, mate_id)

    with open_vcf(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 8:
                continue

            filt = fields[6]
            if filt != "PASS" and filt != ".":
                continue

            info = fields[7]

            # Check SVTYPE=BND
            if "SVTYPE=BND" not in info:
                continue

            record_id = fields[2]
            chrom = fields[0]
            pos = int(fields[1])
            alt = fields[4]

            # Extract MATEID
            mate_match = re.search(r"MATEID=([^;]+)", info)
            if not mate_match:
                continue
            mate_id = mate_match.group(1)

            bnds[record_id] = (chrom, pos, alt, mate_id)

    # Pair BND records
    seen = set()
    sv_events = []

    for record_id, (chrom1, pos1, alt1, mate_id) in bnds.items():
        pair_key = tuple(sorted([record_id, mate_id]))
        if pair_key in seen:
            continue
        seen.add(pair_key)

        if mate_id not in bnds:
            continue

        chrom2, pos2, alt2, _ = bnds[mate_id]
        strand1 = get_bnd_strand(alt1)
        strand2 = get_bnd_strand(alt2)
        sv_type = classify_sv_type(chrom1, chrom2, strand1, strand2)

        size = abs(pos2 - pos1) if chrom1 == chrom2 else 0

        sv_events.append(
            {
                "chrom1": chrom1,
                "pos1": pos1,
                "chrom2": chrom2,
                "pos2": pos2,
                "sv_type": sv_type,
                "size": size,
            }
        )

    return sv_events


def detect_clustering(sv_events):
    """Detect SV clustering based on breakpoint proximity.

    An SV is 'clustered' if any of its breakpoints is within
    CLUSTERING_THRESHOLD_BP of another SV's breakpoint.
    """
    if len(sv_events) <= 1:
        return [False] * len(sv_events)

    # Collect all breakpoints with their SV index
    breakpoints = []
    for i, sv in enumerate(sv_events):
        breakpoints.append((sv["chrom1"], sv["pos1"], i))
        breakpoints.append((sv["chrom2"], sv["pos2"], i))

    # Sort by chrom and position
    breakpoints.sort(key=lambda x: (x[0], x[1]))

    # Find clustered SVs
    clustered_svs = set()
    for j in range(len(breakpoints) - 1):
        chrom_a, pos_a, idx_a = breakpoints[j]
        chrom_b, pos_b, idx_b = breakpoints[j + 1]
        if chrom_a == chrom_b and idx_a != idx_b:
            if abs(pos_b - pos_a) <= CLUSTERING_THRESHOLD_BP:
                clustered_svs.add(idx_a)
                clustered_svs.add(idx_b)

    return [i in clustered_svs for i in range(len(sv_events))]


def classify_sv_32(sv_events):
    """Classify each SV into one of 32 COSMIC SV types.

    Returns a list of COSMIC SV-32 type labels.
    """
    is_clustered = detect_clustering(sv_events)
    labels = []

    for i, sv in enumerate(sv_events):
        cluster_prefix = "clustered" if is_clustered[i] else "non-clustered"
        sv_type = sv["sv_type"]

        if sv_type == "trans":
            labels.append(f"{cluster_prefix}_trans")
        else:
            size_bin = get_size_bin(sv["size"])
            labels.append(f"{cluster_prefix}_{sv_type}_{size_bin}")

    return labels


def build_sv_matrix(sv_labels, cosmic_types, sample_id):
    """Build SV-32 count matrix."""
    counts = {t: 0 for t in cosmic_types}
    for label in sv_labels:
        if label in counts:
            counts[label] += 1
        else:
            print(f"WARNING: SV type '{label}' not in COSMIC-32, skipping", file=sys.stderr)

    df = pd.DataFrame(
        {"Type": list(counts.keys()), sample_id: list(counts.values())}
    )
    df.set_index("Type", inplace=True)
    return df


def fit_nnls(matrix_df, sig_ref_df, sample_id):
    """Fit SV signatures using non-negative least squares.

    Returns dict of {signature_id: proportion} for non-zero exposures.
    """
    # Align matrix rows with reference rows
    common_types = matrix_df.index.intersection(sig_ref_df.index)
    if len(common_types) == 0:
        return {}

    M = matrix_df.loc[common_types, sample_id].values.astype(float)
    S = sig_ref_df.loc[common_types].values.astype(float)

    if M.sum() == 0:
        return {}

    exposures, _ = nnls(S, M)
    total = exposures.sum()
    if total == 0:
        return {}

    proportions = exposures / total
    result = {}
    for sig_id, prop in zip(sig_ref_df.columns, proportions):
        if prop > 0:
            result[sig_id] = round(float(prop), 5)

    return result


def load_metadata(metadata_path, id_col):
    meta = {}
    with open(metadata_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            meta[row[id_col]] = row
    return meta


def write_contribution(proportions, metadata, sample_id, output_path):
    """Write SV signature contribution in cBioPortal format."""
    if not proportions:
        with open(output_path, "w") as fh:
            fh.write(f"ENTITY_STABLE_ID\tNAME\tDESCRIPTION\t{sample_id}\n")
        return

    rows = []
    for sig_id, prop in proportions.items():
        sig_id = sig_id.strip()
        meta = metadata.get(sig_id, {})
        category = meta.get("category", "Unknown")
        etiology = meta.get("etiology", f"Mutational signature {sig_id}")
        main_effect = meta.get("main_effect", "")

        rows.append(
            {
                "ENTITY_STABLE_ID": f"mutational_signatures_contribution_{sig_id}",
                "NAME": f"{sig_id} ({category})" if category else sig_id,
                "DESCRIPTION": f"{etiology} ({main_effect})" if main_effect else etiology,
                sample_id: prop,
            }
        )

    out_df = pd.DataFrame(rows).sort_values("ENTITY_STABLE_ID")
    out_df.to_csv(output_path, sep="\t", index=False)
    print(f"Written {len(rows)} SV contribution rows to {output_path}")


def write_counts(matrix_df, sample_id, output_path):
    """Write SV counts in cBioPortal GENERIC_ASSAY format."""
    rows = []
    for sv_type, count in matrix_df[sample_id].items():
        # Extract base type for CATEGORY colour grouping
        parts = sv_type.split("_")
        # e.g., "clustered_del_1-10Kb" → base_type = "del"
        # e.g., "non-clustered_trans" → base_type = "trans"
        if parts[0] == "non-clustered":
            base_type = parts[1]
        else:
            base_type = parts[1]

        # Sanitize ">" to "gt" for cBioPortal (only alphanumeric, _, - allowed)
        safe_type = sv_type.replace(">", "gt")
        rows.append(
            {
                "ENTITY_STABLE_ID": f"mutational_signatures_matrix_SV_{safe_type}",
                "NAME": safe_type,
                "CATEGORY": base_type.upper(),
                sample_id: int(count),
            }
        )

    out_df = pd.DataFrame(rows).sort_values("ENTITY_STABLE_ID")
    out_df.to_csv(output_path, sep="\t", index=False)
    print(f"Written {len(rows)} SV counts rows to {output_path}")


def main():
    args = parse_args()

    # Load COSMIC SV-32 reference
    print(f"Loading SV signatures DB: {args.signatures_db}")
    sig_ref = pd.read_csv(args.signatures_db, index_col=0)
    sig_ref.columns = [c.strip() for c in sig_ref.columns]
    cosmic_types = list(sig_ref.index)

    # Parse ESVEE VCF
    print(f"Parsing SVs from VCF: {args.vcf}")
    sv_events = parse_bnd_records(args.vcf)
    print(f"Found {len(sv_events)} SV events")

    # Classify into COSMIC SV-32 types
    sv_labels = classify_sv_32(sv_events) if sv_events else []

    # Build count matrix
    matrix_df = build_sv_matrix(sv_labels, cosmic_types, args.sample)
    total_svs = int(matrix_df[args.sample].sum())
    print(f"Total classified SVs: {total_svs}")

    # Write counts
    write_counts(matrix_df, args.sample, args.output_counts)

    # Fit using NNLS
    if total_svs == 0:
        print("WARNING: No SVs found, writing empty contribution file", file=sys.stderr)
        with open(args.output_contrib, "w") as fh:
            fh.write(f"ENTITY_STABLE_ID\tNAME\tDESCRIPTION\t{args.sample}\n")
        return

    print("Fitting SV signatures using NNLS...")
    proportions = fit_nnls(matrix_df, sig_ref, args.sample)

    # Write contribution
    metadata = load_metadata(args.metadata, "sv_id")
    write_contribution(proportions, metadata, args.sample, args.output_contrib)


if __name__ == "__main__":
    main()
