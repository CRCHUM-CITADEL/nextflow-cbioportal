#!/usr/bin/env python3
"""
Extract indels from a somatic VCF, classify into the COSMIC ID-83 scheme
using reference FASTA context (via pysam), fit against COSMIC v3.6 ID-83
reference using SigProfilerAssignment, and produce cBioPortal GENERIC_ASSAY
contribution and counts files.

Inputs:
  --vcf           : PAVE somatic VCF (.vcf.gz or .vcf)
  --fasta         : Reference genome FASTA (indexed with .fai)
  --signatures_db : COSMIC_Human_ID-83_GRCh38_v3.6.csv
  --metadata      : cosmic_id_metadata.tsv (id_id, category, etiology, main_effect)
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
import tempfile

import pandas as pd
import pysam

# The 83 COSMIC ID mutation types in canonical order
ID83_TYPES = [
    "1:Del:C:0", "1:Del:C:1", "1:Del:C:2", "1:Del:C:3", "1:Del:C:4", "1:Del:C:5",
    "1:Del:T:0", "1:Del:T:1", "1:Del:T:2", "1:Del:T:3", "1:Del:T:4", "1:Del:T:5",
    "1:Ins:C:0", "1:Ins:C:1", "1:Ins:C:2", "1:Ins:C:3", "1:Ins:C:4", "1:Ins:C:5",
    "1:Ins:T:0", "1:Ins:T:1", "1:Ins:T:2", "1:Ins:T:3", "1:Ins:T:4", "1:Ins:T:5",
    "2:Del:R:0", "2:Del:R:1", "2:Del:R:2", "2:Del:R:3", "2:Del:R:4", "2:Del:R:5",
    "3:Del:R:0", "3:Del:R:1", "3:Del:R:2", "3:Del:R:3", "3:Del:R:4", "3:Del:R:5",
    "4:Del:R:0", "4:Del:R:1", "4:Del:R:2", "4:Del:R:3", "4:Del:R:4", "4:Del:R:5",
    "5:Del:R:0", "5:Del:R:1", "5:Del:R:2", "5:Del:R:3", "5:Del:R:4", "5:Del:R:5",
    "2:Ins:R:0", "2:Ins:R:1", "2:Ins:R:2", "2:Ins:R:3", "2:Ins:R:4", "2:Ins:R:5",
    "3:Ins:R:0", "3:Ins:R:1", "3:Ins:R:2", "3:Ins:R:3", "3:Ins:R:4", "3:Ins:R:5",
    "4:Ins:R:0", "4:Ins:R:1", "4:Ins:R:2", "4:Ins:R:3", "4:Ins:R:4", "4:Ins:R:5",
    "5:Ins:R:0", "5:Ins:R:1", "5:Ins:R:2", "5:Ins:R:3", "5:Ins:R:4", "5:Ins:R:5",
    "2:Del:M:1",
    "3:Del:M:1", "3:Del:M:2",
    "4:Del:M:1", "4:Del:M:2", "4:Del:M:3",
    "5:Del:M:1", "5:Del:M:2", "5:Del:M:3", "5:Del:M:4", "5:Del:M:5",
]

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--vcf", required=True)
    p.add_argument("--fasta", required=True)
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


def extract_indels_from_vcf(vcf_path):
    """Extract indels from a somatic VCF.

    Returns list of tuples: (chrom, pos, ref, alt) where len(ref) != len(alt).
    Position is 1-based as in VCF.
    """
    indels = []

    with open_vcf(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 7:
                continue

            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3].upper()
            alt = fields[4].upper().split(",")[0]
            filt = fields[6]

            if filt != "PASS" and filt != ".":
                continue

            if not re.match(r"^[ACGT]+$", ref) or not re.match(r"^[ACGT]+$", alt):
                continue

            # Only indels: ref and alt differ in length
            if len(ref) == len(alt):
                continue

            indels.append((chrom, pos, ref, alt))

    return indels


def classify_1bp_del(deleted_base, chrom, pos, ref, fasta):
    """Classify a 1bp deletion into the ID-83 scheme.

    Args:
        deleted_base: the single deleted nucleotide (A/C/G/T)
        chrom: chromosome
        pos: 1-based position of the VCF record
        ref: full REF allele
        fasta: pysam.FastaFile

    Returns: ID-83 type string, e.g. '1:Del:C:3'
    """
    # Normalize to C or T
    if deleted_base in ("G", "A"):
        base = reverse_complement(deleted_base)
    else:
        base = deleted_base

    # The deleted base is at position pos+1 (0-based: pos) in the reference
    # (VCF convention: REF = anchor + deleted bases)
    del_pos_0based = pos  # 0-based position of deleted base

    # Count homopolymer run of the same base at and after the deletion site
    # Fetch context: up to 50bp downstream from the deletion position
    try:
        context = fasta.fetch(chrom, del_pos_0based, del_pos_0based + 50).upper()
    except (ValueError, KeyError):
        context = ""

    # Count how many consecutive bases match the deleted base at the site
    homopolymer_len = 0
    for c in context:
        if c == deleted_base:
            homopolymer_len += 1
        else:
            break

    # Cap at 5
    homopolymer_len = min(homopolymer_len, 5)

    return f"1:Del:{base}:{homopolymer_len}"


def classify_1bp_ins(inserted_base, chrom, pos, fasta):
    """Classify a 1bp insertion into the ID-83 scheme.

    Args:
        inserted_base: the single inserted nucleotide
        chrom: chromosome
        pos: 1-based VCF position (anchor base position)
        fasta: pysam.FastaFile

    Returns: ID-83 type string, e.g. '1:Ins:T:2'
    """
    # Normalize to C or T
    if inserted_base in ("G", "A"):
        base = reverse_complement(inserted_base)
    else:
        base = inserted_base

    # Count homopolymer run of the same base starting after the anchor
    ins_pos_0based = pos  # 0-based position right after anchor

    try:
        context = fasta.fetch(chrom, ins_pos_0based, ins_pos_0based + 50).upper()
    except (ValueError, KeyError):
        context = ""

    homopolymer_len = 0
    for c in context:
        if c == inserted_base:
            homopolymer_len += 1
        else:
            break

    homopolymer_len = min(homopolymer_len, 5)

    return f"1:Ins:{base}:{homopolymer_len}"


def count_repeat_units(indel_seq, chrom, pos_after_anchor, indel_len, fasta, is_del):
    """Count how many times the indel sequence repeats at the genomic position.

    For deletions: count repeats of the deleted sequence starting right after
    the deletion site.
    For insertions: count repeats of the inserted sequence starting at the
    insertion site.

    Returns: number of additional repeat units (0 if no repeat context).
    """
    unit = indel_seq.upper()
    unit_len = len(unit)

    # Fetch downstream context (enough for up to 6 repeat units)
    fetch_len = unit_len * 7
    try:
        if is_del:
            # For deletions, look at sequence after the deleted region
            start = pos_after_anchor + indel_len
        else:
            # For insertions, look at sequence starting at insertion point
            start = pos_after_anchor
        context = fasta.fetch(chrom, start, start + fetch_len).upper()
    except (ValueError, KeyError):
        return 0

    repeat_count = 0
    for i in range(0, len(context) - unit_len + 1, unit_len):
        if context[i : i + unit_len] == unit:
            repeat_count += 1
        else:
            break

    return repeat_count


def find_microhomology(chrom, del_start, del_end, fasta):
    """Find microhomology length at a deletion breakpoint.

    Microhomology: bases at the 3' end of the deleted sequence that match
    bases immediately after the deletion.

    Args:
        chrom: chromosome
        del_start: 0-based start of deleted region
        del_end: 0-based end of deleted region (exclusive)
        fasta: pysam.FastaFile

    Returns: microhomology length (0 if none)
    """
    del_len = del_end - del_start
    try:
        deleted_seq = fasta.fetch(chrom, del_start, del_end).upper()
        downstream = fasta.fetch(chrom, del_end, del_end + del_len).upper()
    except (ValueError, KeyError):
        return 0

    mh_len = 0
    for i in range(len(deleted_seq)):
        if i < len(downstream) and deleted_seq[i] == downstream[i]:
            mh_len += 1
        else:
            break

    return mh_len


def classify_indel(chrom, pos, ref, alt, fasta):
    """Classify a single indel into the COSMIC ID-83 scheme.

    Args:
        chrom, pos, ref, alt: VCF fields (pos is 1-based)
        fasta: pysam.FastaFile

    Returns: ID-83 type string
    """
    ref_len = len(ref)
    alt_len = len(alt)

    if ref_len > alt_len:
        # Deletion
        # Strip common prefix (anchor base)
        prefix_len = 0
        for i in range(min(ref_len, alt_len)):
            if ref[i] == alt[i]:
                prefix_len += 1
            else:
                break

        deleted_seq = ref[prefix_len:]
        indel_len = len(deleted_seq)
        # 0-based position of first deleted base
        del_start_0based = pos - 1 + prefix_len

        if indel_len == 1:
            return classify_1bp_del(deleted_seq, chrom, del_start_0based, ref, fasta)

        # Clamp indel size to 5 for classification
        size = min(indel_len, 5)

        # Check for repeat context
        repeat_units = count_repeat_units(
            deleted_seq, chrom, del_start_0based, indel_len, fasta, is_del=True
        )

        if repeat_units > 0:
            return f"{size}:Del:R:{min(repeat_units, 5)}"

        # Check for microhomology
        del_end_0based = del_start_0based + indel_len
        mh_len = find_microhomology(chrom, del_start_0based, del_end_0based, fasta)

        if mh_len > 0:
            # Cap microhomology at indel_len
            mh_len = min(mh_len, indel_len)
            # For the M category, max is size (capped at 5)
            mh_len = min(mh_len, size)
            return f"{size}:Del:M:{mh_len}"

        # No repeat, no microhomology: classify as R:0
        return f"{size}:Del:R:0"

    else:
        # Insertion
        prefix_len = 0
        for i in range(min(ref_len, alt_len)):
            if ref[i] == alt[i]:
                prefix_len += 1
            else:
                break

        inserted_seq = alt[prefix_len:]
        indel_len = len(inserted_seq)
        # 0-based position right after the anchor
        ins_pos_0based = pos - 1 + prefix_len

        if indel_len == 1:
            return classify_1bp_ins(inserted_seq, chrom, pos - 1 + prefix_len, fasta)

        # Clamp indel size to 5
        size = min(indel_len, 5)

        # Check for repeat context
        repeat_units = count_repeat_units(
            inserted_seq, chrom, ins_pos_0based, indel_len, fasta, is_del=False
        )

        if repeat_units > 0:
            return f"{size}:Ins:R:{min(repeat_units, 5)}"

        # No repeat: classify as R:0 (insertions don't have M category)
        return f"{size}:Ins:R:0"


def build_id_matrix(classifications, sample_id):
    """Build an ID-83 count matrix DataFrame."""
    counts = {t: 0 for t in ID83_TYPES}
    unmapped = 0
    for id_type in classifications:
        if id_type in counts:
            counts[id_type] += 1
        else:
            print(f"WARNING: ID type '{id_type}' not in ID-83, skipping", file=sys.stderr)
            unmapped += 1

    if unmapped > 0:
        print(f"WARNING: {unmapped} indels could not be mapped to ID-83 types", file=sys.stderr)

    df = pd.DataFrame(
        {"MutationType": list(counts.keys()), sample_id: list(counts.values())}
    )
    df.set_index("MutationType", inplace=True)
    return df


def prepare_signatures_db(cosmic_csv_path, tmp_dir):
    """Convert COSMIC ID-83 CSV to TSV for SigProfilerAssignment."""
    df = pd.read_csv(cosmic_csv_path, index_col=0)
    df.columns = [c.strip() for c in df.columns]
    # Rename index to match SigProfilerAssignment expectations
    df.index.name = "Mutation Types"
    out_path = os.path.join(tmp_dir, "COSMIC_ID83.tsv")
    df.to_csv(out_path, sep="\t", index=True)
    return out_path


def run_cosmic_fit_id(matrix_tsv_path, sig_db_tsv_path, output_dir):
    """Run SigProfilerAssignment for ID-83."""
    from SigProfilerAssignment import Analyzer

    Analyzer.cosmic_fit(
        samples=matrix_tsv_path,
        output=output_dir,
        signature_database=sig_db_tsv_path,
        genome_build="GRCh38",
        input_type="matrix",
        context_type="ID",
        make_plots=False,
        export_probabilities=False,
        collapse_to_SBS96=False,
        exome=False,
    )

    activities_path = os.path.join(
        output_dir,
        "Assignment_Solution",
        "Activities",
        "Assignment_Solution_Activities.txt",
    )
    if not os.path.exists(activities_path):
        sys.exit(f"ERROR: Expected output not found: {activities_path}")
    return activities_path


def load_metadata(metadata_path, id_col):
    """Load metadata TSV -> dict keyed by id column."""
    meta = {}
    with open(metadata_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            meta[row[id_col]] = row
    return meta


def write_contribution(activities_path, metadata, sample_id, output_path):
    """Convert Activities to cBioPortal contribution format."""
    acts = pd.read_csv(activities_path, sep="\t", index_col=0)

    if sample_id not in acts.index:
        sys.exit(f"ERROR: sample '{sample_id}' not in Activities: {list(acts.index)}")

    row = acts.loc[sample_id]
    total = row.sum()

    if total == 0:
        print(f"WARNING: Zero ID activity for {sample_id}", file=sys.stderr)
        with open(output_path, "w") as fh:
            fh.write(f"ENTITY_STABLE_ID\tNAME\tDESCRIPTION\t{sample_id}\n")
        return

    proportions = row / total
    proportions = proportions[proportions > 0]

    rows = []
    for sig_id, prop in proportions.items():
        sig_id = sig_id.strip()
        meta = metadata.get(sig_id, {})
        category = meta.get("category", "Unknown")
        etiology = meta.get("etiology", f"Mutational signature {sig_id}")
        main_effect = meta.get("main_effect", "")

        rows.append(
            {
                "ENTITY_STABLE_ID": f"mutational_signature_contribution_{sig_id}",
                "NAME": f"{sig_id} ({category})" if category else sig_id,
                "DESCRIPTION": f"{etiology} ({main_effect})" if main_effect else etiology,
                sample_id: round(prop, 5),
            }
        )

    out_df = pd.DataFrame(rows).sort_values("ENTITY_STABLE_ID")
    out_df.to_csv(output_path, sep="\t", index=False)
    print(f"Written {len(rows)} ID contribution rows to {output_path}")


def write_counts(matrix_df, sample_id, output_path):
    """Write ID counts in cBioPortal GENERIC_ASSAY format."""
    rows = []
    for id_type, count in matrix_df[sample_id].items():
        # Entity ID: replace colons with underscores
        entity_id = f"mutational_signatures_matrix_{id_type.replace(':', '_')}"
        rows.append(
            {
                "ENTITY_STABLE_ID": entity_id,
                "NAME": id_type,
                sample_id: int(count),
            }
        )

    out_df = pd.DataFrame(rows).sort_values("ENTITY_STABLE_ID")
    out_df.to_csv(output_path, sep="\t", index=False)
    print(f"Written {len(rows)} ID counts rows to {output_path}")


def main():
    args = parse_args()

    fasta = pysam.FastaFile(args.fasta)

    with tempfile.TemporaryDirectory() as tmp_dir:
        # Prepare COSMIC reference and get type list
        print(f"Preparing ID signatures DB: {args.signatures_db}")
        sig_db_tsv = prepare_signatures_db(args.signatures_db, tmp_dir)

        # Extract indels from VCF
        print(f"Extracting indels from VCF: {args.vcf}")
        indels = extract_indels_from_vcf(args.vcf)
        print(f"Found {len(indels)} indels")

        # Classify each indel into ID-83 scheme
        print("Classifying indels into ID-83 channels...")
        classifications = []
        for chrom, pos, ref, alt in indels:
            id_type = classify_indel(chrom, pos, ref, alt, fasta)
            classifications.append(id_type)

        # Build count matrix
        matrix_df = build_id_matrix(classifications, args.sample)
        total_id = int(matrix_df[args.sample].sum())
        print(f"Total classified indels: {total_id}")

        # Write counts file
        write_counts(matrix_df, args.sample, args.output_counts)

        if total_id == 0:
            print("WARNING: No indels found, writing empty contribution file", file=sys.stderr)
            with open(args.output_contrib, "w") as fh:
                fh.write(f"ENTITY_STABLE_ID\tNAME\tDESCRIPTION\t{args.sample}\n")
            fasta.close()
            return

        if total_id < 50:
            print(
                f"WARNING: Only {total_id} indels for {args.sample}. "
                "Results may be unreliable.",
                file=sys.stderr,
            )

        # Write matrix for SigProfilerAssignment
        matrix_tsv = os.path.join(tmp_dir, "id_matrix.tsv")
        matrix_df.to_csv(matrix_tsv, sep="\t")

        # Run SigProfilerAssignment
        spa_outdir = os.path.join(tmp_dir, "spa_output")
        os.makedirs(spa_outdir)
        print("Running SigProfilerAssignment for ID...")
        activities_path = run_cosmic_fit_id(matrix_tsv, sig_db_tsv, spa_outdir)

        # Write contribution file
        metadata = load_metadata(args.metadata, "id_id")
        write_contribution(activities_path, metadata, args.sample, args.output_contrib)

    fasta.close()


if __name__ == "__main__":
    main()
