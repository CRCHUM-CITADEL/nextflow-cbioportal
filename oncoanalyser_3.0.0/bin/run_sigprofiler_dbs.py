#!/usr/bin/env python3
"""
Extract doublet base substitutions (DBS) from a somatic VCF, fit against
COSMIC v3.6 DBS-78 reference using SigProfilerAssignment, and produce
cBioPortal GENERIC_ASSAY contribution and counts files.

Inputs:
  --vcf           : PAVE somatic VCF (.vcf.gz or .vcf)
  --signatures_db : COSMIC_Human_DBS-78_GRCh38_v3.6.csv (extracted from zip)
  --metadata      : cosmic_dbs_metadata.tsv (dbs_id, category, etiology, main_effect)
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

# The 10 canonical DBS reference dinucleotides used by COSMIC
# Dinucleotides not in this set must be reverse-complemented
CANONICAL_REFS = {"AC", "AT", "CC", "CG", "CT", "GC", "TA", "TC", "TG", "TT"}

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1]


def canonical_dbs(ref_di, alt_di):
    """Return canonical DBS type string (e.g., 'AC>CA').

    If ref_di is not in the canonical set, reverse complement both.
    """
    if ref_di in CANONICAL_REFS:
        return f"{ref_di}>{alt_di}"
    rc_ref = reverse_complement(ref_di)
    rc_alt = reverse_complement(alt_di)
    return f"{rc_ref}>{rc_alt}"


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
    """Open a VCF file, handling gzip compression."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def extract_dbs_from_vcf(vcf_path):
    """Extract DBS mutations from a somatic VCF.

    Strategy:
    1. Check for MNV records: REF and ALT are both 2bp (direct DBS).
    2. Check for adjacent SNV pairs: consecutive positions on the same chrom.

    Returns a list of canonical DBS type strings (e.g., ['AC>CA', 'AT>CG']).
    """
    snvs = []  # (chrom, pos, ref, alt) for SNVs
    dbs_list = []

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
            alt = fields[4].upper().split(",")[0]  # Take first ALT allele
            filt = fields[6]

            if filt != "PASS" and filt != ".":
                continue

            # Skip complex/symbolic alleles
            if not re.match(r"^[ACGT]+$", ref) or not re.match(r"^[ACGT]+$", alt):
                continue

            # Direct MNV record: both REF and ALT are 2bp
            if len(ref) == 2 and len(alt) == 2:
                dbs_list.append(canonical_dbs(ref, alt))
                continue

            # Collect SNVs for adjacent-pair detection
            if len(ref) == 1 and len(alt) == 1:
                snvs.append((chrom, pos, ref, alt))

    # Sort SNVs and find adjacent pairs
    snvs.sort(key=lambda x: (x[0], x[1]))
    for i in range(len(snvs) - 1):
        chrom1, pos1, ref1, alt1 = snvs[i]
        chrom2, pos2, ref2, alt2 = snvs[i + 1]
        if chrom1 == chrom2 and pos2 == pos1 + 1:
            ref_di = ref1 + ref2
            alt_di = alt1 + alt2
            dbs_list.append(canonical_dbs(ref_di, alt_di))

    return dbs_list


def build_dbs_matrix(dbs_list, cosmic_types, sample_id):
    """Build a DBS-78 count matrix DataFrame."""
    counts = {t: 0 for t in cosmic_types}
    for dbs_type in dbs_list:
        if dbs_type in counts:
            counts[dbs_type] += 1
        else:
            print(f"WARNING: DBS type '{dbs_type}' not in COSMIC-78, skipping", file=sys.stderr)

    df = pd.DataFrame(
        {"MutationType": list(counts.keys()), sample_id: list(counts.values())}
    )
    df.set_index("MutationType", inplace=True)
    return df


def prepare_signatures_db(cosmic_csv_path, tmp_dir):
    """Convert COSMIC DBS-78 CSV to TSV for SigProfilerAssignment."""
    df = pd.read_csv(cosmic_csv_path, index_col=0)
    df.columns = [c.strip() for c in df.columns]
    out_path = os.path.join(tmp_dir, "COSMIC_DBS78.tsv")
    df.to_csv(out_path, sep="\t", index=True, index_label="Type")
    return out_path, list(df.index)


def run_cosmic_fit_dbs(matrix_tsv_path, sig_db_tsv_path, output_dir):
    """Run SigProfilerAssignment for DBS-78."""
    from SigProfilerAssignment import Analyzer

    Analyzer.cosmic_fit(
        samples=matrix_tsv_path,
        output=output_dir,
        signature_database=sig_db_tsv_path,
        genome_build="GRCh38",
        input_type="matrix",
        context_type="DINUC",
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
        print(f"WARNING: Zero DBS activity for {sample_id}", file=sys.stderr)
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
                "ENTITY_STABLE_ID": f"mutational_signatures_contribution_{sig_id}",
                "NAME": f"{sig_id} ({category})" if category else sig_id,
                "DESCRIPTION": f"{etiology} ({main_effect})" if main_effect else etiology,
                sample_id: round(prop, 5),
            }
        )

    out_df = pd.DataFrame(rows).sort_values("ENTITY_STABLE_ID")
    out_df.to_csv(output_path, sep="\t", index=False)
    print(f"Written {len(rows)} DBS contribution rows to {output_path}")


def write_counts(matrix_df, sample_id, output_path):
    """Write DBS counts in cBioPortal GENERIC_ASSAY format."""
    rows = []
    for dbs_type, count in matrix_df[sample_id].items():
        ref = dbs_type.split(">")[0]
        alt = dbs_type.split(">")[1]
        rows.append(
            {
                "ENTITY_STABLE_ID": f"mutational_signatures_matrix_{ref}-{alt}",
                "NAME": dbs_type,
                sample_id: int(count),
            }
        )

    out_df = pd.DataFrame(rows).sort_values("ENTITY_STABLE_ID")
    out_df.to_csv(output_path, sep="\t", index=False)
    print(f"Written {len(rows)} DBS counts rows to {output_path}")


def main():
    args = parse_args()

    with tempfile.TemporaryDirectory() as tmp_dir:
        # Prepare COSMIC reference and get type list
        print(f"Preparing DBS signatures DB: {args.signatures_db}")
        sig_db_tsv, cosmic_types = prepare_signatures_db(args.signatures_db, tmp_dir)

        # Extract DBS from VCF
        print(f"Extracting DBS from VCF: {args.vcf}")
        dbs_list = extract_dbs_from_vcf(args.vcf)
        print(f"Found {len(dbs_list)} DBS mutations")

        # Build count matrix
        matrix_df = build_dbs_matrix(dbs_list, cosmic_types, args.sample)
        total_dbs = int(matrix_df[args.sample].sum())
        print(f"Total DBS counts: {total_dbs}")

        # Write counts file
        write_counts(matrix_df, args.sample, args.output_counts)

        if total_dbs == 0:
            print("WARNING: No DBS found, writing empty contribution file", file=sys.stderr)
            with open(args.output_contrib, "w") as fh:
                fh.write(f"ENTITY_STABLE_ID\tNAME\tDESCRIPTION\t{args.sample}\n")
            return

        # Write matrix for SigProfilerAssignment
        matrix_tsv = os.path.join(tmp_dir, "dbs_matrix.tsv")
        matrix_df.to_csv(matrix_tsv, sep="\t")

        # Run SigProfilerAssignment
        spa_outdir = os.path.join(tmp_dir, "spa_output")
        os.makedirs(spa_outdir)
        print("Running SigProfilerAssignment for DBS...")
        activities_path = run_cosmic_fit_dbs(matrix_tsv, sig_db_tsv, spa_outdir)

        # Write contribution file
        metadata = load_metadata(args.metadata, "dbs_id")
        write_contribution(activities_path, metadata, args.sample, args.output_contrib)


if __name__ == "__main__":
    main()
