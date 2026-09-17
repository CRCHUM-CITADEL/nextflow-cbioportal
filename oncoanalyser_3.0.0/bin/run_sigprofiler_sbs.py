#!/usr/bin/env python3
"""
Run SigProfilerAssignment cosmic_fit on a single sample's trinucleotide counts,
then convert Activities output to cBioPortal GENERIC_ASSAY contribution format.

Inputs:
  --snv_counts    : {subject}-T.sig.snv_counts.csv (BucketName, {subject}-T columns)
  --signatures_db : COSMIC_Human_SBS-96_GRCh38_v3.6.csv (extracted from zip)
  --metadata      : cosmic_sbs_metadata.tsv (sbs_id, category, etiology, main_effect)
  --sample        : sample ID for output column header
  --output        : output path for cBioPortal contribution file

Output:
  {sample}.data_mutational_signatures_contribution_SBS.txt
  Columns: ENTITY_STABLE_ID, NAME, DESCRIPTION, {sample}
"""

import argparse
import csv
import os
import sys
import tempfile

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--snv_counts", required=True, help="Path to .sig.snv_counts.csv")
    p.add_argument("--signatures_db", required=True, help="Path to COSMIC SBS-96 CSV")
    p.add_argument("--metadata", required=True, help="Path to cosmic_sbs_metadata.tsv")
    p.add_argument("--sample", required=True, help="Sample ID for output column")
    p.add_argument("--output", required=True, help="Output file path")
    return p.parse_args()


def bucket_to_sigprofiler(bucket):
    """Convert HMFTOOLS BucketName to SigProfiler context notation.

    C>A_ACA  ->  A[C>A]A
    """
    parts = bucket.split("_")
    if len(parts) != 2 or len(parts[1]) != 3:
        raise ValueError(f"Unexpected BucketName format: {bucket!r}")
    mutation, context = parts[0], parts[1]
    return f"{context[0]}[{mutation}]{context[2]}"


def convert_snv_counts_to_matrix(snv_counts_path, sample_id):
    """Read sig.snv_counts.csv and return a SigProfiler-compatible DataFrame."""
    df = pd.read_csv(snv_counts_path, index_col=0)
    if df.shape[1] != 1:
        raise ValueError(
            f"Expected 1 sample column in {snv_counts_path}, got {df.shape[1]}"
        )
    df.index = [bucket_to_sigprofiler(b) for b in df.index]
    df.index.name = "MutationType"
    df.columns = [sample_id]
    return df


def prepare_signatures_db(cosmic_csv_path, tmp_dir):
    """Convert COSMIC comma-separated CSV to tab-separated TSV for SigProfilerAssignment."""
    df = pd.read_csv(cosmic_csv_path, index_col=0)
    df.columns = [c.strip() for c in df.columns]
    out_path = os.path.join(tmp_dir, "COSMIC_SBS96.tsv")
    df.to_csv(out_path, sep="\t", index=True, index_label="Type")
    return out_path


def run_cosmic_fit(matrix_tsv_path, sig_db_tsv_path, output_dir):
    """Run SigProfilerAssignment cosmic_fit. Returns path to Activities file."""
    from SigProfilerAssignment import Analyzer

    Analyzer.cosmic_fit(
        samples=matrix_tsv_path,
        output=output_dir,
        signature_database=sig_db_tsv_path,
        genome_build="GRCh38",
        input_type="matrix",
        context_type="96",
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


def load_metadata(metadata_path):
    """Load cosmic_sbs_metadata.tsv -> dict keyed by sbs_id."""
    meta = {}
    with open(metadata_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            meta[row["sbs_id"]] = row
    return meta


def activities_to_cbioportal(activities_path, metadata, sample_id, output_path):
    """Convert SigProfilerAssignment Activities to cBioPortal contribution format."""
    acts = pd.read_csv(activities_path, sep="\t", index_col=0)

    if sample_id not in acts.index:
        sys.exit(
            f"ERROR: sample '{sample_id}' not found in Activities rows: "
            f"{list(acts.index)}"
        )

    row = acts.loc[sample_id]
    total = row.sum()

    if total == 0:
        print(
            f"WARNING: Total mutation count is 0 for {sample_id}. "
            "Writing empty contribution file.",
            file=sys.stderr,
        )
        with open(output_path, "w") as fh:
            fh.write(f"ENTITY_STABLE_ID\tNAME\tDESCRIPTION\t{sample_id}\n")
        return

    proportions = row / total
    proportions = proportions[proportions > 0]

    rows = []
    for sbs_id, prop in proportions.items():
        sbs_id = sbs_id.strip()
        meta = metadata.get(sbs_id, {})
        category = meta.get("category", "Unknown")
        etiology = meta.get("etiology", f"Mutational signature {sbs_id}")
        main_effect = meta.get("main_effect", "")

        entity_id = f"mutational_signatures_contribution_{sbs_id}"
        name = f"{sbs_id} ({category})" if category else sbs_id
        description = f"{etiology} ({main_effect})" if main_effect else etiology

        rows.append(
            {
                "ENTITY_STABLE_ID": entity_id,
                "NAME": name,
                "DESCRIPTION": description,
                sample_id: round(prop, 5),
            }
        )

    out_df = pd.DataFrame(rows)
    out_df.sort_values("ENTITY_STABLE_ID", inplace=True)
    out_df.to_csv(output_path, sep="\t", index=False)
    print(f"Written {len(rows)} signature rows to {output_path}")


def main():
    args = parse_args()

    with tempfile.TemporaryDirectory() as tmp_dir:
        # Step 1: Convert SNV counts to SigProfiler matrix
        print(f"Converting SNV counts: {args.snv_counts}")
        matrix_df = convert_snv_counts_to_matrix(args.snv_counts, args.sample)
        total_muts = int(matrix_df[args.sample].sum())
        print(f"Total mutations: {total_muts}")

        if total_muts < 50:
            print(
                f"WARNING: Only {total_muts} mutations for {args.sample}. "
                "Results may be unreliable.",
                file=sys.stderr,
            )

        matrix_tsv = os.path.join(tmp_dir, "matrix.tsv")
        matrix_df.to_csv(matrix_tsv, sep="\t")

        # Step 2: Prepare COSMIC reference (CSV -> TSV)
        print(f"Preparing signatures DB: {args.signatures_db}")
        sig_db_tsv = prepare_signatures_db(args.signatures_db, tmp_dir)

        # Step 3: Run SigProfilerAssignment
        spa_outdir = os.path.join(tmp_dir, "spa_output")
        os.makedirs(spa_outdir)
        print("Running SigProfilerAssignment cosmic_fit...")
        activities_path = run_cosmic_fit(matrix_tsv, sig_db_tsv, spa_outdir)

        # Step 4: Convert to cBioPortal format
        metadata = load_metadata(args.metadata)
        activities_to_cbioportal(activities_path, metadata, args.sample, args.output)


if __name__ == "__main__":
    main()
