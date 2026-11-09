#!/usr/bin/env python3
"""
Write cBioPortal gene panel files from the pipeline's gene panel TSV (see filter_gene_panel.py).

Usage: write_gene_panels.py <gene_panel.tsv> <linking_file> <study_id> [--cnv-failed ID ...]

Writes to the current directory:
  gene_panels/gene_panel_<study_id>_<type>.txt  panel definitions (stable_id, description,
                                                 gene_list), one per type; import each with
                                                 importGenePanel.pl before the study
  data_gene_panel_matrix.txt                     every linking-file sample x profile -> panel
  meta_gene_panel_matrix.txt

The matrix marks every sample as profiled with the panel, so cBioPortal's frequency denominator
is the whole cohort instead of only the samples with data. Samples that failed CNV depth QC
(--cnv-failed, anonymised sample_id) get NA for cna: not profiled, as in the CNA case lists.
"""

import argparse
import os

from filter_gene_panel import read_panel

# panel type -> (cBioPortal profile stable_id from format_meta.py, description label)
PROFILES = {
    "mutation": ("mutations", "mutations"),
    "cna": ("cna", "copy-number alterations"),
    "sv": ("structural_variants", "structural variants / fusions"),
}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("gene_panel")
    parser.add_argument("linking_file")
    parser.add_argument("study_id")
    parser.add_argument("--cnv-failed", nargs="*", default=[])
    args = parser.parse_args()

    panel = read_panel(args.gene_panel)
    missing = [t for t in PROFILES if not panel[t]]
    if missing:
        raise SystemExit(f"ERROR: {args.gene_panel} lists no genes for: {', '.join(missing)}")

    panel_ids = {t: f"{args.study_id}_{t}" for t in PROFILES}
    os.makedirs("gene_panels", exist_ok=True)
    for t, (_, label) in PROFILES.items():
        genes = sorted(panel[t])
        with open(f"gene_panels/gene_panel_{panel_ids[t]}.txt", "w") as fh:
            fh.write(f"stable_id: {panel_ids[t]}\n")
            fh.write(f"description: {args.study_id} panel for {label} ({len(genes)} genes)\n")
            fh.write("gene_list:\t" + "\t".join(genes) + "\n")

    failed = set(args.cnv_failed)
    samples = []
    with open(args.linking_file) as fh:
        next(fh, None)
        for line in fh:
            cols = line.rstrip("\r\n").split("\t")
            if len(cols) >= 2 and cols[1].strip():
                samples.append((cols[0].strip(), cols[1].strip()))

    with open("data_gene_panel_matrix.txt", "w") as fh:
        fh.write("SAMPLE_ID\t" + "\t".join(p for p, _ in PROFILES.values()) + "\n")
        for anon_id, sample_id in samples:
            row = [sample_id]
            for t in PROFILES:
                row.append("NA" if t == "cna" and anon_id in failed else panel_ids[t])
            fh.write("\t".join(row) + "\n")

    with open("meta_gene_panel_matrix.txt", "w") as fh:
        fh.write(
            f"cancer_study_identifier: {args.study_id}\n"
            "genetic_alteration_type: GENE_PANEL_MATRIX\n"
            "datatype: GENE_PANEL_MATRIX\n"
            "data_filename: data_gene_panel_matrix.txt\n"
        )

    print(f"gene panel matrix: {len(samples)} samples, {len(failed & {a for a, _ in samples})} with cna = NA")


if __name__ == "__main__":
    main()
