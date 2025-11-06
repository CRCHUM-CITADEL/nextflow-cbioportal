#!/usr/bin/env Rscript

# Script to format Dragen TPM values with gene symbols
# This script converts a Dragen file to a gene symbol and TPM value file

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(optparse)
})

# Define command line options
option_list <- list(
  # Input files
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Path to input Dragen expression file [required]"),
  make_option(c("-g", "--gene_map"), type="character", default=NULL,
              help="Path to gene_id_to_name.tsv file [required]"),
  make_option(c("-s", "--sample_id"), type="character", default="Sample1",
              help="Sample identifier [default=%default]"),

  # Output file
  make_option(c("-o", "--output"), type="character", default="expression_data.txt",
              help="Path to output expression data file [default=%default]"),

  # Processing options
  make_option("--min_tpm", type="double", default=0,
              help="Minimum TPM threshold (0 = no filtering) [default=%default]"),
  make_option("--trim_version", action="store_true", default=TRUE,
              help="Trim version numbers from Ensembl IDs [default=%default]")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list=option_list,
                           description="Converts Dragen TPM values to gene symbol format",
                           usage="usage: %prog [options]")
opt <- parse_args(opt_parser)

# Check for required arguments
if (is.null(opt$input) || is.null(opt$gene_map)) {
  print_help(opt_parser)
  stop("Both --input and --gene_map arguments are required.", call.=FALSE)
}

# Display configuration
cat("Configuration:\n")
cat("- Input file: ", opt$input, "\n")
cat("- Gene map file: ", opt$gene_map, "\n")
cat("- Sample ID: ", opt$sample_id, "\n")
cat("- Output file: ", opt$output, "\n")
cat("- Min TPM threshold: ", opt$min_tpm, "\n")
cat("- Trim Ensembl versions: ", opt$trim_version, "\n")
cat("\n")

# Step 1: Read input files
cat("Reading expression data...\n")
expr_data <- read_tsv(opt$input, show_col_types = FALSE)

cat("Reading gene ID to name mapping...\n")
# Read the gene_id_to_name.tsv file without headers
gene_map <- read_tsv(opt$gene_map,
                     col_names = c("ensembl_id", "gene_symbol"),
                     show_col_types = FALSE)

# Step 2: Trim Ensembl IDs if necessary
if (opt$trim_version) {
  cat("Trimming version numbers from Ensembl IDs...\n")
  expr_data <- expr_data %>%
    mutate(trimmed_ensembl = str_replace(Name, "\\.[0-9]+$", ""))
} else {
  expr_data$trimmed_ensembl <- expr_data$Name
}

# Step 3: Map Ensembl IDs to gene symbols
cat("Mapping Ensembl IDs to gene symbols...\n")
expression_with_symbols <- expr_data %>%
  inner_join(gene_map, by = c("trimmed_ensembl" = "ensembl_id"))

cat("- Successfully mapped", nrow(expression_with_symbols[!is.na(expression_with_symbols$gene_symbol),]), "out of", nrow(expr_data), "genes\n")

# Step 4: Format output
cat("Formatting output...\n")
formatted_output <- expression_with_symbols %>%
  select(gene_symbol, TPM) %>%
  rename(Hugo_Symbol = gene_symbol, !!opt$sample_id := TPM)

# Step 5: Handle duplicate Hugo symbols (keeping highest expression)
cat("Resolving duplicate Hugo symbols...\n")
duplicate_count <- formatted_output %>%
  count(Hugo_Symbol) %>%
  filter(n > 1) %>%
  nrow()

cat("- Found", duplicate_count, "duplicate gene symbols\n")

formatted_output <- formatted_output %>%
  group_by(Hugo_Symbol) %>%
  arrange(desc(!!sym(opt$sample_id))) %>%
  slice(1) %>%
  ungroup()

# Step 6: Apply minimum TPM threshold if specified
if (opt$min_tpm > 0) {
  cat("Applying minimum TPM threshold of", opt$min_tpm, "...\n")
  pre_filter_count <- nrow(formatted_output)
  formatted_output <- formatted_output %>%
    filter(!!sym(opt$sample_id) >= opt$min_tpm)
  cat("- Removed", pre_filter_count - nrow(formatted_output), "genes below threshold\n")
}

# Step 7: Write data to file
cat("Writing expression data to", opt$output, "...\n")
write_tsv(formatted_output, opt$output)

# Summary
cat("\nComplete!\n")
cat("- Output file:", opt$output, "\n")
cat("- Total genes in output:", nrow(formatted_output), "\n")
