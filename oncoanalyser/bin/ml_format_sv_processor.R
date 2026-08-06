#!/usr/bin/env Rscript

# Gene Fusion Annotation Script for Deep Learning Data Prep
# Vectorized, ML-optimized pipeline to create a sparse binary matrix of validated, high-confidence fusions.

library(tidyverse)
library(readxl)

#' Process fusion data for deep learning input
#' @param input_file Path to raw fusion data
#' @param output_file Path to save the processed matrix
#' @param known_fusions Vector of validated fusion strings (Site1-Site2)
#' @param min_support_reads Minimum combined read count to trust the fusion
#' @export
process_fusion_data <- function(input_file, output_file, known_fusions, min_support_reads = 1) {

  message("Reading raw fusion data...")
  data <- read_tsv(input_file, show_col_types = FALSE)

  # 1. Capture ALL unique samples to prevent dropping patients
  all_samples <- tibble(Sample_Id = unique(data$Sample_Id))
  initial_sample_count <- nrow(all_samples)

  # 2. Strict QC and Database Filtering
  message(sprintf("Filtering for RNA support and >= %d supporting reads...", min_support_reads))

  valid_fusions_data <- data %>%
    mutate(Fusion = paste0(Site1_Hugo_Symbol, "-", Site2_Hugo_Symbol)) %>%
    # Convert read counts to numeric and handle potential NAs
    mutate(
      Split_Reads = replace_na(as.numeric(Tumor_Split_Read_Count), 0),
      Paired_Reads = replace_na(as.numeric(Tumor_Paired_End_Read_Count), 0),
      Total_Support = Split_Reads + Paired_Reads
    ) %>%
    # Apply Quality Control Gates
    filter(RNA_Support == "Yes" | RNA_Support == "TRUE") %>%
    filter(Total_Support >= min_support_reads) %>%
    # Cross-reference with COSMIC/ChimerDB
    filter(Fusion %in% known_fusions) %>%
    mutate(Present = 1) %>%
    distinct(Sample_Id, Fusion, Present)

  active_fusions_count <- length(unique(valid_fusions_data$Fusion))

  # 3. Vectorized Matrix Creation
  message("Pivoting to feature matrix...")
  fusion_matrix <- valid_fusions_data %>%
    pivot_wider(names_from = Fusion, values_from = Present, values_fill = 0)

  # 4. Guarantee all original samples are present
  final_matrix <- all_samples %>%
    left_join(fusion_matrix, by = "Sample_Id") %>%
    rename(sample_id = Sample_Id) # Standardize to sample_id for easier merging downstream

  # Clean up sample names to match other modalities
  final_matrix$sample_id <- gsub("\\.", "-", final_matrix$sample_id)

  # Convert NAs to 0
  final_matrix[is.na(final_matrix)] <- 0

  # Save the optimized matrix
  write_tsv(final_matrix, output_file)

  # Summary
  message("\n=== Fusion Processing Summary ===")
  message(sprintf("- Total samples processed: %d", initial_sample_count))
  message(sprintf("- Active, high-confidence fusions retained: %d", active_fusions_count))
  message("- ML Transformations:")
  message(sprintf("  1. QC Gate: Required RNA_Support and >= %d total reads", min_support_reads))
  message("  2. DB Gate: Cross-referenced with COSMIC/ChimerDB")
  message("  3. Matrix Generation: Zero-padded patients with no valid fusions")
  message("=================================\n")
}

# ============================================================================
# CONFIGURATION & INPUT VALIDATION
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Usage: Rscript fusion_annotator.R <raw_fusions.tsv> <chimerDB.xlsx> <cosmic_v103.tsv> [min_reads]")
}

input_file <- args[1]
known_fusions_file <- args[2]
cosmic_data_file <- args[3]
# Allow user to override the read threshold via command line, default to 1
min_reads <- if (length(args) >= 4) as.numeric(args[4]) else 1

if (!file.exists(input_file)) stop("Error: Input file missing: ", input_file)
if (!file.exists(known_fusions_file)) stop("Error: ChimerDB file missing: ", known_fusions_file)
if (!file.exists(cosmic_data_file)) stop("Error: COSMIC file missing: ", cosmic_data_file)

cat("High-Fidelity Fusion Annotation Pipeline\n")
cat("========================================\n")

# 1. Load Databases
cat("Loading reference databases...\n")
chimerKB_db <- read_excel(known_fusions_file)
chimerKB_fusions <- unique(chimerKB_db$Fusion_pair)

print(head(chimerKB_db))
print(head(chimerKB_fusions))

cosmic_data_df <- read_tsv(cosmic_data_file, show_col_types = FALSE)
print(head(cosmic_data_df))
cosmic_fusions <- cosmic_data_df %>%
  filter(!is.na(COSMIC_FUSION_ID)) %>%
  mutate(Fusion = paste0(FIVE_PRIME_GENE_SYMBOL, "-", THREE_PRIME_GENE_SYMBOL)) %>%
  pull(Fusion) %>%
  unique()

all_known_fusions <- unique(c(cosmic_fusions, chimerKB_fusions))
print("all known fusions")
print(head(all_known_fusions))

# 2. Process Data
output_file <- "filtered_fusions_matrix.tsv"
process_fusion_data(input_file, output_file, all_known_fusions, min_reads)
