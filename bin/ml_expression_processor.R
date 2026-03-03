#!/usr/bin/env Rscript

library(tidyverse)

#' Standardize expression matrix efficiently
#' @param expression_mat Numeric matrix of expression values
#' @return Standardized expression matrix
#' @keywords internal
standardize_expression <- function(expression_mat) {
  # Transpose, use highly optimized base scale(), then transpose back.
  scaled_mat <- scale(t(expression_mat))
  
  # Handle potential NaN values introduced if a gene had 0 variance 
  scaled_mat[is.na(scaled_mat)] <- 0 
  
  return(t(scaled_mat))
}

#' Process TPM expression data for Deep Learning
#' @param input_file Path to input TPM TSV
#' @param min_tpm Minimum TPM value threshold (default: 1)
#' @param min_samples Minimum number of samples where gene should be expressed (default: 3)
#' @export
process_expression_data <- function(input_file, min_tpm = 1, min_samples = 3) {
  output_file_log2 <- "expression_processed_log2.tsv"
  output_file_standardized <- "expression_processed_standardized.tsv"

  message("Reading expression data...")
  expression_data <- read_tsv(input_file, show_col_types = FALSE)

  initial_genes <- nrow(expression_data)
  initial_samples <- ncol(expression_data) - 1 

  message("Processing and normalizing...")
  expression_processed <- expression_data %>%
    # NEW: Drop any corrupted/blank Gene IDs right away
    filter(!is.na(Gene_ID) & Gene_ID != "") %>%
    # 1. Clean Gene IDs
    mutate(Gene_ID = sub("\\.[0-9]+$", "", Gene_ID)) %>%
    # 2. Smart Deduplication
    mutate(row_mean = rowMeans(select(., -Gene_ID), na.rm = TRUE)) %>%
    arrange(desc(row_mean)) %>%
    distinct(Gene_ID, .keep_all = TRUE) %>%
    select(-row_mean) %>%
    # 3. Filter low-expression genes (NEW: Safely handle NAs with na.rm = TRUE)
    filter(rowSums(select(., -Gene_ID) >= min_tpm, na.rm = TRUE) >= min_samples) %>%
    # 4. ML Best Practice: log2(x + 1)
    mutate(across(-Gene_ID, ~log2(.x + 1)))

  # Convert to matrix for fast operations
  Gene_IDs <- expression_processed$Gene_ID
  expr_mat <- as.matrix(expression_processed %>% select(-Gene_ID))
  rownames(expr_mat) <- Gene_IDs

  # 5. Fast Zero-Variance Filtering (NEW: Safely drop NA variances)
  gene_vars <- apply(expr_mat, 1, var, na.rm = TRUE)
  expr_mat <- expr_mat[!is.na(gene_vars) & gene_vars > 0, , drop = FALSE]
  Gene_IDs <- rownames(expr_mat)

  filtered_genes <- nrow(expr_mat)
  removed_genes <- initial_genes - filtered_genes

  # 6. Standardize (Z-score across samples)
  message("Applying Z-score standardization...")
  std_mat <- standardize_expression(expr_mat)

  # 7. Format for ML Dataloaders (Samples as Rows, Genes as Columns)
  expression_final <- as.data.frame(round(t(expr_mat), 4)) %>%
    rownames_to_column("sample_id")
    
  standardized_final <- as.data.frame(round(t(std_mat), 4)) %>%
    rownames_to_column("sample_id")

  # Clean up sample names
  expression_final$sample_id <- gsub("\\.", "-", expression_final$sample_id)
  standardized_final$sample_id <- gsub("\\.", "-", standardized_final$sample_id)

  message("Saving outputs...")
  write_tsv(expression_final, output_file_log2)
  write_tsv(standardized_final, output_file_standardized)

  # Summary
  message("\n=== Expression Processing Summary ===")
  message(sprintf("- Initial genes: %d", initial_genes))
  message(sprintf("- Removed %d genes (duplicates/low expression/zero var)", removed_genes))
  message(sprintf("- Retained %d biologically active genes", filtered_genes))
  message(sprintf("- Number of samples: %d", initial_samples))
  message("- Transformations:")
  message("  1. Smart Deduplication (Highest mean retained)")
  message("  2. log2(TPM + 1) -> Preserves true biological zeros")
  message("  3. Z-score standard scaling")
  message("=====================================\n")
}

args = commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Usage : Rscript ml_expression_processor.R [path_to_expression_data]")
}

input_file <- args[1]
if (!file.exists(input_file)) {
  stop("Error: Input file does not exist: ", input_file)
}

process_expression_data(input_file)
