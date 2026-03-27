#!/usr/bin/env Rscript

library(tidyverse)

#' Process copy number variation data for deep learning input
#' @param input_file Path to input CNV data
#' @param output_file Path to output processed data
#' @return invisible data frame
#' @export
process_cnv_data <- function(input_file, output_file) {
  # Read data
  message("Reading CNV data...")
  cnv_data <- read_tsv(input_file, show_col_types = FALSE)
  
  initial_genes <- nrow(cnv_data)
  initial_samples <- ncol(cnv_data) - 1
  
  # 1. Clean Gene IDs and elegantly handle duplicates
  cnv_cleaned <- cnv_data %>%
    mutate(Gene_ID = sub("\\.[0-9]+$", "", Gene_ID)) %>%
    group_by(Gene_ID) %>%
    # If multiple transcripts collapse to the same gene, take the first one 
    # to guarantee unique row names without crashing
    summarise(across(everything(), first), .groups = "drop")
    
  dedup_genes <- nrow(cnv_cleaned)
  
  # Convert to base R matrix for fast mathematical operations
  cnv_matrix <- cnv_cleaned %>%
    column_to_rownames("Gene_ID") %>%
    as.matrix()
  
  # Count missing values for the summary (NAs are retained)
  na_count <- sum(is.na(cnv_matrix))
  
  # 2. Filter uninformative (Zero-Variance or All-NA) genes
  # Keep only genes that have at least one CNV event (!= 0) in the cohort.
  # na.rm = TRUE ensures that NAs are ignored during this check, so a gene
  # with a mix of 0s and NAs will be correctly dropped, while a gene with 
  # at least one valid 1, -1, 2, or -2 will be kept.
  print(head(cnv_matrix))
  informative_mask <- rowSums(cnv_matrix != 0, na.rm = TRUE) > 0
  
  # drop = FALSE prevents R from coercing to a vector if only 1 gene remains
  print("drop false")
  cnv_filtered <- cnv_matrix[informative_mask, , drop = FALSE]
  
  final_genes <- nrow(cnv_filtered)
  
  # 3. Transpose for DNN (Samples as rows, Genes as columns)
  print("transpose")
  cnv_dnn <- t(cnv_filtered) %>%
    as.data.frame() %>%
    rownames_to_column("Sample_ID")
    
  print("writing")
  # 4. Save output using the actual provided argument
  write.table(cnv_dnn,
              file = output_file,
              sep = "\t",
              row.names = FALSE,
              col.names = TRUE,
              quote = FALSE)

  # write_tsv(cnv_dnn, output_file)
  
  # Print processing summary
  message("\nCNV Processing Summary:")
  message(sprintf("- Initial genes: %d", initial_genes))
  message(sprintf("- Unique genes after Ensembl ID stripping: %d", dedup_genes))
  message(sprintf("- Missing values retained (NAs): %d", na_count))
  message(sprintf("- Uninformative (all-zero or all-NA) genes removed: %d", dedup_genes - final_genes))
  message(sprintf("- Final features (genes) retained: %d", final_genes))
  message(sprintf("- Number of samples: %d", initial_samples))
  message(sprintf("- Output successfully saved to: %s\n", output_file))
  
  # return(invisible(cnv_dnn))
}

args <- commandArgs(trailingOnly = TRUE)

# Enforce exactly 1 argument
if (length(args) != 1) {
  stop("Usage: Rscript ml_cnv_processor.R <input_file>")
}

input <- args[1]

if (!file.exists(input)) {
  stop("Error: input file does not exist: ", input)
}

process_cnv_data(input, "cnv_genes.tsv")
