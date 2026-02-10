#!/usr/bin/env Rscript

library(tidyverse)

#' Process copy number variation data for both standard and DNN input formats
#' @param base_dir Base directory containing the data, defaults to "data/GDC_TCGA"
#' @return List containing both processed CNV data frames (standard and DNN format)
#' @export
process_cnv_data <- function(input_file, output_file) { 
  # Read and process CNV data
  cnv_data <- read_tsv(input_file, show_col_types = FALSE)
  
  # Record initial dimensions
  initial_genes <- nrow(cnv_data)
  initial_samples <- ncol(cnv_data) - 1
  
  # Process the data - common steps
  cnv_data <- cnv_data %>%
    # Remove version numbers from Ensembl IDs
    mutate(Gene_ID = sub("\\.[0-9]+$", "", Gene_ID)) %>%
    # Remove genes that are NA for all samples
    filter(rowSums(!is.na(select(., -Gene_ID))) > 0)
  
  # Create DNN format (samples as rows)
  cnv_data_t <- cnv_data %>%
    column_to_rownames("Gene_ID") %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("Sample_ID")
 
  write_tsv(cnv_data_t, "cnv_genes.tsv")

}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Usage: Rscript ml_cnv_processor.R <input> <output>")
}

input <- args[1]

if (!file.exists(input)) {
  stop("Error: input file does not exist: ", input)
}

output <- args[2]

process_cnv_data(input, output)
