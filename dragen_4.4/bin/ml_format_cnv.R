#!/usr/bin/env Rscript

# Load required libraries
library(tidyr)
library(dplyr)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
#args <- "/project/60005/shared/sub_projects/cbioportal_nextflow/nextflow-cbioportal/output/MoHQ-CM-3/data_cna_long.txt"

# Check if file path argument was provided
if (length(args) == 0) {
  stop("Error: Please provide the path to the input file as an argument.\nUsage: Rscript script.R /path/to/data_cna_long.txt")
}

input_file <- args[1]

# Verify the filename is exactly "data_cna_long.txt"
filename <- basename(input_file)
if (filename != "data_cna_long.txt") {
  stop("Error: Input file must be named 'data_cna_long.txt'. Provided file: ", filename)
}

# Check if file exists
if (!file.exists(input_file)) {
  stop("Error: File not found: ", input_file)
}

cat("Reading file:", input_file, "\n")

# Read the input file
data <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Reshape data from long to wide format
# Sample_Id values become column names
# Hugo_Symbol values become row names
wide_data <- data %>%
  select(Hugo_Symbol, Sample_Id, Value) %>%
  pivot_wider(names_from = Sample_Id, values_from = Value)

# Convert to data frame and set Hugo_Symbol as row names
wide_df <- as.data.frame(wide_data)
colnames(wide_df)[1] <- "Gene_ID"

# Write output to file
output_file <- "cnv_gene.tsv"
write.table(wide_df, output_file, sep = "\t", quote = FALSE, row.names = FALSE)

# Display dimensions and preview
cat("Output dimensions:", dim(wide_df)[1], "genes x", dim(wide_df)[2], "samples\n")
cat("\nFirst few rows and columns:\n")
print(head(wide_df[, 1:min(5, ncol(wide_df))]))

cat("\nOutput saved to:", output_file, "\n")
