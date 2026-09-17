#!/usr/bin/env Rscript

# Load required libraries
library(tidyr)
library(dplyr)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
#args <- "/project/60005/shared/sub_projects/cbioportal_nextflow/nextflow-cbioportal/output/MoHQ-CM-3/data_expression.txt"

# Check if file path argument was provided
if (length(args) == 0) {
  stop("Error: Please provide the path to the input file as an argument.\nUsage: Rscript script.R /path/to/data_expression.txt")
}

input_file <- args[1]

# Verify the filename is exactly "data_expression.txt"
filename <- basename(input_file)
if (filename != "data_expression.txt") {
  stop("Error: Input file must be named 'data_expression.txt'. Provided file: ", filename)
}

# Check if file exists
if (!file.exists(input_file)) {
  stop("Error: File not found: ", input_file)
}

cat("Reading file:", input_file, "\n")

# Read the input file
data <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

data <- data[,-2]
colnames(data)[1] <- "Gene_ID"

# Write output to file
output_file <- "expression_tpm.tsv"
write.table(data, output_file, sep = "\t", quote = FALSE, row.names = FALSE)

# Display dimensions and preview
cat("Output dimensions:", dim(data)[1], "genes x", dim(data)[2], "samples\n")
cat("\nFirst few rows and columns:\n")
print(head(data[, 1:min(5, ncol(data))]))

cat("\nOutput saved to:", output_file, "\n")
