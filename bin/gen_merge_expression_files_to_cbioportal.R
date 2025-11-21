#!/usr/bin/env Rscript
# Script to merge multiple single-sample expression files into one matrix for cBioPortal
# Takes individual sample expression files and combines them into a matrix with genes as rows

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(optparse)
  library(purrr)
  library(stringr)
})

# Define command line options
option_list <- list(
  # Input parameters
  make_option(c("-i", "--input_files"), type="character", default=NULL,
              help="Comma-separated list of expression file paths [required]"),
  make_option(c("--sample_ids"), type="character", default=NULL,
              help="Comma-separated list of sample IDs (overrides extraction from filenames)"),
  make_option("--use_file_numbers", action="store_true", default=FALSE,
              help="Use sequential numbers as sample IDs instead of extracting from filenames [default=%default]"),

  # Output parameters
  make_option(c("-o", "--output_file"), type="character", default="merged_expression.txt",
              help="Output file for merged expression matrix [default=%default]"),
  make_option(c("-m", "--output_meta"), type="character", default="meta_merged_expression.txt",
              help="Output metadata file [default=%default]"),

  # Processing options
  make_option(c("-f", "--fill_missing"), type="character", default="NA",
              help="Value to use for missing genes (NA or 0) [default=%default]"),
  make_option(c("-s", "--study_id"), type="character", default="my_study",
              help="Cancer study identifier for metadata file [default=%default]"),
  make_option("--sample_column", action="store_true", default=FALSE,
              help="Extract sample names from column headers instead of filenames [default=%default]"),
  make_option("--sample_pattern", type="character", default="(.*?)_expression\\.tsv$",
              help="Regex pattern to extract sample ID from filename [default=%default]"),
  make_option("--strict", action="store_true", default=FALSE,
              help="Exit with error if gene lists are not identical [default=%default]")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list=option_list,
                           description="Merges multiple expression files into a single matrix for cBioPortal",
                           usage="usage: %prog [options]")
opt <- parse_args(opt_parser)

# Check for required arguments
if (is.null(opt$input_files)) {
  print_help(opt_parser)
  stop("--input_files argument is required.", call.=FALSE)
}

# Validate options
if (!opt$fill_missing %in% c("NA", "0")) {
  warning("fill_missing should be 'NA' or '0'. Using 'NA'.")
  opt$fill_missing <- "NA"
}

# Parse input file list
expression_files <- strsplit(opt$input_files, ",")[[1]]
expression_files <- trimws(expression_files)  # Remove any whitespace

# Validate that files exist
missing_files <- expression_files[!file.exists(expression_files)]
if (length(missing_files) > 0) {
  stop("The following input files do not exist:\n  ",
       paste(missing_files, collapse="\n  "), call.=FALSE)
}

if (length(expression_files) == 0) {
  stop("No expression files provided.", call.=FALSE)
}

# Display configuration
cat("Configuration:\n")
cat("- Number of input files: ", length(expression_files), "\n")
cat("- Output matrix file: ", opt$output_file, "\n")
cat("- Output meta file: ", opt$output_meta, "\n")
cat("- Fill missing values with: ", opt$fill_missing, "\n")
cat("- Cancer study ID: ", opt$study_id, "\n")
cat("\n")

cat("Found", length(expression_files), "expression files:\n")
for (file in expression_files) {
  cat("  -", basename(file), "\n")
}
cat("\n")

# Step 2: Get sample IDs (either from filenames, manually specified, or use numbers)
if (!is.null(opt$sample_ids)) {
  # Option 1: Use manually specified sample IDs
  cat("Using manually specified sample IDs\n")
  sample_ids <- strsplit(opt$sample_ids, ",")[[1]]
  sample_ids <- trimws(sample_ids)  # Remove any whitespace

  # Check if count matches
  if (length(sample_ids) != length(expression_files)) {
    stop("Number of specified sample IDs (", length(sample_ids),
         ") does not match number of files (", length(expression_files), ")", call.=FALSE)
  }

  cat("Sample IDs:\n")
  for (i in 1:length(sample_ids)) {
    cat(sprintf("  - File %d: %s → %s\n", i, basename(expression_files[i]), sample_ids[i]))
  }

} else if (opt$use_file_numbers) {
  # Option 2: Use sequential numbers as sample IDs
  cat("Using sequential numbers as sample IDs\n")
  sample_ids <- paste0("Sample", seq_along(expression_files))

  cat("Sample IDs:\n")
  for (i in 1:length(sample_ids)) {
    cat(sprintf("  - File %d: %s → %s\n", i, basename(expression_files[i]), sample_ids[i]))
  }

} else {
  # Option 3: Simply use everything before "_expression.tsv" as the sample ID
  cat("Using filename prefix before '_expression.tsv' as sample ID\n")

  sample_ids <- basename(expression_files) %>%
    str_replace("\\.tpm\\.tsv$", "")

  cat("Sample IDs:\n")
  for (i in 1:length(sample_ids)) {
    cat(sprintf("  - %s → %s\n", basename(expression_files[i]), sample_ids[i]))
  }
}

# Ensure sample IDs are unique
if (length(unique(sample_ids)) != length(sample_ids)) {
  stop("Extracted sample IDs are not unique. Check your sample_pattern.", call.=FALSE)
}

# Step 3: Read all expression files and create a merged data frame
cat("Reading and merging expression files...\n")

# Function to read a single expression file
read_expr_file <- function(file_path, sample_id) {
  df <- read_tsv(file_path, show_col_types = FALSE)

  # If using column header for sample name
  if (opt$sample_column) {
    # Keep only Hugo_Symbol and the sample column (the second column)
    df <- df %>%
      select(Hugo_Symbol, Entrez_Gene_Id, 3)

    # Rename the second column to the sample ID from filename
    colnames(df)[3] <- sample_id
  } else {
    # Assume standard format with Hugo_Symbol and sample ID as column names
    # Rename the second column to the sample ID from filename
    colnames(df)[3] <- sample_id
  }

  return(df)
}

# Read all files
expression_list <- map2(expression_files, sample_ids, read_expr_file)

# Check if all samples have identical gene lists
cat("Checking gene list consistency across samples...\n")
gene_lists <- map(expression_list, function(df) df$Hugo_Symbol)
all_identical <- all(map_lgl(gene_lists, function(genes) setequal(genes, gene_lists[[1]])))

if (!all_identical) {
  cat("\n")
  cat("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
  cat("!!! WARNING: Gene lists are NOT identical across all samples                  !!!\n")
  cat("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
  cat("\n")

  # Count genes in each sample
  gene_counts <- map_int(gene_lists, length)
  names(gene_counts) <- sample_ids

  cat("Gene counts by sample:\n")
  for (i in 1:length(gene_counts)) {
    cat(sprintf("  - %-20s: %d genes\n", names(gene_counts)[i], gene_counts[i]))
  }

  # Find common genes across all samples
  common_genes <- Reduce(intersect, gene_lists)
  cat("\nCommon genes across all samples:", length(common_genes), "\n")

  # Find unique genes (in any sample but not all)
  all_genes <- unique(unlist(gene_lists))
  unique_genes <- setdiff(all_genes, common_genes)
  cat("Genes present in some but not all samples:", length(unique_genes), "\n")

  # Show example of unique genes
  if (length(unique_genes) > 0) {
    sample_unique <- unique_genes[1:min(5, length(unique_genes))]
    cat("Examples of genes not present in all samples:\n")

    for (gene in sample_unique) {
      present_in <- sample_ids[map_lgl(gene_lists, function(genes) gene %in% genes)]
      missing_from <- setdiff(sample_ids, present_in)

      cat(sprintf("  - %s: present in %d samples, missing from %d samples\n",
                 gene, length(present_in), length(missing_from)))

      if (length(missing_from) <= 5) {
        cat("     Missing from:", paste(missing_from, collapse=", "), "\n")
      }
    }
  }

  cat("\nPossible causes:\n")
  cat("  - Different gene filtering criteria used across samples\n")
  cat("  - Some samples processed with different annotation versions\n")
  cat("  - Missing data in some samples\n")

  cat("\nRecommendation:\n")
  cat("  - Check that all samples were processed with the same parameters\n")
  cat("  - Verify that annotation files are consistent\n")
  cat("  - Consider re-running the Kallisto-to-cBioPortal script with identical settings\n")
  cat("\nProceeding with merge operation using all genes (union)...\n\n")
} else {
  cat("All samples have identical gene lists with", length(gene_lists[[1]]), "genes each\n")
}

# Merge all data frames by Hugo_Symbol
cat("Performing merge operation...\n")
merged_expression <- expression_list[[1]]

for (i in 2:length(expression_list)) {
  merged_expression <- merged_expression %>%
    full_join(expression_list[[i]], by = c("Hugo_Symbol","Entrez_Gene_Id"))

  cat("  - Merged", i, "of", length(expression_list), "files\n")
}

# After the merge loop, reorder columns to ensure deterministic behaviour
sample_ids_sorted <- sort(sample_ids)
merged_expression <- merged_expression %>%
  select(Hugo_Symbol, Entrez_Gene_Id, all_of(sample_ids_sorted))

# Step 4: Replace missing values
if (opt$fill_missing == "0") {
  cat("Replacing missing values with 0 in numeric columns...\n")
  # Replace NA values only in numeric columns (all except first and second column)
  for (col in 3:ncol(merged_expression)) {
    col_name <- names(merged_expression)[col]
    merged_expression[[col_name]][is.na(merged_expression[[col_name]])] <- 0
  }
} else {
  cat("Keeping missing values as NA...\n")
}

# Get some statistics
total_genes <- nrow(merged_expression)
total_samples <- ncol(merged_expression) - 1
cat("Merged matrix has", total_genes, "genes and", total_samples, "samples\n")

# Step 5: Write merged expression matrix to file
cat("Writing merged expression matrix to", opt$output_file, "...\n")
write_tsv(merged_expression, opt$output_file)

# Summary
cat("\nComplete!\n")
cat("- Merged expression matrix file:", opt$output_file, "\n")
cat("- Total genes:", total_genes, "\n")
cat("- Total samples:", total_samples, "\n")

