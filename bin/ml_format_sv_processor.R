# Gene Fusion Annotation Script for cBioPortal Data
# This script filters known cancer fusions and creates a sample x fusion matrix

# Load required libraries
library(tidyverse)
library(readxl)

#' Validate fusion data structure and content
#' @param fusion_data Processed fusion data frame
#' @return Logical indicating if validation passed (invisible)
#' @keywords internal
validate_fusion_data <- function(fusion_data) {
  # Check if data frame is empty
  if (nrow(fusion_data) == 0) {
    stop("No samples remained after processing")
  }

  # Check if sample column exists
  if (!"sample" %in% colnames(fusion_data)) {
    stop("Sample column not found in fusion data")
  }

  # Check for duplicate samples
  if (any(duplicated(fusion_data$sample))) {
    stop("Duplicate samples found in fusion data")
  }

  # Check if any non-numeric values in fusion columns
  non_numeric_check <- fusion_data %>%
    select(-sample) %>%
    sapply(function(x) all(is.numeric(x)))

  if (!all(non_numeric_check)) {
    stop("Non-numeric values found in fusion data")
  }

  invisible(TRUE)
}

#' Encode fusions using binary representation
#' @keywords internal
encode_binary_fusions <- function(data, samples, known_fusions) {
  # Create sparse matrix of fusion presence/absence
  fusion_matrix <- matrix(0, nrow = length(samples), ncol = length(known_fusions),
                          dimnames = list(samples, known_fusions))

  # Fill matrix with 1's where fusions exist
  for (i in seq_len(nrow(data))) {
    if (data$Fusion[i] %in% known_fusions) {
      fusion_matrix[data$Sample_Id[i], data$Fusion[i]] <- 1
    }
  }
#data <- data[data$Fusion %in% known_fusions,]

  as_tibble(fusion_matrix, rownames = "sample")
}


#' Process fusion data for deep learning input
#' @param base_dir Base directory containing the data
#' @return List containing different fusion encoding matrices
#' @export
process_fusion_data <- function(input_file, output_file, known_fusions ) {

  # Read fusion data
  data <- read_tsv(input_file, show_col_types = FALSE)
  data$Fusion <- paste0(data$Site1_Hugo_Symbol,"-",data$Site2_Hugo_Symbol)

  # Get unique samples and genes
  samples <- unique(data$Sample_Id)

  # Create different encodings
  binary_matrix <- encode_binary_fusions(data, samples, known_fusions)

  # Validate matrix
  validate_fusion_data(binary_matrix)

  # Save binary matrix
  write_tsv(binary_matrix, output_file) 

}

# ============================================================================
# CONFIGURATION & INPUT VALIDATION
# ============================================================================

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Usage: Rscript fusion_annotator.R <path_to_data_fusions.txt> <known_chimerDB_xlsx> <cosmic v103 GRCh38>")
}

input_file <- args[1]
known_fusions <- args[2] # Get known fusions from ChimerDB (https://www.kobic.re.kr/chimerdb/download)
cosmic_data <- args[3] # Get known fusions from Cosmic v103 GRCh38

# Check that file exists
if (!file.exists(input_file)) {
  stop("Error: Input file does not exist: ", input_file)
}

if (!file.exists(known_fusions)) {
  stop("Error : known fusions file does not exist: ",known_fusions)
}

if (!file.exists(cosmic_data)) {
  stop("Error : cosmic data file does not exist: ",cosmic_data)
}

cat("Gene Fusion Annotation Pipeline\n")
cat("================================\n\n")
cat("Input file:", input_file, "\n")

# Output file (same directory as input)
output_file <- "filtered_fusions_matrix.tsv"

cosmic_fusions <- cosmic_data[,!is.null(cosmic_data$COSMIC_FUSION_ID)]
cosmic_fusions <- unique(paste0(cosmic_db[,"FIVE_PRIME_GENE_SYMBOL"],"-",cosmic_data[,"THREE_PRIME_GENE_SYMBOL"]))

chimerKB_fusions <- unique(chimerKB_db$Fusion_pair)

all_fusions <- c(cosmic_fusions, chimerKB_fusions)


