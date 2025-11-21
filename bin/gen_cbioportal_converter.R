#!/usr/bin/env Rscript

# CNVkit/Dragen CNV to cBioPortal Converter
# This script converts CNVkit or Dragen CNV output to cBioPortal format

library(optparse)
library(data.table)

# Parse command line arguments
option_list <- list(
  make_option(c("-v", "--vcf"), type="character", help="Input VCF file from CNVkit or Dragen"),
  make_option(c("-t", "--tsv"), type="character", help="Input TSV file with gene-level CNV data"),
  make_option(c("-s", "--sample_id"), type="character", help="Tumor sample id"),
  #make_option(c("-s", "--study_id"), type="character", help="Cancer study identifier"),
  make_option(c("-o", "--output_dir"), type="character", default=".", help="Output directory [default: %default]"),
  make_option(c("--vcf_source"), type="character", default="auto", help="VCF source: 'cnvkit', 'dragen', or 'auto' for automatic detection [default: %default]")
)

opt_parser <- OptionParser(option_list=option_list,
                          description="Convert CNVkit or Dragen CNV output to cBioPortal formats")
opt <- parse_args(opt_parser)

# Check required arguments
if (is.null(opt$vcf) || is.null(opt$tsv) || is.null(opt$sample_id)) {
  cat("ERROR: Missing required arguments\n")
  print_help(opt_parser)
  quit(status = 1)
}

# Create output directory if it doesn't exist
if (!dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir, recursive = TRUE)
}

# Function to convert fold change to discrete copy number value
# MODIFIED to classify any gain with fold_change >= 2 as high-level amplification (2)
fold_change_to_discrete <- function(fold_change) {
  # Estimate copy number (assuming fold_change is relative to diploid state)
  estimated_cn <- 2 * fold_change

  if (estimated_cn < 0.5) {
    return(-2)  # homozygous deletion
  } else if (estimated_cn < 1.3) {
    return(-1)  # hemizygous deletion
  } else if (estimated_cn <= 2.7) {
    return(0)   # neutral
  } else if (estimated_cn < 4) {  # Changed from <= 5 to < 4 (fold_change < 2)
    return(1)   # gain
  } else {
    return(2)   # high-level amplification (4 copies and more)
  }
}

# Function to convert fold change to log2 ratio
fold_change_to_log2 <- function(fold_change) {
  return(log2(fold_change))
}

# Function to calculate fold change from log2 ratio
log2_to_fold_change <- function(log2_ratio) {
  return(2^log2_ratio)
}

# Function to detect VCF source
detect_vcf_source <- function(vcf_file) {
  vcf_lines <- readLines(vcf_file, n = 100)  # Read first 100 lines

  # Check for DRAGEN specific header
  if (any(grepl("##DRAGENCommandLine", vcf_lines, fixed = TRUE)) ||
      any(grepl("##source=DRAGEN_CNV", vcf_lines, fixed = TRUE))) {
    return("dragen")
  }

  # Check for CNVkit specific header
  if (any(grepl("##source=CNVkit", vcf_lines, fixed = TRUE))) {
    return("cnvkit")
  }

  # Check INFO fields
  if (any(grepl("##INFO=<ID=FOLD_CHANGE,", vcf_lines, fixed = TRUE))) {
    return("cnvkit")
  }

  if (any(grepl("##FORMAT=<ID=SM,", vcf_lines, fixed = TRUE))) {
    return("dragen")
  }

  # Default to CNVkit if cannot determine
  cat("WARNING: Could not automatically determine VCF source. Defaulting to CNVkit format.\n")
  return("cnvkit")
}

# Determine VCF source
if (opt$vcf_source == "auto") {
  vcf_source <- detect_vcf_source(opt$vcf)
  cat("Detected VCF source:", vcf_source, "\n")
} else {
  vcf_source <- opt$vcf_source
  cat("Using user-specified VCF source:", vcf_source, "\n")
}

# Read the TSV file with gene-level CNV data
cat("Reading TSV file:", opt$tsv, "\n")
gene_data <- fread(opt$tsv, header = TRUE)

# Read the VCF file for segment-level data
cat("Reading VCF file:", opt$vcf, "\n")
vcf_lines <- readLines(opt$vcf)

# Extract the sample ID from VCF
col_names_line <- grep("^#CHROM", vcf_lines, value = TRUE)
col_names <- strsplit(col_names_line, "\t")[[1]]
sample_id <- col_names[length(col_names)]
cat("Sample ID:", sample_id, "\n")

# Process segment data from VCF
data_lines <- vcf_lines[!grepl("^#", vcf_lines)]
vcf_data <- data.table()

if (vcf_source == "cnvkit") {
  # Process CNVkit VCF
  for (line in data_lines) {
    fields <- strsplit(line, "\t")[[1]]

    # Skip incomplete lines
    if (length(fields) < 8) next

    # Extract fields
    chrom <- gsub("^chr", "", fields[1])  # Remove 'chr' prefix if present
    pos <- as.numeric(fields[2])
    info_str <- fields[8]

    # Parse INFO field
    info_parts <- strsplit(info_str, ";")[[1]]
    info <- list()

    for (part in info_parts) {
      if (grepl("=", part)) {
        kv <- strsplit(part, "=")[[1]]
        info[[kv[1]]] <- kv[2]
      } else {
        info[[part]] <- TRUE
      }
    }

    # Extract required values
    if (!all(c("END", "FOLD_CHANGE", "SVTYPE") %in% names(info))) next

    end <- as.numeric(info[["END"]])
    fold_change <- as.numeric(info[["FOLD_CHANGE"]])
    log2_ratio <- fold_change_to_log2(fold_change)
    svtype <- info[["SVTYPE"]]
    probes <- if ("PROBES" %in% names(info)) as.numeric(info[["PROBES"]]) else 0

    # Add to segment data
    vcf_data <- rbind(vcf_data, data.table(
      ID = opt$sample_id,
      chrom = chrom,
      loc.start = pos,
      loc.end = end,
      num.mark = probes,
      seg.mean = log2_ratio
    ))
  }
} else if (vcf_source == "dragen") {
  # Process Dragen VCF
  for (line in data_lines) {
    #print(line)
    fields <- strsplit(line, "\t")[[1]]

    # Skip incomplete lines
    if (length(fields) < 10) next

    # Extract fields
    chrom <- gsub("^chr", "", fields[1])  # Remove 'chr' prefix if present
    pos <- as.numeric(fields[2])
    info_str <- fields[8]
    format_str <- fields[9]
    sample_str <- fields[10]

    # Parse INFO field
    info_parts <- strsplit(info_str, ";")[[1]]
    info <- list()

    for (part in info_parts) {
      if (grepl("=", part)) {
        kv <- strsplit(part, "=")[[1]]
        info[[kv[1]]] <- kv[2]
      } else {
        info[[part]] <- TRUE
      }
    }

    # Parse FORMAT field
    format_keys <- strsplit(format_str, ":")[[1]]
    sample_values <- strsplit(sample_str, ":")[[1]]
    format_data <- setNames(as.list(sample_values), format_keys)

    # Extract required values
    if (!("END" %in% names(info))) next

    end <- as.numeric(info[["END"]])
    svtype <- if ("SVTYPE" %in% names(info)) info[["SVTYPE"]] else "REF"

    # Skip REF segments (optional, comment out if you want to include them)
    if (svtype == "REF") next

    # Extract SM (segment mean) value from FORMAT field if available, or use CN/CNF for fold change calculation
     if ("SM" %in% names(format_data)) {
      sm <- as.numeric(format_data[["SM"]])
      seg_mean <- fold_change_to_log2(sm)
    } else {
      # Skip if we can't determine the copy number
      next
    }

    # Get number of bins if available
    probes <- if ("BC" %in% names(format_data)) as.numeric(format_data[["BC"]]) else 0

    # Add to segment data
    vcf_data <- rbind(vcf_data, data.table(
      ID = opt$sample_id,
      chrom = chrom,
      loc.start = pos,
      loc.end = end,
      num.mark = probes,
      seg.mean = seg_mean
    ))
  }
}

# Create discrete CNA data from the TSV file
cat("Creating discrete CNA data...\n")

# Calculate discrete values based on fold change
gene_data$discrete_cna <- sapply(gene_data$fold_change, fold_change_to_discrete)
gene_data$log2_ratio <- sapply(gene_data$fold_change, fold_change_to_log2)

# Create cBioPortal format data for discrete CNA
discrete_cna <- data.table(
  Hugo_Symbol = gene_data$gene_symbol
)
discrete_cna[[sample_id]] <- gene_data$discrete_cna

# Create cBioPortal format data for continuous (log2) CNA
log2_cna <- data.table(
  Hugo_Symbol = gene_data$gene_symbol
)
log2_cna[[sample_id]] <- gene_data$log2_ratio

# Create DISCRETE_LONG format
cat("Creating DISCRETE_LONG format...\n")
discrete_long <- data.table(
  Hugo_Symbol = gene_data$gene_symbol,
  Entrez_Gene_Id = gene_data$entrez_ncbi_id,
  Sample_Id = opt$sample_id,
  Value = gene_data$discrete_cna
)

# If there's a duplicate Entrez_Gene_Id, take the first one.
discrete_long <- discrete_long[!duplicated(discrete_long, by = c("Sample_Id", "Entrez_Gene_Id"))]

# Write segment file
seg_file <- file.path(opt$output_dir, paste0(opt$sample_id,"_data_cna_hg38.seg"))
cat("Writing segment file:", seg_file, "\n")
write.table(vcf_data, seg_file, sep = "\t", quote = FALSE, row.names = FALSE)


# Write discrete CNA file
discrete_file <- file.path(opt$output_dir, paste0(opt$sample_id,"_data_cna.txt"))
cat("Writing discrete CNA file:", discrete_file, "\n")
write.table(discrete_cna, discrete_file, sep = "\t", quote = FALSE, row.names = FALSE)


# Write DISCRETE_LONG format CNA file
discrete_long_file <- file.path(opt$output_dir, paste0(opt$sample_id,"_data_cna_long.txt"))
cat("Writing DISCRETE_LONG CNA file:", discrete_long_file, "\n")
write.table(discrete_long, discrete_long_file, sep = "\t", quote = FALSE, row.names = FALSE)

log2_file <- file.path(opt$output_dir, paste0(opt$sample_id,"_data_log2CNA.txt"))
cat("Writing log2 CNA file:", log2_file, "\n")
write.table(log2_cna, log2_file, sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nConversion complete!\n")
cat("Output files:\n")
cat("- ", seg_file, "\n")
cat("- ", discrete_file, "\n")
cat("- ", discrete_long_file, "\n")
cat("- ", log2_file, "\n")
