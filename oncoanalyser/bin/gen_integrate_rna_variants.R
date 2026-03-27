#!/usr/bin/env Rscript

# Final fixed script with direct position+allele matching approach
# Based on debug findings showing 371 exact matches

library(data.table)
library(dplyr)
library(optparse)

# Set up command line arguments
option_list <- list(
  make_option(c("-d", "--dna"), type="character", help="DNA MAF file path"),
  make_option(c("-r", "--rna"), type="character", help="RNA VCF file path"),
  make_option(c("-o", "--output"), type="character", help="Output MAF file path"),
  make_option(c("--min_depth"), type="integer", default=10, help="Minimum RNA read depth to consider a variant expressed [default: %default]"),
  make_option(c("--min_vaf"), type="double", default=0.03, help="Minimum RNA VAF to consider a variant expressed [default: %default]")
)

opt_parser <- OptionParser(option_list=option_list,
                          usage="Usage: %prog [options] or %prog -d dna.maf -r rna.vcf -o output.maf --min_depth=10 --min_vaf=0.03",
                          description="Integrate RNA-seq variant information into DNA MAF file")
opt <- parse_args(opt_parser)

# Check for positional arguments for backward compatibility
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 3 && is.null(opt$dna) && is.null(opt$rna) && is.null(opt$output)) {
  opt$dna <- args[1]
  opt$rna <- args[2]
  opt$output <- args[3]

  if (length(args) >= 4) {
    opt$min_depth <- as.integer(args[4])
  }

  if (length(args) >= 5) {
    opt$min_vaf <- as.numeric(args[5])
  }
}

# Check that required arguments are provided
if (is.null(opt$dna) || is.null(opt$rna) || is.null(opt$output)) {
  stop("Missing required arguments. Use --help for usage information.")
}

# Set the file paths and thresholds
dna_maf_file <- opt$dna
rna_vcf_file <- opt$rna
output_maf_file <- opt$output
MIN_RNA_DEPTH <- opt$min_depth
MIN_RNA_VAF <- opt$min_vaf

message("Parameters:")
message("  DNA MAF file: ", dna_maf_file)
message("  RNA VCF file: ", rna_vcf_file)
message("  Output file: ", output_maf_file)
message("  Min RNA depth: ", MIN_RNA_DEPTH)
message("  Min RNA VAF: ", MIN_RNA_VAF)

# Read DNA MAF file
message("Reading DNA MAF file...")
dna_maf <- fread(dna_maf_file, sep="\t", header=TRUE, stringsAsFactors=FALSE,
                 quote="", fill=TRUE, skip="Hugo_Symbol", data.table=FALSE)

# Read RNA VCF file
message("Reading RNA VCF file...")
vcf_lines <- readLines(rna_vcf_file)
header_line_idx <- grep("^#CHROM", vcf_lines)
vcf_header <- vcf_lines[header_line_idx]
vcf_data <- vcf_lines[(header_line_idx+1):length(vcf_lines)]

# Parse VCF header and data
vcf_columns <- unlist(strsplit(vcf_header, "\t"))
vcf_data_table <- data.table::fread(
  text = paste(vcf_data, collapse = "\n"),
  sep = "\t",
  header = FALSE,
  col.names = gsub("^#", "", vcf_columns),
  stringsAsFactors = FALSE
)

message("Processing ", nrow(vcf_data_table), " RNA variants...")

# Function to extract FORMAT field values
extract_format_value <- function(format_def, format_val, field_name) {
  # Split both fields
  format_fields <- unlist(strsplit(format_def, ":"))
  sample_values <- unlist(strsplit(format_val, ":"))

  # Find the position of the desired field
  field_pos <- which(format_fields == field_name)

  # Return NA if field not found or if the sample data doesn't contain the field
  if (length(field_pos) == 0 || field_pos > length(sample_values)) {
    return(NA)
  }

  # Return the value
  return(sample_values[field_pos])
}

# Extract basic information from VCF
rna_variants <- vcf_data_table %>%
  rename(FORMAT_COL = ncol(vcf_data_table)) %>%
  mutate(
    # Prepare for matching with both chromosome formats
    chr_raw = CHROM,
    Chromosome = gsub("^chr", "", CHROM),
    Start_Position = as.numeric(POS),
    Reference_Allele = REF,
    Variant_Allele = ALT
  )

# Extract depth, VAF, and ref/alt counts
message("Extracting RNA values...")
rna_variants$t_depth_rna <- sapply(1:nrow(rna_variants), function(i) {
  dp_str <- extract_format_value(rna_variants$FORMAT[i], rna_variants$FORMAT_COL[i], "DP")
  return(as.numeric(dp_str))
})

rna_variants$t_vaf_rna <- sapply(1:nrow(rna_variants), function(i) {
  af_str <- extract_format_value(rna_variants$FORMAT[i], rna_variants$FORMAT_COL[i], "AF")
  return(as.numeric(af_str))
})

rna_variants$ad_str <- sapply(1:nrow(rna_variants), function(i) {
  ad_str <- extract_format_value(rna_variants$FORMAT[i], rna_variants$FORMAT_COL[i], "AD")
  return(ad_str)
})

# Process AD (allele depth) strings to get ref/alt counts
rna_variants <- rna_variants %>%
  mutate(
    ad_split = strsplit(ad_str, ","),
    t_ref_count_rna = sapply(ad_split, function(x) if(length(x) >= 1) as.numeric(x[1]) else NA),
    t_alt_count_rna = sapply(ad_split, function(x) if(length(x) >= 2) as.numeric(x[2]) else NA)
  ) %>%
  select(chr_raw, Chromosome, Start_Position, Reference_Allele, Variant_Allele,
         t_depth_rna, t_vaf_rna, t_ref_count_rna, t_alt_count_rna)

# Also normalize chromosome in DNA data
message("Preparing DNA data for matching...")
dna_maf <- dna_maf %>%
  mutate(
    chr_raw = Chromosome,
    Chromosome = gsub("^chr", "", Chromosome)
  )

# DIRECT POSITION+ALLELE MATCHING APPROACH
message("Matching variants by position and alleles...")

# Create matching table with key RNA data
rna_match <- rna_variants %>%
  select(Chromosome, Start_Position, Reference_Allele, Variant_Allele, t_depth_rna, t_vaf_rna, t_ref_count_rna, t_alt_count_rna)

# Match using normalized chromosome (without chr prefix)
message("Matching using normalized chromosomes...")
matched_variants <- merge(
  dna_maf,
  rna_match,
  by.x = c("Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"),
  by.y = c("Chromosome", "Start_Position", "Reference_Allele", "Variant_Allele"),
  all.x = TRUE
)

# Set the expression flag
message("Setting expression flags...")
matched_variants <- matched_variants %>%
  mutate(
    Flag_RNA_Expressed = case_when(
      !is.na(t_depth_rna) & !is.na(t_vaf_rna) &
        as.numeric(t_depth_rna) >= MIN_RNA_DEPTH &
        as.numeric(t_vaf_rna) >= MIN_RNA_VAF ~ 1,
      TRUE ~ NA_real_
    )
  )

# Add RNA info to Mutation_Status
message("Adding RNA information to Mutation_Status column...")
matched_variants$rna_info <- sapply(1:nrow(matched_variants), function(i) {
  depth <- matched_variants$t_depth_rna[i]
  vaf <- matched_variants$t_vaf_rna[i]
  expressed <- matched_variants$Flag_RNA_Expressed[i]
  if(is.na(expressed)) {
  	expressed=0
  }

  if (!is.na(depth)) {
    paste0(
      "RNA(d=", depth,
      ",v=", ifelse(is.na(vaf), "NA", round(as.numeric(vaf), 2)),
      ",filt=", ifelse(expressed == 1, "YES", "NO"), ")"
    )
  } else {
    "RNA(data=NO)"
  }
})

# Update Mutation_Status column
if ("Mutation_Status" %in% colnames(matched_variants)) {
  message("Found existing Mutation_Status column, appending RNA information")
  matched_variants$Mutation_Status <- ifelse(
    is.na(matched_variants$Mutation_Status) | matched_variants$Mutation_Status == "",
    matched_variants$rna_info,
    paste0(matched_variants$Mutation_Status, ";", matched_variants$rna_info)
  )
} else {
  matched_variants$Mutation_Status <- matched_variants$rna_info
}

# Remove temporary columns
matched_variants$rna_info <- NULL
matched_variants$chr_raw <- NULL

# Write the extended MAF file
message("Writing integrated MAF file...")
# Get the header lines from the original MAF file
maf_header_lines <- readLines(dna_maf_file, n = 2)

# Create output columns list (keep original plus new columns)
output_cols <- c(names(dna_maf), "t_depth_rna", "t_vaf_rna", "t_ref_count_rna", "t_alt_count_rna", "Flag_RNA_Expressed")
output_cols <- unique(output_cols[!output_cols %in% c("chr_raw")]) # Remove any duplicates and temp columns

# Add Mutation_Status if present and not already included
if ("Mutation_Status" %in% names(matched_variants) && !("Mutation_Status" %in% output_cols)) {
  output_cols <- c(output_cols, "Mutation_Status")
}

# Make sure all columns exist in the data frame
output_cols <- output_cols[output_cols %in% names(matched_variants)]

# Write the new header
new_header <- c(maf_header_lines[1], paste(output_cols, collapse = "\t"))
writeLines(new_header, output_maf_file)

# Append data
fwrite(
  matched_variants[, output_cols, drop = FALSE],
  output_maf_file,
  sep = "\t",
  append = TRUE,
  quote = FALSE,
  col.names = FALSE,
  row.names = FALSE
)

message("Done! Integrated MAF file saved to: ", output_maf_file)

# Summary statistics
message("\nSummary statistics:")
message("Total DNA variants: ", nrow(dna_maf))
message("Total RNA variants with depth: ", sum(!is.na(rna_variants$t_depth_rna)))
message("Total RNA variants with VAF: ", sum(!is.na(rna_variants$t_vaf_rna)))

expressed_variants <- sum(matched_variants$Flag_RNA_Expressed == 1, na.rm = TRUE)
no_rna_data_variants <- sum(is.na(matched_variants$Flag_RNA_Expressed))

message("DNA variants expressed in RNA: ", expressed_variants, " (",
        round(expressed_variants/nrow(dna_maf)*100, 1), "%)")
message("DNA variants with no RNA data: ", no_rna_data_variants, " (",
        round(no_rna_data_variants/nrow(dna_maf)*100, 1), "%)")
