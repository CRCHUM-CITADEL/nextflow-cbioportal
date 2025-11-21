#!/usr/bin/env Rscript

# Load necessary libraries
library(dplyr)
library(stringr)
library(readr)
library(optparse)

# Set up command line argument options
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Input fusion candidates file path", metavar="FILE"),
  make_option(c("-o", "--output"), type="character", default="data_sv.txt",
              help="Output file path [default= %default]", metavar="FILE"),
  make_option(c("-s", "--sample"), type="character", default=NULL,
              help="Sample ID to use (if not provided, will extract from input filename)")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check if input file was provided
if (is.null(opt$input)) {
  stop("Input file must be specified. Use -i or --input option.")
}

# Function to convert RNA fusion data to cBioPortal SV format
convert_to_cbioportal <- function(input_file, output_file, sample_id = NULL) {
  # Read the RNA fusion file
  fusions <- read_tsv(input_file)

  # Get sample ID from filename if not provided
  if (is.null(sample_id)) {
    sample_id <- str_replace(basename(input_file), "\\.RNA\\.fusion_candidates\\.final\\.txt$", "-1RT")
  }

  # Initialize output dataframe
  cbioportal_data <- data.frame(
    Sample_Id = rep(sample_id, nrow(fusions)),
    SV_Status = rep("SOMATIC", nrow(fusions)),
    Site1_Hugo_Symbol = NA,
    Site1_Ensembl_Transcript_Id = NA,
    Site1_Entrez_Gene_Id = NA,
    Site1_Region_Number = NA,
    Site1_Region = NA,
    Site2_Hugo_Symbol = NA,
    Site2_Ensembl_Transcript_Id = NA,
    Site2_Entrez_Gene_Id = NA,
    Site2_Region_Number = NA,
    Site2_Region = NA,
    Site2_Effect_On_Frame = NA,
    NCBI_Build = rep("GRCh38", nrow(fusions)),
    Class = NA,
    DNA_Support = rep("No", nrow(fusions)),
    RNA_Support = rep("Yes", nrow(fusions)),
    Tumor_Variant_Count = NA,
    Connection_Type = rep("5to3", nrow(fusions)),
    Breakpoint_Type = rep("PRECISE", nrow(fusions)),
    Event_Info = NA,
    Annotation = NA,
    Site1_Chromosome = NA,
    Site1_Position = NA,
    Site2_Chromosome = NA,
    Site2_Position = NA,
    Tumor_Split_Read_Count = NA,
    Tumor_Paired_End_Read_Count = NA,
    stringsAsFactors = FALSE
  )

  # Process each fusion
  for (i in 1:nrow(fusions)) {
    # Handle column name with possible # prefix
    fusion_column <- if("#FusionGene" %in% names(fusions)) "#FusionGene" else names(fusions)[1]

    # Split fusion gene name
    genes <- str_split(fusions[[fusion_column]][i], "--", simplify = TRUE)
    gene1 <- genes[1]
    gene2 <- genes[2]

    # Extract chromosome and position information
    left_bp <- str_split(fusions$LeftBreakpoint[i], ":", simplify = TRUE)
    right_bp <- str_split(fusions$RightBreakpoint[i], ":", simplify = TRUE)

    # Remove "chr" prefix if present
    chr1 <- str_remove(left_bp[1], "^chr")
    pos1 <- left_bp[2]
    chr2 <- str_remove(right_bp[1], "^chr")
    pos2 <- right_bp[2]

    # Set SV class
    sv_class <- "SV"

    # Fill in data
    cbioportal_data$Site1_Hugo_Symbol[i] <- gene1
    cbioportal_data$Site2_Hugo_Symbol[i] <- gene2
    cbioportal_data$Class[i] <- sv_class
    cbioportal_data$Event_Info[i] <- paste0("RNA-seq Fusion: ", gene1, "-", gene2)
    cbioportal_data$Annotation[i] <- paste0(gene1, " - ", gene2, " fusion")
    cbioportal_data$Site1_Chromosome[i] <- chr1
    cbioportal_data$Site1_Position[i] <- pos1
    cbioportal_data$Site2_Chromosome[i] <- chr2
    cbioportal_data$Site2_Position[i] <- pos2
    cbioportal_data$Tumor_Split_Read_Count[i] <- fusions$NumSplitReads[i]
    cbioportal_data$Tumor_Paired_End_Read_Count[i] <- fusions$NumPairedReads[i]
  }

  # Write output to file
  #print(nrow(cbioportal_data))
  #print(head(cbioportal_data))
  #print(nrow(unique(cbioportal_data)))
  #print(colnames(cbioportal_data))
  write_tsv(cbioportal_data[!duplicated(cbioportal_data[, 20:24]),], output_file)

  # Print summary
  cat(sprintf("Processed %d fusions from %s\n", nrow(fusions), input_file))
  cat(sprintf("Output written to %s\n", output_file))

  return(cbioportal_data)
}

# Run the conversion with command line arguments
result <- convert_to_cbioportal(opt$input, opt$output, opt$sample)
