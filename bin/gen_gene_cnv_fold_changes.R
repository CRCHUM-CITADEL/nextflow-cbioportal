#!/usr/bin/env Rscript

# Dragen CNV to CNVkit Format Converter
# Maps Ensembl IDs to gene symbols and formats Dragen CNV output to match CNVkit format

library(data.table)
library(optparse)

# Parse command line arguments
option_list <- list(
  make_option(c("-v", "--vcf"), type="character", default=NULL,
              help="Input Dragen CNV VCF file"),
  make_option(c("-a", "--annotation"), type="character", default=NULL,
              help="Input gene annotation TSV file (biomart format)"),
  make_option(c("-o", "--output"), type="character", default="cnv_fold_changes_per_gene_dragen.tsv",
              help="Output file name [default: %default]")
)

opt_parser <- OptionParser(option_list=option_list,
                          description="Map Ensembl IDs in Dragen CNV VCF to gene symbols")
opt <- parse_args(opt_parser)

# Check required arguments
if (is.null(opt$vcf) || is.null(opt$annotation)) {
  cat("ERROR: Both VCF and annotation files are required\n")
  print_help(opt_parser)
  quit(status = 1)
}

# Read the annotation file
cat(sprintf("Reading annotation file: %s\n", opt$annotation))
annotations <- fread(opt$annotation, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Clean up the chromosome format in annotations and ensure numeric columns are numeric
annotations$chr <- gsub("^chr", "", annotations$chr)
annotations$start <- as.numeric(annotations$start)
annotations$stop <- as.numeric(annotations$end)

# Define canonical chromosomes
canonical_chr <- c(as.character(1:22), "X", "Y", "MT", "M")

# Filter to keep only genes on canonical chromosomes
annotations <- annotations[annotations$chr %in% canonical_chr, ]

# Create an ordered factor for chromosome sorting
chr_order <- c(as.character(1:22), "X", "Y", "MT", "M")
annotations$chr_factor <- factor(annotations$chr, levels = chr_order, ordered = TRUE)

# Output statistics about annotations
cat("\nAnnotation Statistics:\n")
cat(sprintf("Total annotations: %d\n", nrow(annotations)))
cat(sprintf("Unique Ensembl IDs: %d\n", length(unique(annotations$ensembl_id))))
cat(sprintf("Unique Entrez IDs: %d\n", length(unique(annotations$entrez_ncbi_id))))
cat(sprintf("Unique gene symbols: %d\n", length(unique(annotations$gene_symbol))))

# Read the VCF file
cat(sprintf("\nReading VCF file: %s\n", opt$vcf))
vcf_lines <- readLines(opt$vcf)
header_lines <- vcf_lines[grep("^#", vcf_lines)]
data_lines <- vcf_lines[!grepl("^#", vcf_lines)]

# Check if we have data lines
if (length(data_lines) == 0) {
  cat("ERROR: No data lines found in VCF file\n")
  quit(status = 1)
}

# Process VCF data lines
vcf_data <- fread(text = data_lines, header = FALSE, sep = "\t")
colnames(vcf_data) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE")

cat("Extracting CNV regions...\n")
vcf_data$chr <- gsub("^chr", "", vcf_data$CHROM)
vcf_data$start <- as.numeric(vcf_data$POS)
vcf_data$end <- sapply(vcf_data$INFO, function(info) {
  end_match <- regmatches(info, regexpr("END=[0-9]+", info))
  as.numeric(gsub("END=", "", end_match))
})

# Determine if segment is a DEL or DUP based on ID field
vcf_data$svtype <- sapply(vcf_data$ID, function(id) {
  if (grepl("LOSS", id)) {
    return("DEL")
  } else if (grepl("GAIN", id)) {
    return("DUP")
  } else if (grepl("REF", id)) {
    return(NA)  # Skip reference segments
  } else {
    return(NA)  # Skip unknown types
  }
})

# Filter out rows with NA svtype (reference regions)
vcf_data <- vcf_data[!is.na(vcf_data$svtype), ]

# Extract CN (copy number) from SAMPLE field
vcf_data$copy_number <- sapply(1:nrow(vcf_data), function(i) {
  format_fields <- unlist(strsplit(vcf_data$FORMAT[i], ":"))
  sample_values <- unlist(strsplit(vcf_data$SAMPLE[i], ":"))
  cn_index <- which(format_fields == "CNF")
  if (length(cn_index) > 0 && cn_index <= length(sample_values)) {
    return(as.numeric(sample_values[cn_index]))
  } else {
    return(NA)
  }
})

# Calculate fold change from copy number
#vcf_data$fold_change <- vcf_data$copy_number / 2

vcf_data$fold_change <- sapply(1:nrow(vcf_data), function(i) {
  format_fields <- unlist(strsplit(vcf_data$FORMAT[i], ":"))
  sample_values <- unlist(strsplit(vcf_data$SAMPLE[i], ":"))
  foldchange_index <- which(format_fields == "SM")
  return(as.numeric(sample_values[foldchange_index]))
})

# Map each CNV to overlapping genes
cat("Mapping CNVs to genes...\n")
result <- data.table()
print(str(annotations))
print("Sample CNV:")
print(head(vcf_data, 1))
print("Sample annotation:")
print(head(annotations, 1))
result <- data.table()
print(str(annotations))


for (i in 1:nrow(vcf_data)) {
  cnv <- vcf_data[i, ]


  # Print progress every 10 CNVs
  if (i %% 10 == 0) {
    cat(sprintf("Processing CNV %d of %d...\n", i, nrow(vcf_data)))
  }
  
  # Debug: print first CNV details
  if (i == 1) {
    cat("First CNV details:\n")
    cat("  chr:", cnv$chr, "\n")
    cat("  start:", cnv$start, "\n")
    cat("  end:", cnv$end, "\n")
    cat("Matching genes on same chr:", nrow(annotations[annotations$chr == cnv$chr,]), "\n")
  }

  # Find genes overlapping with this CNV
  overlapping_genes <- annotations[
    annotations$chr == cnv$chr &
    annotations$start <= cnv$end &
    annotations$stop >= cnv$start,
  ]

  if (nrow(overlapping_genes) > 0) {
    # Add to results with gene information
    for (j in 1:nrow(overlapping_genes)) {
      gene <- overlapping_genes[j, ]

      cnv_result <- data.table(
        chr = cnv$chr,
        start = cnv$start,
        end = cnv$end,
        svtype = cnv$svtype,
        fold_change = cnv$fold_change,
        ensembl_id = gene$ensembl_id,
	entrez_ncbi_id = gene$entrez_ncbi_id,
        gene_symbol = gene$gene_symbol,
        gene_chr = gene$chr,
        gene_start = gene$start,
        gene_stop = gene$stop,
        gene_strand = gene$strand,
        gene_description = gene$description
      )

      # Add gene_biotype if present in annotations
      if(cnv_result$gene_symbol!="") {
      result <- rbind(result, cnv_result)
      }
    }
  }
}

# After collecting all results, find the most significant alteration for each gene
cat("\nPost-processing results to keep only the most significant alteration per gene...\n")

# Create a copy of the original results
all_results <- copy(result)

# Empty the result table
result <- data.table()

# Group by gene symbol and find the most significant alteration
gene_symbols <- unique(all_results$gene_symbol)
cat(sprintf("Processing %d unique genes to find most significant alterations...\n", length(gene_symbols)))

for (symbol in gene_symbols) {
  # Get all rows for this gene
  gene_rows <- all_results[all_results$gene_symbol == symbol, ]

  # If only one entry, keep it
  if (nrow(gene_rows) == 1) {
    result <- rbind(result, gene_rows)
    next
  }

  # Find the most significant alteration based on fold change
  # For gains (fold_change > 1), higher fold_change is more significant
  # For losses (fold_change < 1), lower fold_change is more significant
  gains <- gene_rows[gene_rows$fold_change > 1, ]
  losses <- gene_rows[gene_rows$fold_change < 1, ]

  # If we have both gains and losses, take the most extreme one
  if (nrow(gains) > 0 && nrow(losses) > 0) {
    # Find max gain and min loss
    max_gain <- gains[which.max(gains$fold_change), ]
    min_loss <- losses[which.min(losses$fold_change), ]

    # Compare deviations from normal
    gain_deviation <- max_gain$fold_change - 1
    loss_deviation <- 1 - min_loss$fold_change

    if (loss_deviation > gain_deviation) {
      result <- rbind(result, min_loss)
    } else {
      result <- rbind(result, max_gain)
    }
  }
  # If we only have gains, take the highest
  else if (nrow(gains) > 0) {
    max_gain <- gains[which.max(gains$fold_change), ]
    result <- rbind(result, max_gain)
  }
  # If we only have losses, take the lowest
  else if (nrow(losses) > 0) {
    min_loss <- losses[which.min(losses$fold_change), ]
    result <- rbind(result, min_loss)
  }
}

# Ensure all required columns exist and are in correct order
output_columns <- c("chr", "start", "end", "svtype", "fold_change",
                   "ensembl_id", "gene_symbol", "gene_chr", "gene_start", "gene_stop",
                   "gene_strand", "gene_description", "entrez_ncbi_id")

# Add missing columns
for (col in output_columns) {
  if (!(col %in% colnames(result))) {
    cat(sprintf("Adding missing column: %s\n", col))
    result[[col]] <- NA
  }
}

# Reorder columns to match the CNVkit output
result <- result[, ..output_columns]

# Output statistics
cat("\nMapping Results:\n")
cat(sprintf("Total CNVs in VCF: %d\n", nrow(vcf_data)))
cat(sprintf("Original gene-CNV mappings: %d\n", nrow(all_results)))
cat(sprintf("Final unique gene mappings: %d\n", nrow(result)))

# Write results to output file
cat(sprintf("\nWriting results to: %s\n", opt$output))
fwrite(result, opt$output, sep="\t", quote=FALSE)

cat("Done!\n")
