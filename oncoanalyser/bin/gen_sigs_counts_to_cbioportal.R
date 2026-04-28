#!/usr/bin/env Rscript
# Convert HMFTOOLS SIGS SNV counts to per-sample cBioPortal GENERIC_ASSAY format.
# Input: {subject}-T.sig.snv_counts.csv (columns: BucketName, {sample_id})
# Output: {sample}.data_mutational_signatures_counts_SBS.txt (ENTITY_STABLE_ID, NAME, DESCRIPTION, {sample})
# Context names are transformed from C>A_ACA to A[C>A]A format.

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
})

option_list <- list(
    make_option(c("-i", "--input"),  type = "character", help = "Path to .sig.snv_counts.csv file [required]"),
    make_option(c("-s", "--sample"), type = "character", help = "Sample ID for output column header [required]"),
    make_option(c("-o", "--output"), type = "character", help = "Output file path [required]")
)

opt <- parse_args(OptionParser(option_list = option_list))

for (arg in c("input", "sample", "output")) {
    if (is.null(opt[[arg]])) stop(paste0("Missing required argument: --", arg))
}

# Transform C>A_ACA → A[C>A]A
transform_context <- function(x) {
    parts <- strsplit(x, "_")[[1]]   # e.g. c("C>A", "ACA")
    if (length(parts) != 2) stop(paste0("Unexpected BucketName format: ", x))
    mut <- parts[1]                   # "C>A"
    ctx <- parts[2]                   # "ACA"
    if (nchar(ctx) != 3) stop(paste0("Expected 3-character context, got: ", ctx))
    paste0(substr(ctx, 1, 1), "[", mut, "]", substr(ctx, 3, 3))
}

cat("Reading SNV counts file:", opt$input, "\n")
counts <- fread(opt$input, sep = ",", header = TRUE)

if (ncol(counts) < 2) stop("Expected at least 2 columns in SNV counts file")
if (colnames(counts)[1] != "BucketName") stop("Expected first column to be 'BucketName'")

cat("Found", nrow(counts), "trinucleotide contexts\n")

entity_ids <- sapply(counts$BucketName, transform_context, USE.NAMES = FALSE)

out <- data.table(
    ENTITY_STABLE_ID = entity_ids,
    NAME             = entity_ids,
    DESCRIPTION      = paste0("Trinucleotide context: ", entity_ids)
)
out[[opt$sample]] <- counts[[2]]

write.table(out, opt$output, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
cat("Written", nrow(out), "rows to", opt$output, "\n")
