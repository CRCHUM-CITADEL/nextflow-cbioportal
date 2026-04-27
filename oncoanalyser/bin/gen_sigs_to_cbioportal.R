#!/usr/bin/env Rscript
# Convert HMFTOOLS SIGS allocation output to per-sample cBioPortal GENERIC_ASSAY format.
# Input: {subject}-T.sig.allocation.tsv (columns: signature, allocation, percent)
# Output: {sample}.data_sigs.txt (ENTITY_STABLE_ID, NAME, DESCRIPTION, {sample})

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
})

option_list <- list(
    make_option(c("-i", "--input"),  type = "character", help = "Path to .sig.allocation.tsv file [required]"),
    make_option(c("-s", "--sample"), type = "character", help = "Sample ID for output column header [required]"),
    make_option(c("-o", "--output"), type = "character", help = "Output file path [required]")
)

opt <- parse_args(OptionParser(option_list = option_list))

for (arg in c("input", "sample", "output")) {
    if (is.null(opt[[arg]])) stop(paste0("Missing required argument: --", arg))
}

cat("Reading SIGS allocation file:", opt$input, "\n")
sigs <- fread(opt$input, sep = "\t", header = TRUE)

# Validate expected columns
expected_cols <- c("signature", "allocation", "percent")
missing_cols <- setdiff(expected_cols, colnames(sigs))
if (length(missing_cols) > 0) {
    stop("Missing columns in input: ", paste(missing_cols, collapse = ", "))
}

cat("Found", nrow(sigs), "allocated signatures\n")

# Build cBioPortal GENERIC_ASSAY format (one row per signature)
out <- data.table(
    ENTITY_STABLE_ID = sigs$signature,
    NAME             = sigs$signature,
    DESCRIPTION      = paste0("Mutational signature ", sigs$signature)
)
out[[opt$sample]] <- sigs$percent

write.table(out, opt$output, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
cat("Written", nrow(out), "rows to", opt$output, "\n")
