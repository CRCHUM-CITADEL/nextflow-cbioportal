#!/usr/bin/env Rscript
# Merge per-sample SIGS counts cBioPortal files into a single wide-format matrix.
# Each input file has columns: ENTITY_STABLE_ID, NAME, DESCRIPTION, {sample_id}
# Output: data_mutational_signatures_counts_SBS.txt with all samples as columns, contexts as rows.
# Missing contexts for a given sample are filled with 0.

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
})

option_list <- list(
    make_option(c("-i", "--input_files"), type = "character",
                help = "Comma-separated list of per-sample counts files [required]"),
    make_option(c("-o", "--output_file"), type = "character",
                default = "data_mutational_signatures_counts_SBS.txt",
                help = "Output file path [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input_files)) stop("--input_files is required")

input_files <- trimws(strsplit(opt$input_files, ",")[[1]])
cat("Merging", length(input_files), "counts file(s)\n")

# Known metadata column names (non-sample columns)
KNOWN_META_COLS <- c("ENTITY_STABLE_ID", "NAME", "DESCRIPTION", "CATEGORY", "URL")

read_counts_file <- function(path) {
    dt <- fread(path, sep = "\t", header = TRUE)
    if (!"ENTITY_STABLE_ID" %in% colnames(dt)) stop("Missing ENTITY_STABLE_ID in ", path)
    # Coerce metadata columns to character — fread infers logical for header-only files
    meta_cols <- intersect(KNOWN_META_COLS, colnames(dt))
    for (col in meta_cols) dt[[col]] <- as.character(dt[[col]])
    dt
}

tables <- lapply(input_files, read_counts_file)

# Auto-detect merge keys: metadata columns present in the first file
meta_keys <- intersect(KNOWN_META_COLS, colnames(tables[[1]]))
cat("Merge keys:", paste(meta_keys, collapse = ", "), "\n")

merged <- tables[[1]]
if (length(tables) > 1) {
    for (i in 2:length(tables)) {
        merged <- merge(merged, tables[[i]],
                        by = meta_keys,
                        all = TRUE)
        cat("  Merged file", i, "of", length(tables), "\n")
    }
}

# Fill missing counts with 0
sample_cols <- setdiff(colnames(merged), meta_keys)
for (col in sample_cols) {
    merged[[col]][is.na(merged[[col]])] <- 0
}

# Sort sample columns alphabetically for deterministic output
sample_cols_sorted <- sort(sample_cols)
setcolorder(merged, c(meta_keys, sample_cols_sorted))

# Sort rows by ENTITY_STABLE_ID
setorder(merged, ENTITY_STABLE_ID)

cat("Output:", nrow(merged), "contexts x", length(sample_cols_sorted), "samples\n")
write.table(merged, opt$output_file, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
cat("Written to", opt$output_file, "\n")
