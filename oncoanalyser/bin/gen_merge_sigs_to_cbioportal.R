#!/usr/bin/env Rscript
# Merge per-sample SIGS cBioPortal files into a single wide-format matrix.
# Each input file has columns: ENTITY_STABLE_ID, NAME, DESCRIPTION, {sample_id}
# Output: data_sigs.txt with all samples as columns, signatures as rows.
# Missing signatures for a given sample are filled with 0.

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
})

option_list <- list(
    make_option(c("-i", "--input_files"), type = "character",
                help = "Comma-separated list of per-sample data_sigs files [required]"),
    make_option(c("-o", "--output_file"), type = "character", default = "data_sigs.txt",
                help = "Output file path [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input_files)) stop("--input_files is required")

input_files <- trimws(strsplit(opt$input_files, ",")[[1]])
cat("Merging", length(input_files), "signature file(s)\n")

# Read each per-sample file; track the sample column (4th column)
read_sigs_file <- function(path) {
    dt <- fread(path, sep = "\t", header = TRUE)
    expected <- c("ENTITY_STABLE_ID", "NAME", "DESCRIPTION")
    missing <- setdiff(expected, colnames(dt))
    if (length(missing) > 0) stop("Missing columns in ", path, ": ", paste(missing, collapse = ", "))
    dt
}

tables <- lapply(input_files, read_sigs_file)

# Merge all tables by ENTITY_STABLE_ID, NAME, DESCRIPTION
merged <- tables[[1]]
if (length(tables) > 1) {
    for (i in 2:length(tables)) {
        merged <- merge(merged, tables[[i]],
                        by = c("ENTITY_STABLE_ID", "NAME", "DESCRIPTION"),
                        all = TRUE)
        cat("  Merged file", i, "of", length(tables), "\n")
    }
}

# Fill missing signature contributions with 0
sample_cols <- setdiff(colnames(merged), c("ENTITY_STABLE_ID", "NAME", "DESCRIPTION"))
for (col in sample_cols) {
    merged[[col]][is.na(merged[[col]])] <- 0
}

# Sort sample columns alphabetically for deterministic output
sample_cols_sorted <- sort(sample_cols)
setcolorder(merged, c("ENTITY_STABLE_ID", "NAME", "DESCRIPTION", sample_cols_sorted))

# Sort rows by ENTITY_STABLE_ID
setorder(merged, ENTITY_STABLE_ID)

cat("Output:", nrow(merged), "signatures x", length(sample_cols_sorted), "samples\n")
write.table(merged, opt$output_file, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
cat("Written to", opt$output_file, "\n")
