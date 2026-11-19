#!/usr/bin/env Rscript
# MERGE_TIMELINE: union the per-EVENT_TYPE parts; column order is a contract with combine_cbioportal_outputs.py

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--inputs", type="character", default=NULL,
              help="Comma-separated list of per-event timeline TSVs (missing/absent ones are skipped) [REQUIRED]", metavar="FILE,FILE,..."),
  make_option(c("-o", "--output"), type="character", default="data_timeline.txt",
              help="Output file path [default= %default]", metavar="FILE")
)

opt_parser <- OptionParser(
  option_list = option_list,
  usage = "Usage: %prog --inputs a.txt,b.txt,... [options]",
  description = "Union the per-EVENT_TYPE cBioPortal timeline slices into one data_timeline.txt"
)
opt <- parse_args(opt_parser)

if (is.null(opt$inputs)) stop("Error: --inputs argument is required")

# Bind data frames with different columns; missing columns filled with NA.
rbind_fill <- function(df_list) {
  df_list <- Filter(function(d) !is.null(d) && nrow(d) > 0, df_list)
  if (length(df_list) == 0) return(NULL)
  all_cols <- unique(unlist(lapply(df_list, names)))
  aligned  <- lapply(df_list, function(d) {
    missing <- setdiff(all_cols, names(d))
    for (col in missing) d[[col]] <- NA
    d[, all_cols, drop = FALSE]
  })
  do.call(rbind, aligned)
}

input_files <- strsplit(opt$inputs, ",")[[1]]
parts <- lapply(input_files, function(f) {
    # Parts with no rows are not staged (optional output); skip any empty file anyway
  if (!file.exists(f) || file.info(f)$size == 0) return(NULL)
  read.csv(f, sep = "\t", header = TRUE, stringsAsFactors = FALSE, colClasses = "character")
})

combined <- rbind_fill(parts)
if (!is.null(combined) && nrow(combined) > 0) {
    # Stable sort by a fixed EVENT_TYPE order so the output does not depend on part arrival order
  event_order <- factor(combined$EVENT_TYPE,
                         levels = c("SURGERY", "TREATMENT", "STATUS", "SPECIMEN", "LAB_TEST"))
  combined <- combined[order(event_order), , drop = FALSE]

  # Put common columns first, then the rest alphabetically
  common <- c("PATIENT_ID", "START_DATE", "STOP_DATE", "EVENT_TYPE")
  rest   <- sort(setdiff(names(combined), common))
  combined <- combined[, c(common, rest), drop = FALSE]
  write.table(combined, opt$output, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  cat(sprintf("Wrote %d row(s) to %s\n", nrow(combined), opt$output))
} else {
    # No file written: the output is optional, so the group gets no data_timeline.txt
  cat("No timeline data to write.\n")
}
