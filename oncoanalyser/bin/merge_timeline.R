#!/usr/bin/env Rscript
#
# MERGE_TIMELINE — unions the per-EVENT_TYPE timeline slices (from
# gen_timeline_{surgery,treatment,status,specimen,lab_test}.R) into one
# data_timeline.txt.
#
# Column ordering is a hard contract with bin/combine_cbioportal_outputs.py's
# "union-append" merge strategy for data_timeline.txt (see CLAUDE.md): the 4
# common columns first, then the rest alphabetically. Do not change this order
# without updating that script's expectations too.

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
  # Each per-event script always creates its output file, even when it has no
  # rows to write (a zero-byte file) — Nextflow outputs here are declared
  # non-optional because relying on an optional output that is genuinely
  # absent has proven unreliable on some Nextflow versions. So a missing file
  # would be a real error, but a zero-byte one is the normal "no data" case
  # and must not be passed to read.csv() (which errors on an empty file).
  if (!file.exists(f)) stop(paste("Error: File does not exist:", f))
  if (file.info(f)$size == 0) return(NULL)
  read.csv(f, sep = "\t", header = TRUE, stringsAsFactors = FALSE, colClasses = "character")
})

combined <- rbind_fill(parts)
if (!is.null(combined) && nrow(combined) > 0) {
  # Stable-sort by EVENT_TYPE in the fixed order the original monolithic
  # gen_timeline.R always processed sections (a)-(e) in, so the merged row
  # order — and therefore the file's content/md5 — is deterministic
  # regardless of which order the five part-files physically arrive in
  # (mix()/groupTuple() in the subworkflow makes no ordering guarantee).
  # R's default order() is a stable radix sort for this size of input, so
  # rows within one EVENT_TYPE keep their original relative order.
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
  # Always create the output file, even when empty. The Nextflow module
  # declares this output non-optional for the same reason as the per-event
  # scripts above.
  file.create(opt$output)
  cat("No timeline data to write.\n")
}
