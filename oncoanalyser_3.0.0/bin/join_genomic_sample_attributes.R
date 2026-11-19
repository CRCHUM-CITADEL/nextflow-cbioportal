#!/usr/bin/env Rscript
#
# Append genomic sample attributes (TMB/TML, HRD from GENOMIC_SAMPLE_ATTRIBUTES) to a
# cBioPortal data_clinical_sample.txt, matched on SAMPLE_ID. cBioPortal accepts one
# SAMPLE_ATTRIBUTES file per study, so they cannot ship as a second file.
# Both inputs: 4 '#' metadata rows, the attribute-ID row, then data rows.

suppressPackageStartupMessages(library(optparse))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--clinical",   type = "character", help = "data_clinical_sample.txt"),
  make_option("--attributes", type = "character", help = "Merged genomic sample attributes"),
  make_option("--output",     type = "character", help = "Joined data_clinical_sample.txt")
)))
for (o in c("clinical", "attributes", "output")) {
  if (is.null(opt[[o]])) stop(sprintf("--%s is required", o))
}

read_cbio <- function(path) {
  lines <- readLines(path)
  if (length(lines) < 5 || !all(startsWith(lines[1:4], "#"))) {
    stop(sprintf("%s is not a cBioPortal clinical file (4 '#' rows + header)", path))
  }
  split <- function(l) strsplit(l, "\t", fixed = TRUE)[[1]]
  header <- split(lines[5])
  meta <- lapply(lines[1:4], function(l) split(sub("^#", "", l)))
  rows <- lapply(lines[-(1:5)], function(l) {
    v <- strsplit(l, "\t", fixed = TRUE)[[1]]
    c(v, rep("", length(header) - length(v)))  # strsplit drops trailing empty fields
  })
  rows <- rows[nzchar(lines[-(1:5)])]
  list(header = header, meta = meta, rows = rows)
}

clin  <- read_cbio(opt$clinical)
attrs <- read_cbio(opt$attributes)

if (!"SAMPLE_ID" %in% clin$header)  stop("clinical file has no SAMPLE_ID column")
if (!"SAMPLE_ID" %in% attrs$header) stop("attributes file has no SAMPLE_ID column")

new_cols <- setdiff(attrs$header, c("PATIENT_ID", "SAMPLE_ID"))
clash <- intersect(new_cols, clin$header)
if (length(clash)) stop(sprintf("columns already in the clinical file: %s", paste(clash, collapse = ", ")))
idx <- match(new_cols, attrs$header)

by_sample <- list()
for (r in attrs$rows) by_sample[[r[match("SAMPLE_ID", attrs$header)]]] <- r[idx]

sid <- match("SAMPLE_ID", clin$header)
clin_ids <- vapply(clin$rows, function(r) r[sid], "")
unmatched <- setdiff(names(by_sample), clin_ids)
if (length(unmatched)) {
  message(sprintf("WARNING: %d genomic sample(s) have no clinical sample row and are dropped: %s",
                  length(unmatched), paste(unmatched, collapse = ", ")))
}

out <- c(
  vapply(seq_along(clin$meta), function(i) {
    paste0("#", paste(c(clin$meta[[i]], attrs$meta[[i]][idx]), collapse = "\t"))
  }, ""),
  paste(c(clin$header, new_cols), collapse = "\t"),
  vapply(clin$rows, function(r) {
    v <- by_sample[[r[sid]]]
    if (is.null(v)) v <- rep("NA", length(new_cols))
    paste(c(r, v), collapse = "\t")
  }, "")
)
writeLines(out, opt$output)
cat(sprintf("Joined %d attribute column(s) onto %d sample row(s)\n", length(new_cols), length(clin$rows)))
