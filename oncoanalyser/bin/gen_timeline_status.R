#!/usr/bin/env Rscript
#
# GENERATE_TIMELINE_STATUS — emits STATUS timeline events (disease status at
# each follow-up visit) from follow_ups.csv. Sibling of gen_timeline_specimen.R
# (both read follow_ups.csv independently — SPECIMEN needs the whole follow-up
# frame to find the nearest visit per specimen). One of five per-EVENT_TYPE
# scripts merged by merge_timeline.R into data_timeline.txt.

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--sample_registrations", type="character", default=NULL,
              help="sample registrations CSV [OPTIONAL]", metavar="FILE"),
  make_option("--follow_ups",           type="character", default=NULL,
              help="follow-ups CSV [OPTIONAL]",            metavar="FILE"),
  make_option("--genomic_subjects",     type="character", default=NULL,
              help="TSV with subject_id and sample_id columns; filters output to genomic subjects [OPTIONAL]", metavar="FILE"),
  make_option("--primary_site_map",     type="character", default=NULL,
              help="MOHCCN primary site mapping CSV (text -> ICD-O code) [OPTIONAL]",  metavar="FILE"),
  make_option(c("-o", "--output"), type="character", default="data_timeline_status.txt",
              help="Output file path [default= %default]", metavar="FILE"),
  make_option("--common_lib", type="character", default="clinical_common.R",
              help="Path to clinical_common.R [default= %default]", metavar="FILE")
)

opt_parser <- OptionParser(
  option_list = option_list,
  usage = "Usage: %prog [options]",
  description = "Generate the STATUS slice of the cBioPortal timeline from ARGO follow-ups"
)
opt <- parse_args(opt_parser)
source(opt$common_lib)

for (f in Filter(Negate(is.null), opt[setdiff(names(opt), c("help", "output", "common_lib"))])) {
  if (is.character(f) && nchar(f) > 0 && !file.exists(f)) {
    stop(paste("Error: File does not exist:", f))
  }
}

cat("=== Parameters ===\n")
for (nm in c("sample_registrations", "follow_ups", "genomic_subjects", "primary_site_map")) {
  cat(sprintf("%-22s %s\n", paste0(nm, ":"), ifelse(is.null(opt[[nm]]), "NULL", opt[[nm]])))
}
cat("==================\n\n")

valid_patients <- character(0)
if (!is.null(opt$sample_registrations)) {
  cat("Reading sample registrations...\n")
  reg <- read.csv(opt$sample_registrations, header = TRUE, stringsAsFactors = FALSE)
  sample_reg <- reg[is_cbio_sample(reg), ]
  valid_patients <- unique(sample_reg$submitter_donor_id)
}

if (!is.null(opt$genomic_subjects)) {
  cat("Reading genomic subjects for filtering...\n")
  genomic        <- read_genomic_subjects(opt$genomic_subjects)
  valid_patients <- genomic$subject_id
}

filter_patients <- function(df, id_col = "submitter_donor_id") {
  if (length(valid_patients) > 0) df[df[[id_col]] %in% valid_patients, ] else df
}

# ps_map is REVERSED (ICD-O code -> label): anatomic_site_progression_or_recurrence
# is already an ICD-O code.
ps_map_raw <- read_mohccn_map(opt$primary_site_map)   # label -> code
ps_map     <- reverse_mohccn_map(ps_map_raw)

fu_data <- NULL
if (!is.null(opt$follow_ups)) {
  cat("Reading follow-ups...\n")
  fu_raw           <- read.csv(opt$follow_ups, header = TRUE, stringsAsFactors = FALSE)
  fu_raw$start_day <- sapply(fu_raw$date_of_followup, parse_day_interval)
  fu_data          <- fu_raw[!is.na(fu_raw$start_day), ]
  fu_data          <- filter_patients(fu_data)
}

status_out <- NULL
if (!is.null(fu_data) && nrow(fu_data) > 0) {
  site_raw   <- if ("anatomic_site_progression_or_recurrence" %in% names(fu_data))
    clean_json_array(fu_data$anatomic_site_progression_or_recurrence) else NA_character_
  method_raw <- if ("method_of_progression_status" %in% names(fu_data))
    clean_json_array(fu_data$method_of_progression_status) else NA_character_

  status_out <- data.frame(
    PATIENT_ID            = fu_data$submitter_donor_id,
    START_DATE            = fu_data$start_day,
    STOP_DATE             = "",
    EVENT_TYPE            = "STATUS",
    STATUS                = if ("disease_status_at_followup" %in% names(fu_data)) fu_data$disease_status_at_followup else NA_character_,
    RELAPSE_TYPE          = if ("relapse_type" %in% names(fu_data)) fu_data$relapse_type else NA_character_,
    RELAPSE_DATE          = sapply(fu_data$date_of_relapse, parse_day_interval, USE.NAMES = FALSE),
    METHOD_OF_PROGRESSION = method_raw,
    SITE                  = site_raw,
    SITE_LABEL            = apply_ps_label(site_raw, ps_map),
    stringsAsFactors = FALSE
  )
  status_out$START_DATE[is.na(status_out$START_DATE)] <- 0
}

if (!is.null(status_out) && nrow(status_out) > 0) {
  write.table(status_out, opt$output, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  cat(sprintf("Wrote %d row(s) to %s\n", nrow(status_out), opt$output))
} else {
  cat("No STATUS timeline data to write.\n")
}
