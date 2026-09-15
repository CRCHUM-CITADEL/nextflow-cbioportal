#!/usr/bin/env Rscript
#
# GENERATE_TIMELINE_SURGERY — emits SURGERY timeline events from treatments.csv
# rows flagged as surgical. Sibling of gen_timeline_treatment.R (both read and
# filter treatments.csv independently — cheap, and keeps each script single-purpose).
# One of five per-EVENT_TYPE scripts merged by merge_timeline.R into data_timeline.txt.

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--sample_registrations", type="character", default=NULL,
              help="sample registrations CSV [OPTIONAL]", metavar="FILE"),
  make_option("--treatments",           type="character", default=NULL,
              help="treatments CSV [OPTIONAL]",            metavar="FILE"),
  make_option("--surgeries",            type="character", default=NULL,
              help="surgeries CSV [OPTIONAL]",             metavar="FILE"),
  make_option("--genomic_subjects",     type="character", default=NULL,
              help="TSV with subject_id and sample_id columns; filters output to genomic subjects [OPTIONAL]", metavar="FILE"),
  make_option("--primary_site_map",     type="character", default=NULL,
              help="MOHCCN primary site mapping CSV (text -> ICD-O code) [OPTIONAL]",  metavar="FILE"),
  make_option("--treatment_intent_map", type="character", default=NULL,
              help="MOHCCN treatment intent mapping CSV [OPTIONAL]",       metavar="FILE"),
  make_option(c("-o", "--output"), type="character", default="data_timeline_surgery.txt",
              help="Output file path [default= %default]", metavar="FILE"),
  make_option("--common_lib", type="character", default="clinical_common.R",
              help="Path to clinical_common.R [default= %default]", metavar="FILE")
)

opt_parser <- OptionParser(
  option_list = option_list,
  usage = "Usage: %prog [options]",
  description = "Generate the SURGERY slice of the cBioPortal timeline from ARGO treatments"
)
opt <- parse_args(opt_parser)
source(opt$common_lib)

for (f in Filter(Negate(is.null), opt[setdiff(names(opt), c("help", "output", "common_lib"))])) {
  if (is.character(f) && nchar(f) > 0 && !file.exists(f)) {
    stop(paste("Error: File does not exist:", f))
  }
}

cat("=== Parameters ===\n")
for (nm in c("sample_registrations", "treatments", "surgeries",
             "genomic_subjects", "primary_site_map", "treatment_intent_map")) {
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

# ps_map is REVERSED (ICD-O code -> label): surgery_site is already an ICD-O code.
ps_map_raw <- read_mohccn_map(opt$primary_site_map)   # label -> code
ps_map     <- reverse_mohccn_map(ps_map_raw)
ti_map     <- read_mohccn_map(opt$treatment_intent_map)

surg_out <- NULL
if (!is.null(opt$treatments)) {
  cat("Reading treatments...\n")
  treat <- read.csv(opt$treatments, header = TRUE, stringsAsFactors = FALSE)
  treat <- treat[!is.na(treat$is_primary_treatment) & treat$is_primary_treatment == "Yes", ]
  treat$start_day        <- sapply(treat$treatment_start_date, parse_day_interval)
  treat$stop_day         <- sapply(treat$treatment_end_date,   parse_day_interval)
  treat$treatment_type_c <- sapply(treat$treatment_type,       clean_json_array)
  treat <- filter_patients(treat)

  surg_treat <- treat[grepl("Surgery", treat$treatment_type_c, fixed = TRUE), ]

  if (nrow(surg_treat) > 0) {
    surg_df <- data.frame(
      PATIENT_ID       = surg_treat$submitter_donor_id,
      START_DATE       = surg_treat$start_day,
      STOP_DATE        = surg_treat$stop_day,
      EVENT_TYPE       = "SURGERY",
      tx_id            = surg_treat$submitter_treatment_id,
      SUBTYPE          = NA_character_,
      SITE             = NA_character_,
      TREATMENT_INTENT = if ("treatment_intent" %in% names(surg_treat)) surg_treat$treatment_intent else NA_character_,
      stringsAsFactors = FALSE
    )

    if (!is.null(opt$surgeries)) {
      cat("Reading surgeries...\n")
      surg <- read.csv(opt$surgeries, header = TRUE, stringsAsFactors = FALSE)
      surg_sub <- surg[, intersect(c("submitter_treatment_id", "surgery_type", "surgery_site"), names(surg)), drop = FALSE]
      surg_df <- merge(surg_df, surg_sub, by.x = "tx_id", by.y = "submitter_treatment_id", all.x = TRUE)
      if ("surgery_type" %in% names(surg_df)) surg_df$SUBTYPE <- surg_df$surgery_type
      if ("surgery_site" %in% names(surg_df)) surg_df$SITE   <- surg_df$surgery_site
    }

    surg_df$SITE_LABEL           <- apply_ps_label(surg_df$SITE, ps_map)
    surg_df$TREATMENT_INTENT_CODE <- apply_mohccn_map(surg_df$TREATMENT_INTENT, ti_map)
    surg_out <- surg_df[, c("PATIENT_ID", "START_DATE", "STOP_DATE", "EVENT_TYPE", "SUBTYPE",
                            "SITE", "SITE_LABEL", "TREATMENT_INTENT", "TREATMENT_INTENT_CODE")]
    surg_out$START_DATE[is.na(surg_out$START_DATE)] <- 0
  }
}

# Always create the output file, even when empty. The Nextflow module declares
# this output non-optional: relying on an optional output that is genuinely
# absent has proven unreliable on some Nextflow versions (see merge_timeline.R
# and CLAUDE.md).
if (!is.null(surg_out) && nrow(surg_out) > 0) {
  write.table(surg_out, opt$output, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  cat(sprintf("Wrote %d row(s) to %s\n", nrow(surg_out), opt$output))
} else {
  file.create(opt$output)
  cat("No SURGERY timeline data to write.\n")
}
