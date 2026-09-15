#!/usr/bin/env Rscript
#
# GENERATE_TIMELINE_TREATMENT — emits TREATMENT timeline events (systemic
# therapy) from treatments.csv rows flagged as systemic therapy. Sibling of
# gen_timeline_surgery.R (both read and filter treatments.csv independently).
# One of five per-EVENT_TYPE scripts merged by merge_timeline.R into data_timeline.txt.

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--sample_registrations", type="character", default=NULL,
              help="sample registrations CSV [OPTIONAL]", metavar="FILE"),
  make_option("--treatments",           type="character", default=NULL,
              help="treatments CSV [OPTIONAL]",            metavar="FILE"),
  make_option("--systemic_therapies",   type="character", default=NULL,
              help="systemic therapies CSV [OPTIONAL]",    metavar="FILE"),
  make_option("--genomic_subjects",     type="character", default=NULL,
              help="TSV with subject_id and sample_id columns; filters output to genomic subjects [OPTIONAL]", metavar="FILE"),
  make_option("--treatment_intent_map", type="character", default=NULL,
              help="MOHCCN treatment intent mapping CSV [OPTIONAL]",       metavar="FILE"),
  make_option(c("-o", "--output"), type="character", default="data_timeline_treatment.txt",
              help="Output file path [default= %default]", metavar="FILE"),
  make_option("--common_lib", type="character", default="clinical_common.R",
              help="Path to clinical_common.R [default= %default]", metavar="FILE")
)

opt_parser <- OptionParser(
  option_list = option_list,
  usage = "Usage: %prog [options]",
  description = "Generate the TREATMENT slice of the cBioPortal timeline from ARGO treatments"
)
opt <- parse_args(opt_parser)
source(opt$common_lib)

for (f in Filter(Negate(is.null), opt[setdiff(names(opt), c("help", "output", "common_lib"))])) {
  if (is.character(f) && nchar(f) > 0 && !file.exists(f)) {
    stop(paste("Error: File does not exist:", f))
  }
}

cat("=== Parameters ===\n")
for (nm in c("sample_registrations", "treatments", "systemic_therapies",
             "genomic_subjects", "treatment_intent_map")) {
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

ti_map <- read_mohccn_map(opt$treatment_intent_map)

syst_out <- NULL
if (!is.null(opt$treatments)) {
  cat("Reading treatments...\n")
  treat <- read.csv(opt$treatments, header = TRUE, stringsAsFactors = FALSE)
  treat <- treat[!is.na(treat$is_primary_treatment) & treat$is_primary_treatment == "Yes", ]
  treat$start_day        <- sapply(treat$treatment_start_date, parse_day_interval)
  treat$stop_day         <- sapply(treat$treatment_end_date,   parse_day_interval)
  treat$treatment_type_c <- sapply(treat$treatment_type,       clean_json_array)
  treat <- filter_patients(treat)

  syst_treat <- treat[grepl("Systemic therapy", treat$treatment_type_c, fixed = TRUE), ]

  if (nrow(syst_treat) > 0) {
    tx_base <- data.frame(
      PATIENT_ID       = syst_treat$submitter_donor_id,
      START_DATE       = syst_treat$start_day,
      STOP_DATE        = syst_treat$stop_day,
      tx_id            = syst_treat$submitter_treatment_id,
      TREATMENT_INTENT = if ("treatment_intent" %in% names(syst_treat)) syst_treat$treatment_intent else NA_character_,
      stringsAsFactors = FALSE
    )

    if (!is.null(opt$systemic_therapies)) {
      cat("Reading systemic therapies...\n")
      syst <- read.csv(opt$systemic_therapies, header = TRUE, stringsAsFactors = FALSE)
      syst_sub <- syst[, intersect(c("submitter_treatment_id", "drug_name", "systemic_therapy_type"), names(syst)), drop = FALSE]
      merged <- merge(tx_base, syst_sub, by.x = "tx_id", by.y = "submitter_treatment_id", all.x = TRUE)
      syst_out <- data.frame(
        PATIENT_ID            = merged$PATIENT_ID,
        START_DATE            = merged$START_DATE,
        STOP_DATE             = merged$STOP_DATE,
        EVENT_TYPE            = "TREATMENT",
        TREATMENT_TYPE        = if ("systemic_therapy_type" %in% names(merged)) merged$systemic_therapy_type else NA_character_,
        AGENT                 = if ("drug_name"             %in% names(merged)) merged$drug_name             else NA_character_,
        TREATMENT_INTENT      = merged$TREATMENT_INTENT,
        TREATMENT_INTENT_CODE = apply_mohccn_map(merged$TREATMENT_INTENT, ti_map),
        stringsAsFactors = FALSE
      )
    } else {
      syst_out <- data.frame(
        PATIENT_ID            = tx_base$PATIENT_ID,
        START_DATE            = tx_base$START_DATE,
        STOP_DATE             = tx_base$STOP_DATE,
        EVENT_TYPE            = "TREATMENT",
        TREATMENT_TYPE        = NA_character_,
        AGENT                 = NA_character_,
        TREATMENT_INTENT      = tx_base$TREATMENT_INTENT,
        TREATMENT_INTENT_CODE = apply_mohccn_map(tx_base$TREATMENT_INTENT, ti_map),
        stringsAsFactors = FALSE
      )
    }
    syst_out$START_DATE[is.na(syst_out$START_DATE)] <- 0
  }
}

# Always create the output file, even when empty. The Nextflow module declares
# this output non-optional: relying on an optional output that is genuinely
# absent has proven unreliable on some Nextflow versions (see merge_timeline.R
# and CLAUDE.md).
if (!is.null(syst_out) && nrow(syst_out) > 0) {
  write.table(syst_out, opt$output, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  cat(sprintf("Wrote %d row(s) to %s\n", nrow(syst_out), opt$output))
} else {
  file.create(opt$output)
  cat("No TREATMENT timeline data to write.\n")
}
