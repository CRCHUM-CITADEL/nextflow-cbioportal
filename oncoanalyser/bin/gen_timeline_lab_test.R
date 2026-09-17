#!/usr/bin/env Rscript
#
# GENERATE_TIMELINE_LAB_TEST — pivots biomarker results into cBioPortal LAB_TEST
# timeline events. One of five per-EVENT_TYPE scripts (see gen_timeline_surgery.R,
# gen_timeline_treatment.R, gen_timeline_status.R, gen_timeline_specimen.R) whose
# outputs are unioned by merge_timeline.R into data_timeline.txt.
#
# Split out of gen_timeline.R section (e); fully independent of the other four —
# needs only --biomarkers plus the patient-eligibility inputs shared by all five.

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--sample_registrations", type="character", default=NULL,
              help="sample registrations CSV [OPTIONAL]", metavar="FILE"),
  make_option("--biomarkers",           type="character", default=NULL,
              help="biomarkers CSV [OPTIONAL]",            metavar="FILE"),
  make_option("--genomic_subjects",     type="character", default=NULL,
              help="TSV with subject_id and sample_id columns; filters output to genomic subjects [OPTIONAL]", metavar="FILE"),
  make_option(c("-o", "--output"), type="character", default="data_timeline_lab_test.txt",
              help="Output file path [default= %default]", metavar="FILE"),
  make_option("--common_lib", type="character", default="clinical_common.R",
              help="Path to clinical_common.R [default= %default]", metavar="FILE")
)

opt_parser <- OptionParser(
  option_list = option_list,
  usage = "Usage: %prog [options]",
  description = "Generate the LAB_TEST slice of the cBioPortal timeline from ARGO biomarkers"
)
opt <- parse_args(opt_parser)
source(opt$common_lib)

# File existence checks
for (f in Filter(Negate(is.null), opt[setdiff(names(opt), c("help", "output", "common_lib"))])) {
  if (is.character(f) && nchar(f) > 0 && !file.exists(f)) {
    stop(paste("Error: File does not exist:", f))
  }
}

cat("=== Parameters ===\n")
for (nm in c("sample_registrations", "biomarkers", "genomic_subjects")) {
  cat(sprintf("%-22s %s\n", paste0(nm, ":"), ifelse(is.null(opt[[nm]]), "NULL", opt[[nm]])))
}
cat("==================\n\n")

# ── Valid patient IDs from sample registrations (same eligibility rule as clin_format.R) ──
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

# ── (e) Lab test timeline (biomarkers pivoted to long format) ─────────────────
lab_out <- NULL
if (!is.null(opt$biomarkers)) {
  cat("Reading biomarkers...\n")
  biom <- read.csv(opt$biomarkers, header = TRUE, stringsAsFactors = FALSE)
  biom$start_day <- sapply(biom$test_date, parse_day_interval)
  biom <- biom[!is.na(biom$start_day), ]
  biom <- filter_patients(biom)

  if (nrow(biom) > 0) {
    # Columns to pivot: name -> corresponding _not_available column (or NULL)
    biomarker_defs <- list(
      ca125               = "ca125_not_available",
      cea                 = "cea_not_available",
      psa_level           = "psa_level_not_available",
      er_percent_positive = "er_percent_positive_not_available",
      pr_percent_positive = "pr_percent_positive_not_available",
      er_status           = NULL,
      pr_status           = NULL,
      her2_ihc_status     = NULL,
      her2_ish_status     = NULL,
      hpv_ihc_status      = NULL,
      hpv_pcr_status      = NULL,
      hpv_strain          = NULL
    )

    lab_rows <- lapply(seq_len(nrow(biom)), function(i) {
      row <- biom[i, , drop = FALSE]
      result_rows <- lapply(names(biomarker_defs), function(col) {
        if (!(col %in% names(row))) return(NULL)
        val <- row[[col]]
        if (is.na(val) || nchar(trimws(as.character(val))) == 0) return(NULL)
        not_avail_col <- biomarker_defs[[col]]
        if (!is.null(not_avail_col) && not_avail_col %in% names(row)) {
          if (!is.na(row[[not_avail_col]]) && row[[not_avail_col]] == "True") return(NULL)
        }
        data.frame(
          PATIENT_ID = row$submitter_donor_id,
          START_DATE = row$start_day,
          STOP_DATE  = "",
          EVENT_TYPE = "LAB_TEST",
          TEST       = toupper(col),
          RESULT     = as.character(val),
          stringsAsFactors = FALSE
        )
      })
      do.call(rbind, Filter(Negate(is.null), result_rows))
    })

    lab_out <- do.call(rbind, Filter(Negate(is.null), lab_rows))
    if (!is.null(lab_out)) lab_out$START_DATE[is.na(lab_out$START_DATE)] <- 0
  }
}

if (!is.null(lab_out) && nrow(lab_out) > 0) {
  write.table(lab_out, opt$output, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  cat(sprintf("Wrote %d row(s) to %s\n", nrow(lab_out), opt$output))
} else {
  cat("No LAB_TEST timeline data to write.\n")
}

