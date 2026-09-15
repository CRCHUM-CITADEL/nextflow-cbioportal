#!/usr/bin/env Rscript
#
# GENERATE_TIMELINE_SPECIMEN — emits SPECIMEN timeline events (collection date,
# site, tissue source, and SAMPLE_TYPE derived from the nearest follow-up visit)
# from specimens.csv. Sibling of gen_timeline_status.R (both read follow_ups.csv;
# here it is only used to find the nearest visit per specimen, not to emit its
# own events). One of five per-EVENT_TYPE scripts merged by merge_timeline.R
# into data_timeline.txt.

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--sample_registrations",       type="character", default=NULL,
              help="sample registrations CSV [OPTIONAL]", metavar="FILE"),
  make_option("--specimens",                  type="character", default=NULL,
              help="specimens CSV [OPTIONAL]",             metavar="FILE"),
  make_option("--follow_ups",                 type="character", default=NULL,
              help="follow-ups CSV [OPTIONAL]",            metavar="FILE"),
  make_option("--genomic_subjects",           type="character", default=NULL,
              help="TSV with subject_id and sample_id columns; filters output to genomic subjects [OPTIONAL]", metavar="FILE"),
  make_option("--primary_site_map",           type="character", default=NULL,
              help="MOHCCN primary site mapping CSV (text -> ICD-O code) [OPTIONAL]",  metavar="FILE"),
  make_option("--specimen_tissue_source_map", type="character", default=NULL,
              help="MOHCCN specimen tissue source mapping CSV [OPTIONAL]", metavar="FILE"),
  make_option(c("-o", "--output"), type="character", default="data_timeline_specimen.txt",
              help="Output file path [default= %default]", metavar="FILE"),
  make_option("--common_lib", type="character", default="clinical_common.R",
              help="Path to clinical_common.R [default= %default]", metavar="FILE")
)

opt_parser <- OptionParser(
  option_list = option_list,
  usage = "Usage: %prog [options]",
  description = "Generate the SPECIMEN slice of the cBioPortal timeline from ARGO specimens"
)
opt <- parse_args(opt_parser)
source(opt$common_lib)

for (f in Filter(Negate(is.null), opt[setdiff(names(opt), c("help", "output", "common_lib"))])) {
  if (is.character(f) && nchar(f) > 0 && !file.exists(f)) {
    stop(paste("Error: File does not exist:", f))
  }
}

cat("=== Parameters ===\n")
for (nm in c("sample_registrations", "specimens", "follow_ups",
             "genomic_subjects", "primary_site_map", "specimen_tissue_source_map")) {
  cat(sprintf("%-22s %s\n", paste0(nm, ":"), ifelse(is.null(opt[[nm]]), "NULL", opt[[nm]])))
}
cat("==================\n\n")

# ── Valid patient IDs + specimen maps from sample registrations ────────────────
valid_patients         <- character(0)
spec_type_map          <- character(0)   # submitter_specimen_id -> specimen_type
spec_tissue_source_map <- character(0)   # submitter_specimen_id -> specimen_tissue_source

if (!is.null(opt$sample_registrations)) {
  cat("Reading sample registrations...\n")
  reg <- read.csv(opt$sample_registrations, header = TRUE, stringsAsFactors = FALSE)
  sample_reg <- reg[is_cbio_sample(reg), ]
  valid_patients <- unique(sample_reg$submitter_donor_id)
  # Deduplicate by specimen_id — multiple sequencing types (DNA, RNA) per specimen produce
  # duplicate rows; take the first entry so the named-vector lookup stays unambiguous.
  # Tumour rows sort first so a normal sharing a specimen ID cannot supply the
  # SPECIMEN event's specimen_type / specimen_tissue_source.
  reg_ord   <- reg[order(reg$tumour_normal_designation != "Tumour"), ]
  reg_dedup <- reg_ord[!duplicated(reg_ord$submitter_specimen_id), ]
  if ("specimen_type" %in% names(reg_dedup)) {
    spec_type_map <- setNames(reg_dedup$specimen_type, reg_dedup$submitter_specimen_id)
  }
  if ("specimen_tissue_source" %in% names(reg_dedup)) {
    spec_tissue_source_map <- setNames(reg_dedup$specimen_tissue_source, reg_dedup$submitter_specimen_id)
  }
}

if (!is.null(opt$genomic_subjects)) {
  cat("Reading genomic subjects for filtering...\n")
  genomic        <- read_genomic_subjects(opt$genomic_subjects)
  valid_patients <- genomic$subject_id
}

filter_patients <- function(df, id_col = "submitter_donor_id") {
  if (length(valid_patients) > 0) df[df[[id_col]] %in% valid_patients, ] else df
}

# ps_map is REVERSED (ICD-O code -> label): specimen_anatomic_location is
# already an ICD-O code.
ps_map_raw <- read_mohccn_map(opt$primary_site_map)   # label -> code
ps_map     <- reverse_mohccn_map(ps_map_raw)
sts_map    <- read_mohccn_map(opt$specimen_tissue_source_map)

# Derive cBioPortal SAMPLE_TYPE from the nearest follow-up disease status.
get_sample_type <- function(donor_id, specimen_day, fu_df) {
  if (is.null(fu_df) || nrow(fu_df) == 0) return(NA_character_)
  patient_fu <- fu_df[fu_df$submitter_donor_id == donor_id & !is.na(fu_df$start_day), ]
  if (nrow(patient_fu) == 0) return(NA_character_)
  nearest <- patient_fu[which.min(abs(patient_fu$start_day - specimen_day)), ]
  status  <- nearest$disease_status_at_followup
  if (is.na(status) || nchar(trimws(as.character(status))) == 0) return(NA_character_)
  if (grepl("metastatic|metastasis", status, ignore.case = TRUE)) return("Metastasis")
  if (grepl("relapse|recurrence|recurred|progression|progressed", status, ignore.case = TRUE)) return("Recurrence")
  return("Primary")
}

fu_data <- NULL
if (!is.null(opt$follow_ups)) {
  cat("Reading follow-ups...\n")
  fu_raw           <- read.csv(opt$follow_ups, header = TRUE, stringsAsFactors = FALSE)
  fu_raw$start_day <- sapply(fu_raw$date_of_followup, parse_day_interval)
  fu_data          <- fu_raw[!is.na(fu_raw$start_day), ]
  fu_data          <- filter_patients(fu_data)
}

spec_out <- NULL
if (!is.null(opt$specimens)) {
  cat("Reading specimens...\n")
  spec <- read.csv(opt$specimens, header = TRUE, stringsAsFactors = FALSE)
  spec$start_day <- sapply(spec$specimen_collection_date, parse_day_interval)
  spec <- spec[!is.na(spec$start_day), ]
  spec <- filter_patients(spec)

  if (nrow(spec) > 0) {
    spec_type_vals <- if (length(spec_type_map) > 0) {
      unname(spec_type_map[spec$submitter_specimen_id])
    } else {
      rep(NA_character_, nrow(spec))
    }
    spec_tissue_source_vals <- if (length(spec_tissue_source_map) > 0) {
      unname(spec_tissue_source_map[spec$submitter_specimen_id])
    } else {
      rep(NA_character_, nrow(spec))
    }

    spec_site_vals <- if ("specimen_anatomic_location" %in% names(spec)) spec$specimen_anatomic_location else NA_character_

    sample_type_vals <- mapply(get_sample_type,
      spec$submitter_donor_id, spec$start_day,
      MoreArgs = list(fu_df = fu_data),
      USE.NAMES = FALSE)

    spec_out <- data.frame(
      PATIENT_ID                  = spec$submitter_donor_id,
      START_DATE                  = spec$start_day,
      STOP_DATE                   = "",
      EVENT_TYPE                  = "SPECIMEN",
      SPECIMEN_SITE               = spec_site_vals,
      SPECIMEN_SITE_LABEL         = apply_ps_label(spec_site_vals, ps_map),
      SPECIMEN_TYPE               = spec_type_vals,
      SAMPLE_TYPE                 = sample_type_vals,
      SPECIMEN_TISSUE_SOURCE      = spec_tissue_source_vals,
      SPECIMEN_TISSUE_SOURCE_CODE = apply_mohccn_map(spec_tissue_source_vals, sts_map),
      stringsAsFactors = FALSE
    )
    spec_out$START_DATE[is.na(spec_out$START_DATE)] <- 0
  }
}

# Always create the output file, even when empty. The Nextflow module declares
# this output non-optional: relying on an optional output that is genuinely
# absent has proven unreliable on some Nextflow versions (see merge_timeline.R
# and CLAUDE.md).
if (!is.null(spec_out) && nrow(spec_out) > 0) {
  write.table(spec_out, opt$output, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  cat(sprintf("Wrote %d row(s) to %s\n", nrow(spec_out), opt$output))
} else {
  file.create(opt$output)
  cat("No SPECIMEN timeline data to write.\n")
}
