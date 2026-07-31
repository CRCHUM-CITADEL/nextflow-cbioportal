#!/usr/bin/env Rscript

library(optparse)

option_list <- list(
  make_option("--sample_registrations", type="character", default=NULL,
              help="sample registrations CSV [OPTIONAL]", metavar="FILE"),
  make_option("--treatments",           type="character", default=NULL,
              help="treatments CSV [OPTIONAL]",            metavar="FILE"),
  make_option("--surgeries",            type="character", default=NULL,
              help="surgeries CSV [OPTIONAL]",             metavar="FILE"),
  make_option("--systemic_therapies",   type="character", default=NULL,
              help="systemic therapies CSV [OPTIONAL]",    metavar="FILE"),
  make_option("--follow_ups",           type="character", default=NULL,
              help="follow-ups CSV [OPTIONAL]",            metavar="FILE"),
  make_option("--specimens",            type="character", default=NULL,
              help="specimens CSV [OPTIONAL]",             metavar="FILE"),
  make_option("--biomarkers",           type="character", default=NULL,
              help="biomarkers CSV [OPTIONAL]",            metavar="FILE"),
  make_option("--genomic_subjects",           type="character", default=NULL,
              help="TSV with subject_id and sample_id columns; filters output to genomic subjects [OPTIONAL]", metavar="FILE"),
  make_option("--primary_site_map",           type="character", default=NULL,
              help="MOHCCN primary site mapping CSV (text → ICD-O code) [OPTIONAL]",  metavar="FILE"),
  make_option("--specimen_tissue_source_map", type="character", default=NULL,
              help="MOHCCN specimen tissue source mapping CSV [OPTIONAL]", metavar="FILE"),
  make_option("--treatment_intent_map",       type="character", default=NULL,
              help="MOHCCN treatment intent mapping CSV [OPTIONAL]",       metavar="FILE")
)

opt_parser <- OptionParser(
  option_list = option_list,
  usage = "Usage: %prog [options]",
  description = "Generate cBioPortal timeline files from ARGO clinical CSVs"
)
opt <- parse_args(opt_parser)

# File existence checks
for (f in Filter(Negate(is.null), opt[names(opt) != "help"])) {
  if (is.character(f) && nchar(f) > 0 && !file.exists(f)) {
    stop(paste("Error: File does not exist:", f))
  }
}

# Read MOHCCN mapping CSV → named vector: plain-text value → ontology/ICD-O code.
# Skips empty rows; reads only first 3 cols (treatment_intent.csv has trailing commas).
read_mohccn_map <- function(filepath) {
  if (is.null(filepath) || !file.exists(filepath)) return(NULL)
  raw <- read.csv(filepath, header=FALSE, skip=1, stringsAsFactors=FALSE, fill=TRUE)
  raw <- raw[, 1:3, drop=FALSE]
  raw <- raw[nchar(trimws(raw[[2]])) > 0, ]
  setNames(trimws(raw[[3]]), trimws(raw[[2]]))
}

# Look up values in a MOHCCN map; returns NA_character_ where there is no match or empty code.
apply_mohccn_map <- function(values, map) {
  if (is.null(map)) return(rep(NA_character_, length(values)))
  result <- map[as.character(values)]
  result[is.na(result) | result == ""] <- NA_character_
  unname(result)
}

cat("=== Parameters ===\n")
for (nm in c("sample_registrations", "treatments", "surgeries",
             "systemic_therapies", "follow_ups", "specimens", "biomarkers",
             "genomic_subjects", "primary_site_map", "specimen_tissue_source_map", "treatment_intent_map")) {
  cat(sprintf("%-22s %s\n", paste0(nm, ":"), ifelse(is.null(opt[[nm]]), "NULL", opt[[nm]])))
}
cat("==================\n\n")

# Extract day_interval from ARGO JSON date strings: {"day_interval": N, "month_interval": M}
parse_day_interval <- function(x) {
  x <- as.character(x)
  m <- regmatches(x, regexpr('"day_interval"\\s*:\\s*(-?[0-9]+)', x))
  if (length(m) == 0 || nchar(m) == 0) return(NA_real_)
  as.numeric(sub('"day_interval"\\s*:\\s*(-?[0-9]+)', '\\1', m))
}

# Flatten a JSON array string: ["Surgery"] -> "Surgery"
clean_json_array <- function(x) {
  x <- gsub('\\[|\\]', '', x)
  x <- gsub('"', '', x)
  trimws(x)
}

# Write a timeline data frame to file; skip if empty
write_timeline <- function(df, filename) {
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  write.table(df, filename, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  cat(sprintf("Wrote %d row(s) to %s\n", nrow(df), filename))
}

# ── Valid patient IDs + specimen maps from sample registrations ────────────────
valid_patients         <- character(0)
spec_type_map          <- character(0)   # submitter_specimen_id -> specimen_type
spec_tissue_source_map <- character(0)   # submitter_specimen_id -> specimen_tissue_source

if (!is.null(opt$sample_registrations)) {
  cat("Reading sample registrations...\n")
  reg <- read.csv(opt$sample_registrations, header = TRUE, stringsAsFactors = FALSE)
  tumour_reg <- reg[
    reg$tumour_normal_designation == "Tumour" &
    reg$sample_type              == "Total DNA" &
    reg$specimen_tissue_source   == "Solid tissue" &
    grepl("-\\d+[DR]T$", reg$submitter_sample_id), ]
  valid_patients <- unique(tumour_reg$submitter_donor_id)
  # Deduplicate by specimen_id — multiple sequencing types (DNA, RNA) per specimen produce
  # duplicate rows; take the first entry so the named-vector lookup stays unambiguous.
  reg_dedup <- reg[!duplicated(reg$submitter_specimen_id), ]
  if ("specimen_type" %in% names(reg_dedup)) {
    spec_type_map <- setNames(reg_dedup$specimen_type, reg_dedup$submitter_specimen_id)
  }
  if ("specimen_tissue_source" %in% names(reg_dedup)) {
    spec_tissue_source_map <- setNames(reg_dedup$specimen_tissue_source, reg_dedup$submitter_specimen_id)
  }
}

# When a genomic linking file is provided, restrict to those subjects (same as clin_format.R)
if (!is.null(opt$genomic_subjects)) {
  cat("Reading genomic subjects for filtering...\n")
  genomic        <- read.table(opt$genomic_subjects, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  valid_patients <- genomic$subject_id
}

# Load MOHCCN mapping tables.
# ps_map is REVERSED (ICD-O code → label) since surgery_site / specimen_anatomic_location
# are already ICD-O codes in the data; we add a human-readable label column.
ps_map_raw <- read_mohccn_map(opt$primary_site_map)   # label → code
ps_map     <- if (!is.null(ps_map_raw)) setNames(names(ps_map_raw), ps_map_raw) else NULL
sts_map    <- read_mohccn_map(opt$specimen_tissue_source_map)
ti_map     <- read_mohccn_map(opt$treatment_intent_map)

# Apply reversed primary-site map: ICD-O code (possibly with sub-code like "C22.0")
# → human-readable label. Strips the decimal suffix before lookup ("C22.0" → "C22").
apply_ps_label <- function(codes, reversed_map) {
  if (is.null(reversed_map)) return(rep(NA_character_, length(codes)))
  prefixes <- sub("\\..*", "", trimws(as.character(codes)))
  result   <- reversed_map[prefixes]
  result[is.na(result) | result == ""] <- NA_character_
  unname(result)
}

# Helper: filter to valid patients when the set is non-empty
filter_patients <- function(df, id_col = "submitter_donor_id") {
  if (length(valid_patients) > 0) df[df[[id_col]] %in% valid_patients, ] else df
}

# Derive cBioPortal SAMPLE_TYPE from the nearest follow-up disease status.
# Colors: Recurrence/Progression → orange; Metastasis → red; otherwise → black.
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

# ── Pre-read follow-ups (used in both section c STATUS and section d SPECIMEN) ─────────
fu_data <- NULL
if (!is.null(opt$follow_ups)) {
  cat("Reading follow-ups...\n")
  fu_raw           <- read.csv(opt$follow_ups, header = TRUE, stringsAsFactors = FALSE)
  fu_raw$start_day <- sapply(fu_raw$date_of_followup, parse_day_interval)
  fu_data          <- fu_raw[!is.na(fu_raw$start_day), ]
  fu_data          <- filter_patients(fu_data)
}

# ── (a) Surgery + (b) Systemic treatment timelines (both sourced from treatments.csv) ──

if (!is.null(opt$treatments)) {
  cat("Reading treatments...\n")
  treat <- read.csv(opt$treatments, header = TRUE, stringsAsFactors = FALSE)
  treat <- treat[!is.na(treat$is_primary_treatment) & treat$is_primary_treatment == "Yes", ]
  treat$start_day        <- sapply(treat$treatment_start_date, parse_day_interval)
  treat$stop_day         <- sapply(treat$treatment_end_date,   parse_day_interval)
  treat$treatment_type_c <- sapply(treat$treatment_type,       clean_json_array)
  treat <- filter_patients(treat)

  # ── (a) Surgery timeline ──────────────────────────────────────────────────────
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
    write_timeline(surg_out, "data_timeline_surgery.txt")
  }

  # ── (b) Systemic treatment timeline ──────────────────────────────────────────
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
    write_timeline(syst_out, "data_timeline_treatment.txt")
  }
}

# ── (c) Status timeline (follow-up visits) ────────────────────────────────────
if (!is.null(fu_data) && nrow(fu_data) > 0) {
  status_out <- data.frame(
    PATIENT_ID = fu_data$submitter_donor_id,
    START_DATE = fu_data$start_day,
    STOP_DATE  = "",
    EVENT_TYPE = "STATUS",
    STATUS     = if ("disease_status_at_followup" %in% names(fu_data)) fu_data$disease_status_at_followup else NA_character_,
    stringsAsFactors = FALSE
  )
  write_timeline(status_out, "data_timeline_status.txt")
}

# ── (d) Specimen timeline ──────────────────────────────────────────────────────
if (!is.null(opt$specimens)) {
  cat("Reading specimens...\n")
  spec <- read.csv(opt$specimens, header = TRUE, stringsAsFactors = FALSE)
  spec$start_day <- sapply(spec$specimen_collection_date, parse_day_interval)
  spec <- spec[!is.na(spec$start_day), ]
  spec <- filter_patients(spec)

  if (nrow(spec) > 0) {
    # Look up specimen_type and specimen_tissue_source from sample_registrations
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
    write_timeline(spec_out, "data_timeline_specimen.txt")
  }
}

# ── (e) Lab test timeline (biomarkers pivoted to long format) ─────────────────
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
    write_timeline(lab_out, "data_timeline_lab_test.txt")
  }
}
