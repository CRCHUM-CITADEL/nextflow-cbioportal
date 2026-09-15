#!/usr/bin/env Rscript
#
# BUILD_CLINICAL_TABLE — reads and merges the ARGO clinical CSVs into one wide
# intermediate table (clinical_merged.tsv), consumed by write_clinical_sample.R
# and write_clinical_patient.R. This is the loaders + merge-cascade half of the
# former clin_format.R; running it once (instead of once per mode) removes the
# duplicated ~95% of the work that used to happen twice per group.
#
# The intermediate table is written with na="NA" and every column forced to
# character on read-back (colClasses="character" in the writers), so no value
# is reformatted (e.g. a numeric "5.0" must not come back as "5") and columns
# absent from a run stay absent — write_cbio_table()'s col_defs filter depends
# on that to reproduce clin_format.R's output exactly.

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--donors", type="character", default=NULL,
              help="donors CSV [REQUIRED]", metavar="FILE"),
  make_option("--primary_diagnoses", type="character", default=NULL,
              help="primary diagnoses CSV [REQUIRED]", metavar="FILE"),
  make_option("--specimens", type="character", default=NULL,
              help="specimens CSV [REQUIRED]", metavar="FILE"),
  make_option("--sample_registrations", type="character", default=NULL,
              help="sample registrations CSV [REQUIRED]", metavar="FILE"),
  make_option("--treatments", type="character", default=NULL,
              help="treatments CSV [OPTIONAL]", metavar="FILE"),
  make_option("--surgeries", type="character", default=NULL,
              help="surgeries CSV [OPTIONAL]", metavar="FILE"),
  make_option("--systemic_therapies", type="character", default=NULL,
              help="systemic therapies CSV [OPTIONAL]", metavar="FILE"),
  make_option("--radiations", type="character", default=NULL,
              help="radiations CSV [OPTIONAL]", metavar="FILE"),
  make_option("--follow_ups", type="character", default=NULL,
              help="follow-ups CSV [OPTIONAL]", metavar="FILE"),
  make_option("--biomarkers", type="character", default=NULL,
              help="biomarkers CSV [OPTIONAL]", metavar="FILE"),
  make_option("--genomic_subjects", type="character", default=NULL,
              help="TSV with subject_id and sample_id columns from genomic pipeline; filters output to these subjects only [OPTIONAL]", metavar="FILE"),
  make_option("--primary_site_map", type="character", default=NULL,
              help="MOHCCN primary site mapping CSV [OPTIONAL]", metavar="FILE"),
  make_option("--specimen_tissue_source_map", type="character", default=NULL,
              help="MOHCCN specimen tissue source mapping CSV [OPTIONAL]", metavar="FILE"),
  make_option("--treatment_intent_map", type="character", default=NULL,
              help="MOHCCN treatment intent mapping CSV [OPTIONAL]", metavar="FILE"),
  make_option(c("-o", "--output"), type="character", default="clinical_merged.tsv",
              help="Output file path [default= %default]", metavar="FILE"),
  make_option("--common_lib", type="character", default="clinical_common.R",
              help="Path to clinical_common.R [default= %default]", metavar="FILE")
)

opt_parser <- OptionParser(
  option_list=option_list,
  usage="Usage: %prog --donors donors.csv --primary_diagnoses diagnoses.csv --specimens specimens.csv --sample_registrations regs.csv [options]",
  description="Merge ARGO clinical data into one intermediate table for the cBioPortal clinical writers"
)
opt <- parse_args(opt_parser)
source(opt$common_lib)

# Required argument validation
if (is.null(opt$donors))               stop("Error: --donors argument is required")
if (is.null(opt$primary_diagnoses))    stop("Error: --primary_diagnoses argument is required")
if (is.null(opt$specimens))            stop("Error: --specimens argument is required")
if (is.null(opt$sample_registrations)) stop("Error: --sample_registrations argument is required")

# File existence checks
required_files <- c(opt$donors, opt$primary_diagnoses, opt$specimens, opt$sample_registrations)
optional_files <- c(opt$treatments, opt$surgeries, opt$systemic_therapies,
                    opt$radiations, opt$follow_ups, opt$biomarkers)
all_optional <- c(optional_files, opt$genomic_subjects)
for (f in c(required_files, all_optional[!sapply(all_optional, is.null)])) {
  if (!file.exists(f)) stop(paste("Error: File does not exist:", f))
}

cat("=== Parameters ===\n")
cat(sprintf("Donors:               %s\n", opt$donors))
cat(sprintf("Primary diagnoses:    %s\n", opt$primary_diagnoses))
cat(sprintf("Specimens:            %s\n", opt$specimens))
cat(sprintf("Sample registrations: %s\n", opt$sample_registrations))
cat(sprintf("Treatments:           %s\n", ifelse(is.null(opt$treatments),        "NULL", opt$treatments)))
cat(sprintf("Surgeries:            %s\n", ifelse(is.null(opt$surgeries),         "NULL", opt$surgeries)))
cat(sprintf("Systemic therapies:   %s\n", ifelse(is.null(opt$systemic_therapies),"NULL", opt$systemic_therapies)))
cat(sprintf("Radiations:           %s\n", ifelse(is.null(opt$radiations),        "NULL", opt$radiations)))
cat(sprintf("Follow-ups:           %s\n", ifelse(is.null(opt$follow_ups),        "NULL", opt$follow_ups)))
cat(sprintf("Biomarkers:           %s\n", ifelse(is.null(opt$biomarkers),        "NULL", opt$biomarkers)))
cat(sprintf("Genomic subjects:     %s\n", ifelse(is.null(opt$genomic_subjects),  "NULL", opt$genomic_subjects)))
cat(sprintf("Output:               %s\n", opt$output))
cat("==================\n\n")

cat("Reading sample registrations...\n")
reg  <- read.csv(opt$sample_registrations, header=TRUE)
link <- reg[is_cbio_sample(reg), ]

# Deduplicate: when multiple samples share the same patient + specimen
# (e.g. -1DT and -2DT for the same biopsy), keep the first by sample ID.
link <- link[order(link$submitter_sample_id), ]
link <- link[!duplicated(link[, c("submitter_donor_id", "submitter_specimen_id")]), ]
link$patient <- link$submitter_donor_id
link$sample  <- link$submitter_sample_id

# Rename sample_type to analyte_type to avoid column name conflict downstream
names(link)[names(link) == "sample_type"] <- "analyte_type"
link <- link[, intersect(c("patient", "sample", "specimen_type", "specimen_tissue_source",
                             "analyte_type", "tumour_normal_designation",
                             "submitter_specimen_id"), names(link))]


# ── Donors ────────────────────────────────────────────────────────────────────
cat("Reading donors...\n")
donors <- read.csv(opt$donors, header=TRUE)
donors$patient   <- donors$submitter_donor_id
donors$sex       <- gsub("Female", "F", gsub("Male", "M", donors$sex_at_birth))
donors$os_status <- ifelse(donors$is_deceased == "Yes", "1:DECEASED", "0:LIVING")
donors$dob_days  <- sapply(donors$date_of_birth, parse_day_interval)
donors$dod_days  <- sapply(donors$date_of_death,  parse_day_interval)
donors <- donors[, intersect(c("patient", "sex", "os_status", "dob_days", "dod_days",
                                 "cause_of_death", "lost_to_followup_reason"), names(donors))]

# ── Primary diagnoses ──────────────────────────────────────────────────────────
cat("Reading primary diagnoses...\n")
diag <- read.csv(opt$primary_diagnoses, header=TRUE)
diag$patient <- diag$submitter_donor_id
diag <- diag[, intersect(c("patient", "cancer_type_code", "clinical_stage_group",
                              "clinical_t_category", "clinical_n_category", "clinical_m_category",
                              "clinical_tumour_staging_system",
                              "pathological_stage_group", "pathological_t_category",
                              "pathological_n_category", "pathological_m_category",
                              "primary_site", "laterality", "basis_of_diagnosis",
                              "submitter_primary_diagnosis_id"), names(diag))]
# Keep every diagnosis: a donor may have more than one, and the right one for a
# sample is the one its sequenced specimen points at (resolved at merge time below).
diag <- diag[order(diag$patient, diag$submitter_primary_diagnosis_id), ]
diag <- diag[!duplicated(diag[, c("patient", "submitter_primary_diagnosis_id")]), ]

# ── Specimens ──────────────────────────────────────────────────────────────────
cat("Reading specimens...\n")
spec <- read.csv(opt$specimens, header=TRUE)
spec <- spec[, intersect(c("submitter_specimen_id", "tumour_histological_type", "tumour_grade",
                              "tumour_grading_system", "specimen_anatomic_location",
                              "specimen_laterality", "specimen_processing", "specimen_storage",
                              "percent_tumour_cells_range",
                              "reference_pathology_confirmed_diagnosis",
                              "submitter_primary_diagnosis_id",
                              "submitter_treatment_id"), names(spec))]

# ── Merge: link ← spec (by specimen ID), then → donors → diagnoses ───────────
m <- merge(link,  spec,   by="submitter_specimen_id", all.x=TRUE)
m <- merge(m,     donors, by="patient", all.x=TRUE)

# Each sample takes the primary diagnosis of its own sequenced specimen, so a donor
# with several diagnoses gets the one that was actually sequenced rather than the
# first one on file. Specimens with no usable diagnosis link fall back to the
# patient's first diagnosis.
if (!("submitter_primary_diagnosis_id" %in% names(m))) m$submitter_primary_diagnosis_id <- NA_character_
diag_keys  <- paste(diag$patient, diag$submitter_primary_diagnosis_id)
first_idx  <- !duplicated(diag$patient)
first_diag <- setNames(diag$submitter_primary_diagnosis_id[first_idx], diag$patient[first_idx])
fallback   <- is.na(m$submitter_primary_diagnosis_id) |
              !(paste(m$patient, m$submitter_primary_diagnosis_id) %in% diag_keys)
m$submitter_primary_diagnosis_id[fallback] <- first_diag[m$patient[fallback]]

m <- merge(m, diag, by=c("patient", "submitter_primary_diagnosis_id"), all.x=TRUE)

# ── Filter to genomic subjects (optional) ────────────────────────────────────
# When the genomic pipeline provides a linking file (subject_id TAB sample_id),
# restrict output to those subjects and adopt the genomic sample_id convention.
# Clinical submitter_sample_id values (e.g. MoHQ-CM-1-100-76994-1DT) differ from
# the genomic sample_id (e.g. MoHQ-CM-1-100-T); the linking file is authoritative.
if (!is.null(opt$genomic_subjects)) {
  cat("Reading genomic subjects for filtering...\n")
  genomic <- read_genomic_subjects(opt$genomic_subjects)
  m <- m[m$patient %in% genomic$subject_id, ]
  if (nrow(m) == 0) stop("Error: No clinical subjects match the genomic subjects file.")
  # Only remap sample IDs for patients with a single clinical sample.
  # Multi-sample patients keep their clinical sample IDs to avoid duplicate SAMPLE_IDs.
  sample_map    <- setNames(genomic$sample_id, genomic$subject_id)
  sample_counts <- table(m$patient)
  single        <- names(sample_counts[sample_counts == 1])
  idx           <- m$patient %in% single
  m$sample[idx] <- sample_map[m$patient[idx]]
}


# ── Treatments (primary treatment only, optional) ─────────────────────────────
if (!is.null(opt$treatments)) {
  cat("Reading treatments...\n")
  treat <- read.csv(opt$treatments, header=TRUE)
  treat <- treat[!is.na(treat$is_primary_treatment) & treat$is_primary_treatment == "Yes", ]
  treat$patient             <- treat$submitter_donor_id
  treat$treatment_type_str  <- sapply(treat$treatment_type, clean_json_array)
  treat <- treat[, intersect(c("patient", "treatment_intent", "treatment_type_str",
                                 "response_to_treatment", "status_of_treatment"), names(treat))]
  treat <- treat[!duplicated(treat$patient), ]
  m <- merge(m, treat, by="patient", all.x=TRUE)
}

# ── Surgeries (optional) ─────────────────────────────────────────────────────
if (!is.null(opt$surgeries)) {
  cat("Reading surgeries...\n")
  surg <- read.csv(opt$surgeries, header=TRUE)
  surg$patient <- surg$submitter_donor_id
  surg <- surg[, intersect(c("patient", "surgery_type", "surgery_location",
                               "residual_tumour_classification", "lymphovascular_invasion",
                               "perineural_invasion", "tumour_focality",
                               "greatest_dimension_tumour"), names(surg))]
  surg <- surg[!duplicated(surg$patient), ]
  m <- merge(m, surg, by="patient", all.x=TRUE)
}

# ── Systemic therapies (optional) ─────────────────────────────────────────────
if (!is.null(opt$systemic_therapies)) {
  cat("Reading systemic therapies...\n")
  syst <- read.csv(opt$systemic_therapies, header=TRUE)
  syst$patient <- syst$submitter_donor_id
  if ("number_of_cycles_not_available" %in% names(syst)) {
    syst$number_of_cycles[syst$number_of_cycles_not_available == "True"] <- NA
  }
  drug_agg <- tapply(syst$drug_name,             syst$patient,
                     function(x) paste(unique(x[!is.na(x)]), collapse=";"))
  type_agg <- tapply(syst$systemic_therapy_type,  syst$patient,
                     function(x) paste(unique(x[!is.na(x)]), collapse=";"))
  syst_first <- syst[!duplicated(syst$patient), c("patient", "number_of_cycles")]
  syst_sum <- data.frame(
    patient               = names(drug_agg),
    drug_name             = as.character(drug_agg),
    systemic_therapy_type = as.character(type_agg),
    stringsAsFactors      = FALSE
  )
  syst_sum <- merge(syst_sum, syst_first, by="patient", all.x=TRUE)
  m <- merge(m, syst_sum, by="patient", all.x=TRUE)
}

# ── Radiations (optional) ────────────────────────────────────────────────────
if (!is.null(opt$radiations)) {
  cat("Reading radiations...\n")
  radi <- read.csv(opt$radiations, header=TRUE)
  radi$patient <- radi$submitter_donor_id
  if ("radiation_therapy_dosage_not_available" %in% names(radi)) {
    radi$radiation_therapy_dosage[radi$radiation_therapy_dosage_not_available == "True"] <- NA
  }
  if ("radiation_therapy_fractions_not_available" %in% names(radi)) {
    radi$radiation_therapy_fractions[radi$radiation_therapy_fractions_not_available == "True"] <- NA
  }
  radi <- radi[, intersect(c("patient", "anatomical_site_irradiated",
                               "radiation_therapy_modality", "radiation_therapy_type",
                               "radiation_therapy_dosage", "radiation_therapy_fractions"),
                             names(radi))]
  radi <- radi[!duplicated(radi$patient), ]
  m <- merge(m, radi, by="patient", all.x=TRUE)
}

# ── Follow-ups (OS for living patients, DFS, last disease status) ─────────────
last_fu_days              <- c()
first_relapse_days        <- c()
last_disease_status       <- c()
first_relapse_type        <- c()
first_relapse_site        <- c()
last_visit_relapse_type   <- c()
last_method_of_progression <- c()

if (!is.null(opt$follow_ups)) {
  cat("Reading follow-ups...\n")
  fu <- read.csv(opt$follow_ups, header=TRUE)
  fu$patient       <- fu$submitter_donor_id
  fu$followup_days <- sapply(fu$date_of_followup, parse_day_interval)
  fu$relapse_days  <- sapply(fu$date_of_relapse,  parse_day_interval)

  fu_dated <- fu[!is.na(fu$followup_days), ]
  if (nrow(fu_dated) > 0) {
    tmp <- tapply(fu_dated$followup_days, fu_dated$patient, function(x) max(x, na.rm=TRUE))
    last_fu_days <- tmp[!is.infinite(tmp) & !is.na(tmp)]
    fu_sorted <- fu_dated[order(fu_dated$patient, fu_dated$followup_days), ]
    fu_last   <- fu_sorted[!duplicated(fu_sorted$patient, fromLast=TRUE), ]
    if ("disease_status_at_followup" %in% names(fu_last)) {
      last_disease_status <- setNames(fu_last$disease_status_at_followup, fu_last$patient)
    }
    if ("relapse_type" %in% names(fu_last)) {
      last_visit_relapse_type <- setNames(fu_last$relapse_type, fu_last$patient)
    }
    if ("method_of_progression_status" %in% names(fu_last)) {
      last_method_of_progression <- setNames(fu_last$method_of_progression_status, fu_last$patient)
    }
  }

  fu_rel <- fu[!is.na(fu$relapse_days), ]
  if (nrow(fu_rel) > 0) {
    fu_rel <- fu_rel[order(fu_rel$patient, fu_rel$relapse_days), ]
    fu_rel <- fu_rel[!duplicated(fu_rel$patient), ]
    first_relapse_days <- setNames(fu_rel$relapse_days, fu_rel$patient)
    if ("relapse_type" %in% names(fu_rel)) {
      first_relapse_type <- setNames(fu_rel$relapse_type, fu_rel$patient)
    }
    if ("anatomic_site_progression_or_recurrence" %in% names(fu_rel)) {
      first_relapse_site <- setNames(fu_rel$anatomic_site_progression_or_recurrence, fu_rel$patient)
    }
  }
}

# ── Biomarkers (most recent per patient, optional) ────────────────────────────
if (!is.null(opt$biomarkers)) {
  cat("Reading biomarkers...\n")
  biom <- read.csv(opt$biomarkers, header=TRUE)
  biom$patient  <- biom$submitter_donor_id
  biom$test_day <- sapply(biom$test_date, parse_day_interval)
  biom <- biom[order(biom$patient, biom$test_day, decreasing=TRUE), ]
  biom <- biom[!duplicated(biom$patient), ]
  biom <- biom[, intersect(c("patient", "er_status", "er_percent_positive",
                               "pr_status", "pr_percent_positive",
                               "her2_ihc_status", "her2_ish_status",
                               "hpv_ihc_status", "hpv_pcr_status", "hpv_strain",
                               "ca125", "cea", "psa_level"), names(biom))]
  m <- merge(m, biom, by="patient", all.x=TRUE)
}

# ── Computed clinical fields ───────────────────────────────────────────────────
# Age at diagnosis: dob_days is negative (days before day-0 diagnosis)
m$age <- round(-m$dob_days / 365.25, 3)

# OS months: deceased from date_of_death; living from last follow-up (or NA)
m$os_months <- NA_real_
deceased_idx <- !is.na(m$os_status) & m$os_status == "1:DECEASED"
m$os_months[deceased_idx] <- round(m$dod_days[deceased_idx] / 30.4375, 3)

living_idx <- !is.na(m$os_status) & m$os_status == "0:LIVING"
if (any(living_idx) && length(last_fu_days) > 0) {
  fu_match <- last_fu_days[m$patient[living_idx]]
  valid <- !is.na(fu_match)
  m$os_months[living_idx][valid] <- round(as.numeric(fu_match[valid]) / 30.4375, 3)
}

# DFS: based on last visit relapse_type — recurred unless empty or "Not available"
m$dfs_status <- NA_character_
m$dfs_months <- NA_real_
if (length(last_visit_relapse_type) > 0) {
  has_fu <- m$patient %in% names(last_visit_relapse_type)
  rtype  <- last_visit_relapse_type[m$patient[has_fu]]
  is_free <- is.na(rtype) | trimws(rtype) == "" | rtype == "Not available"

  m$dfs_status[has_fu][is_free]  <- "0:DiseaseFree"
  m$dfs_status[has_fu][!is_free] <- "1:Recurred/Progressed"

  # DFS_MONTHS: disease-free → last follow-up date
  if (length(last_fu_days) > 0) {
    free_with_fu <- has_fu & m$dfs_status == "0:DiseaseFree" &
                    m$patient %in% names(last_fu_days)
    m$dfs_months[free_with_fu] <- round(
      as.numeric(last_fu_days[m$patient[free_with_fu]]) / 30.4375, 3
    )
  }
  # DFS_MONTHS: recurred → first date_of_relapse, fallback to last follow-up date
  recurred <- has_fu & m$dfs_status == "1:Recurred/Progressed"
  if (length(first_relapse_days) > 0) {
    has_relapse_date <- recurred & m$patient %in% names(first_relapse_days)
    m$dfs_months[has_relapse_date] <- round(
      as.numeric(first_relapse_days[m$patient[has_relapse_date]]) / 30.4375, 3
    )
    no_relapse_date <- recurred & !m$patient %in% names(first_relapse_days) &
                       m$patient %in% names(last_fu_days)
    m$dfs_months[no_relapse_date] <- round(
      as.numeric(last_fu_days[m$patient[no_relapse_date]]) / 30.4375, 3
    )
  } else if (length(last_fu_days) > 0) {
    recurred_with_fu <- recurred & m$patient %in% names(last_fu_days)
    m$dfs_months[recurred_with_fu] <- round(
      as.numeric(last_fu_days[m$patient[recurred_with_fu]]) / 30.4375, 3
    )
  }
} else if (length(last_fu_days) > 0) {
  # No relapse_type column at all but have follow-ups → all disease-free
  has_fu <- m$patient %in% names(last_fu_days)
  m$dfs_status[has_fu] <- "0:DiseaseFree"
  m$dfs_months[has_fu] <- round(
    as.numeric(last_fu_days[m$patient[has_fu]]) / 30.4375, 3
  )
}

# Follow-up derived status fields
m$disease_status_at_followup <- NA_character_
m$relapse_type_val            <- NA_character_
m$relapse_site_val            <- NA_character_
if (length(last_disease_status) > 0) {
  has_status <- m$patient %in% names(last_disease_status)
  m$disease_status_at_followup[has_status] <- last_disease_status[m$patient[has_status]]
}
if (length(first_relapse_type) > 0) {
  has_rtype <- m$patient %in% names(first_relapse_type)
  m$relapse_type_val[has_rtype] <- first_relapse_type[m$patient[has_rtype]]
}
if (length(first_relapse_site) > 0) {
  has_rsite <- m$patient %in% names(first_relapse_site)
  m$relapse_site_val[has_rsite] <- clean_json_array(first_relapse_site[m$patient[has_rsite]])
}
m$relapse_site_label <- NA_character_
m$method_of_progression <- NA_character_
if (length(last_method_of_progression) > 0) {
  has_method <- m$patient %in% names(last_method_of_progression)
  m$method_of_progression[has_method] <- clean_json_array(last_method_of_progression[m$patient[has_method]])
}

# Sample type: map ARGO specimen_type values to cBioPortal vocabulary
sample_type_map <- c(
  "Primary tumour"    = "Primary",
  "Metastatic tumour" = "Metastasis",
  "Recurrent tumour"  = "Recurrence"
)
if ("specimen_type" %in% names(m)) {
  m$sample_type_mapped <- sample_type_map[m$specimen_type]
  m$sample_type_mapped[is.na(m$sample_type_mapped)] <- m$specimen_type[is.na(m$sample_type_mapped)]
} else {
  m$sample_type_mapped <- "Primary"
}

m[m == ""] <- NA

# ── MOHCCN ontology code lookups ──────────────────────────────────────────────
ps_map  <- read_mohccn_map(opt$primary_site_map)
ps_map_rev <- reverse_mohccn_map(ps_map)
sts_map <- read_mohccn_map(opt$specimen_tissue_source_map)
ti_map  <- read_mohccn_map(opt$treatment_intent_map)

m$primary_site_code           <- apply_mohccn_map(m$primary_site,           ps_map)
m$specimen_tissue_source_code <- apply_mohccn_map(m$specimen_tissue_source, sts_map)
m$relapse_site_label          <- apply_ps_label(m$relapse_site_val, ps_map_rev)
if ("treatment_intent" %in% names(m)) {
  m$treatment_intent_code <- apply_mohccn_map(m$treatment_intent, ti_map)
} else {
  m$treatment_intent_code <- NA_character_
}

cat(sprintf("Writing merged clinical table (%d rows, %d columns)...\n", nrow(m), ncol(m)))
write.table(m, opt$output, sep="\t", row.names=FALSE, quote=FALSE, na="NA")
