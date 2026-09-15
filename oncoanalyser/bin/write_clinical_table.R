#!/usr/bin/env Rscript
#
# WRITE_CLINICAL_TABLE — applies the cBioPortal column definitions to the
# intermediate table from build_clinical_table.R and writes either
# data_clinical_sample.txt or data_clinical_patient.txt, selected by --mode.
# This is the col_defs + write_cbio_table() half of the former clin_format.R;
# it does no merging of its own, so running it twice (once per mode) is cheap.

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--input", type="character", default=NULL,
              help="Merged clinical table TSV from BUILD_CLINICAL_TABLE [REQUIRED]", metavar="FILE"),
  make_option(c("-o", "--output"), type="character", default="data_clinical_sample.txt",
              help="Output file path [default= %default]", metavar="FILE"),
  make_option(c("-m", "--mode"), type="character", default="sample",
              help="Between 'sample' or 'patient' mode.")
)

opt_parser <- OptionParser(
  option_list=option_list,
  usage="Usage: %prog --input clinical_merged.tsv --mode sample|patient [options]",
  description="Write a cBioPortal clinical attribute file from the merged clinical table"
)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) stop("Error: --input argument is required")
if (!file.exists(opt$input)) stop(paste("Error: File does not exist:", opt$input))

cat(sprintf("Reading merged clinical table from %s...\n", opt$input))
# na.strings="NA" + colClasses="character": round-trips build_clinical_table.R's
# output exactly, with no numeric reformatting and no column dropped or coerced.
m <- read.table(opt$input, sep="\t", header=TRUE, na.strings="NA",
                 colClasses="character", check.names=FALSE, quote="")

cat(paste("Writing", opt$mode, "mode output...\n"))

# ── Column definition helper ───────────────────────────────────────────────────
# Each entry: list(cbio_id, r_col, display_name, description, datatype, priority)
# Columns whose r_col is absent from m are silently dropped.
write_cbio_table <- function(m, col_defs, output_file) {
  col_defs <- Filter(function(cd) cd[[2]] %in% names(m), col_defs)
  cbio_ids   <- sapply(col_defs, `[[`, 1)
  r_cols     <- sapply(col_defs, `[[`, 2)
  display    <- sapply(col_defs, `[[`, 3)
  desc       <- sapply(col_defs, `[[`, 4)
  types      <- sapply(col_defs, `[[`, 5)
  priorities <- sapply(col_defs, `[[`, 6)
  writeLines(c(
    paste0("#", paste(display,    collapse="\t")),
    paste0("#", paste(desc,       collapse="\t")),
    paste0("#", paste(types,      collapse="\t")),
    paste0("#", paste(priorities, collapse="\t")),
    paste(cbio_ids, collapse="\t")
  ), con=output_file)
  write.table(
    m[, r_cols, drop=FALSE],
    file=output_file, sep="\t", row.names=FALSE, col.names=FALSE,
    quote=FALSE, append=TRUE, na="NA"
  )
}

# ── Write output ───────────────────────────────────────────────────────────────
if (opt$mode == "patient") {
  # Patient-level: only columns sourced from donors.csv plus derived survival fields.
  # All other attributes (diagnosis, specimen, treatment, biomarker) are sample-level.
  col_defs <- list(
    list("PATIENT_ID",                 "patient",                        "Patient Identifier",          "Identifier to uniquely specify a patient.",                           "STRING", "1"),
    list("SEX",                        "sex",                            "Sex",                         "Sex of the patient.",                                                  "STRING", "1"),
    list("AGE",                        "age",                            "Diagnosis Age",               "Age at which a condition or disease was first diagnosed.",             "NUMBER", "1"),
    list("CAUSE_OF_DEATH",             "cause_of_death",                 "Cause of Death",              "Cause of patient death.",                                              "STRING", "1"),
    list("OS_MONTHS",                  "os_months",                      "Overall Survival (Months)",   "Overall survival in months since initial diagnosis.",                  "NUMBER", "1"),
    list("OS_STATUS",                  "os_status",                      "Overall Survival Status",     "Overall patient survival status.",                                     "STRING", "1"),
    list("DFS_MONTHS",                 "dfs_months",                     "Disease Free (Months)",       "Disease free in months since initial treatment.",                      "NUMBER", "1"),
    list("DFS_STATUS",                 "dfs_status",                     "Disease Free Status",         "Disease free status since initial treatment.",                         "STRING", "1"),
    list("DISEASE_STATUS_AT_FOLLOWUP", "disease_status_at_followup",     "Disease Status at Follow-up", "Patient disease status at the most recent follow-up.",                "STRING", "1"),
    list("RELAPSE_TYPE",               "relapse_type_val",               "Relapse Type",                "Type of disease relapse or recurrence.",                               "STRING", "1"),
    list("RELAPSE_SITE",               "relapse_site_val",               "Relapse Site",                "Anatomic site of disease progression or recurrence (ICD-O code).",     "STRING", "1"),
    list("RELAPSE_SITE_LABEL",         "relapse_site_label",             "Relapse Site Label",           "Human-readable label for the anatomic site of progression.",           "STRING", "1"),
    list("METHOD_OF_PROGRESSION_STATUS", "method_of_progression",        "Method of Progression Status", "Method used to assess disease progression status.",                    "STRING", "1")
  )
  m_patient <- m[!duplicated(m$patient), ]
  write_cbio_table(m_patient, col_defs, opt$output)

} else if (opt$mode == "sample") {
  # Sample-level: all columns NOT sourced from donors.csv.
  # Includes diagnosis, staging, specimen, treatment, surgery, radiation, and biomarker fields.
  col_defs <- list(
    list("PATIENT_ID",               "patient",                            "Patient Identifier",          "Identifier to uniquely specify a patient.",                            "STRING", "1"),
    list("SAMPLE_ID",                "sample",                             "Sample Identifier",           "A unique sample identifier.",                                          "STRING", "1"),
    list("TUMOUR_NORMAL_DESIGNATION","tumour_normal_designation",          "Tumour/Normal Designation",   "Whether the sample is from tumour or normal tissue.",                  "STRING", "1"),
    list("SAMPLE_TYPE",              "sample_type_mapped",                 "Sample Type",                 "The type of sample (e.g., Primary, Metastasis, Recurrence).",         "STRING", "1"),
    list("SPECIMEN_TISSUE_SOURCE",      "specimen_tissue_source",             "Specimen Tissue Source",      "Tissue source of the specimen.",                                       "STRING", "1"),
    list("SPECIMEN_TISSUE_SOURCE_CODE","specimen_tissue_source_code",        "Specimen Tissue Source Code", "Ontology code for the specimen tissue source.",                        "STRING", "1"),
    list("CANCER_TYPE_CODE",           "cancer_type_code",                   "Cancer Code",                 "ICD-O-3 topography cancer code.",                                      "STRING", "1"),
    list("PRIMARY_SITE",               "primary_site",                       "Primary Site",                "Primary site of the tumor.",                                           "STRING", "1"),
    list("PRIMARY_SITE_CODE",          "primary_site_code",                  "Primary Site Code",           "ICD-O topography code for the primary site.",                          "STRING", "1"),
    list("LATERALITY",               "laterality",                         "Laterality",                  "Laterality of the primary tumor.",                                     "STRING", "1"),
    list("BASIS_OF_DIAGNOSIS",       "basis_of_diagnosis",                 "Basis of Diagnosis",          "Basis on which the primary diagnosis was made.",                       "STRING", "1"),
    list("CLINICAL_STAGE",           "clinical_stage_group",               "Clinical Stage",              "Clinical stage group of the tumor.",                                   "STRING", "1"),
    list("CLINICAL_T_CATEGORY",      "clinical_t_category",                "Clinical T Category",         "Clinical T category (TNM staging).",                                   "STRING", "1"),
    list("CLINICAL_N_CATEGORY",      "clinical_n_category",                "Clinical N Category",         "Clinical N category (TNM staging).",                                   "STRING", "1"),
    list("CLINICAL_M_CATEGORY",      "clinical_m_category",                "Clinical M Category",         "Clinical M category (TNM staging).",                                   "STRING", "1"),
    list("STAGING_SYSTEM",           "clinical_tumour_staging_system",     "Staging System",              "Tumor staging system used.",                                           "STRING", "1"),
    list("PATHOLOGICAL_STAGE",       "pathological_stage_group",           "Pathological Stage",          "Pathological stage group of the tumor.",                               "STRING", "1"),
    list("PATHOLOGICAL_T_CATEGORY",  "pathological_t_category",            "Pathological T Category",     "Pathological T category (TNM staging).",                               "STRING", "1"),
    list("PATHOLOGICAL_N_CATEGORY",  "pathological_n_category",            "Pathological N Category",     "Pathological N category (TNM staging).",                               "STRING", "1"),
    list("PATHOLOGICAL_M_CATEGORY",  "pathological_m_category",            "Pathological M Category",     "Pathological M category (TNM staging).",                               "STRING", "1"),
    list("TUMOR_TISSUE_SITE",        "specimen_anatomic_location",         "Tumor Tissue Site",           "Anatomic location of the specimen.",                                   "STRING", "1"),
    list("SPECIMEN_LATERALITY",      "specimen_laterality",                "Specimen Laterality",         "Laterality of the specimen collection site.",                          "STRING", "1"),
    list("SPECIMEN_PROCESSING",      "specimen_processing",                "Specimen Processing",         "Method used to process the specimen.",                                 "STRING", "1"),
    list("SPECIMEN_STORAGE",         "specimen_storage",                   "Specimen Storage",            "Method used to store the specimen.",                                   "STRING", "1"),
    list("TUMOR_HISTOLOGICAL_TYPE",  "tumour_histological_type",           "Tumor Histological Type",     "Tumor histological type (ICD-O morphology code).",                    "STRING", "1"),
    list("TUMOR_GRADE",              "tumour_grade",                       "Tumor Grade",                 "Tumor grade classification.",                                          "STRING", "1"),
    list("TUMOR_GRADING_SYSTEM",     "tumour_grading_system",              "Tumor Grading System",        "System used for tumor grading.",                                       "STRING", "1"),
    list("TUMOR_CELLS_RANGE",        "percent_tumour_cells_range",         "Tumor Cells Range",           "Percentage range of tumor cells in the specimen.",                    "STRING", "1"),
    list("TREATMENT_TYPE",           "treatment_type_str",                 "Treatment Type",              "Type of primary treatment received.",                                  "STRING", "1"),
    list("TREATMENT_INTENT",          "treatment_intent",                   "Treatment Intent",            "Intent of the primary treatment (e.g., Curative, Palliative).",       "STRING", "1"),
    list("TREATMENT_INTENT_CODE",    "treatment_intent_code",              "Treatment Intent Code",       "Ontology code for the treatment intent.",                              "STRING", "1"),
    list("TREATMENT_RESPONSE",       "response_to_treatment",              "Treatment Response",          "Patient response to treatment.",                                        "STRING", "1"),
    list("TREATMENT_STATUS",         "status_of_treatment",                "Treatment Status",            "Current status of the treatment.",                                     "STRING", "1"),
    list("SURGERY_TYPE",             "surgery_type",                       "Surgery Type",                "Type of surgical procedure performed.",                                 "STRING", "1"),
    list("SURGERY_LOCATION",         "surgery_location",                   "Surgery Location",            "Location of the surgical procedure.",                                   "STRING", "1"),
    list("RESIDUAL_TUMOR",           "residual_tumour_classification",     "Residual Tumor",              "Residual tumor classification after surgery (R classification).",      "STRING", "1"),
    list("LYMPHOVASCULAR_INVASION",  "lymphovascular_invasion",            "Lymphovascular Invasion",     "Presence of lymphovascular invasion.",                                  "STRING", "1"),
    list("PERINEURAL_INVASION",      "perineural_invasion",                "Perineural Invasion",         "Presence of perineural invasion.",                                     "STRING", "1"),
    list("TUMOR_FOCALITY",           "tumour_focality",                    "Tumor Focality",              "Focality of the tumor (unifocal or multifocal).",                      "STRING", "1"),
    list("DRUG_NAME",                "drug_name",                          "Drug Name",                   "Name(s) of systemic therapy drug(s) administered.",                   "STRING", "1"),
    list("SYSTEMIC_THERAPY_TYPE",    "systemic_therapy_type",              "Systemic Therapy Type",       "Type of systemic therapy (e.g., Chemotherapy, Immunotherapy).",       "STRING", "1"),
    list("NUMBER_OF_CYCLES",         "number_of_cycles",                   "Number of Cycles",            "Number of systemic therapy cycles administered.",                      "NUMBER", "1"),
    list("RADIATION_SITE",           "anatomical_site_irradiated",         "Radiation Site",              "Anatomical site that received radiation therapy.",                     "STRING", "1"),
    list("RADIATION_MODALITY",       "radiation_therapy_modality",         "Radiation Modality",          "Modality of radiation therapy administered.",                           "STRING", "1"),
    list("RADIATION_TYPE",           "radiation_therapy_type",             "Radiation Type",              "Type of radiation therapy administered.",                               "STRING", "1"),
    list("RADIATION_DOSAGE",         "radiation_therapy_dosage",           "Radiation Dosage (Gy)",       "Total radiation therapy dosage in Gray.",                              "NUMBER", "1"),
    list("RADIATION_FRACTIONS",      "radiation_therapy_fractions",        "Radiation Fractions",         "Number of radiation therapy fractions administered.",                  "NUMBER", "1"),
    list("ER_STATUS",                "er_status",                          "ER Status",                   "Estrogen receptor (ER) status.",                                        "STRING", "1"),
    list("ER_PERCENT_POSITIVE",      "er_percent_positive",                "ER Percent Positive",         "Percentage of tumor cells with positive estrogen receptor staining.",  "NUMBER", "1"),
    list("PR_STATUS",                "pr_status",                          "PR Status",                   "Progesterone receptor (PR) status.",                                    "STRING", "1"),
    list("PR_PERCENT_POSITIVE",      "pr_percent_positive",                "PR Percent Positive",         "Percentage of tumor cells with positive progesterone receptor staining.", "NUMBER", "1"),
    list("HER2_IHC_STATUS",          "her2_ihc_status",                    "HER2 IHC Status",             "HER2 status by immunohistochemistry.",                                  "STRING", "1"),
    list("HER2_ISH_STATUS",          "her2_ish_status",                    "HER2 ISH Status",             "HER2 status by in-situ hybridization.",                                 "STRING", "1"),
    list("HPV_IHC_STATUS",           "hpv_ihc_status",                     "HPV IHC Status",              "Human papillomavirus status by immunohistochemistry.",                 "STRING", "1"),
    list("HPV_PCR_STATUS",           "hpv_pcr_status",                     "HPV PCR Status",              "Human papillomavirus status by PCR.",                                   "STRING", "1"),
    list("HPV_STRAIN",               "hpv_strain",                         "HPV Strain",                  "Human papillomavirus strain.",                                          "STRING", "1"),
    list("CA125",                    "ca125",                              "CA-125",                      "CA-125 tumor marker level.",                                            "NUMBER", "1"),
    list("CEA",                      "cea",                                "CEA",                         "Carcinoembryonic antigen (CEA) level.",                                 "NUMBER", "1"),
    list("PSA_LEVEL",                "psa_level",                          "PSA Level",                   "Prostate-specific antigen (PSA) level.",                               "NUMBER", "1")
  )
  write_cbio_table(m, col_defs, opt$output)
}
