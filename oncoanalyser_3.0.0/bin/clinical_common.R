# clinical_common.R — shared helpers for the clinical output-generation scripts
# (build_clinical_table.R, write_clinical_table.R, gen_timeline_*.R, merge_timeline.R).
#
# Extracted from clin_format.R / gen_timeline.R, which carried verbatim or
# near-verbatim copies of every function below. Keeping one copy means the
# eligibility rule, the MOHCCN map format, and the date-parsing regexes stay
# in sync by construction instead of by discipline — see CLAUDE.md.
#
# Staged into each process as an explicit `path()` input (not relied upon via
# $projectDir/bin on PATH), so it participates in the task hash and -resume
# correctly invalidates when it changes.

# Extract day_interval from ARGO JSON date strings: {"day_interval": N, "month_interval": M}
# Diagnosis is day 0; negative values precede diagnosis, positive values follow it.
parse_day_interval <- function(x) {
  x <- as.character(x)
  m <- regmatches(x, regexpr('"day_interval"\\s*:\\s*(-?[0-9]+)', x))
  if (length(m) == 0 || nchar(m) == 0) return(NA_real_)
  as.numeric(sub('"day_interval"\\s*:\\s*(-?[0-9]+)', '\\1', m))
}

# Flatten a JSON array string to a comma-separated string: ["Surgery"] → Surgery
clean_json_array <- function(x) {
  x <- gsub('\\[|\\]', '', x)
  x <- gsub('"', '', x)
  trimws(x)
}

# Read MOHCCN mapping CSV → named vector: plain-text value → ontology/ICD-O code.
# Skips empty rows (primary_site.csv has 900+ trailing empty rows).
# Reads only first 3 cols (treatment_intent.csv has 26 cols due to trailing commas).
#
# Both value and code are trimmed. Before this shared library existed,
# clin_format.R built this map WITHOUT trimming while gen_timeline.R trimmed
# both sides — a whitespace-padded map entry would resolve in the timeline but
# silently become NA in the clinical file. Trimming here fixes that divergence
# for both callers.
read_mohccn_map <- function(filepath) {
  if (is.null(filepath) || !file.exists(filepath)) return(NULL)
  raw <- read.csv(filepath, header=FALSE, skip=1, stringsAsFactors=FALSE, fill=TRUE)
  raw <- raw[, 1:3, drop=FALSE]
  colnames(raw) <- c("full_str", "value", "code")
  raw <- raw[nchar(trimws(raw$value)) > 0, ]
  setNames(trimws(raw$code), trimws(raw$value))
}

# Reverse a label->code map (as returned by read_mohccn_map) into code->label,
# for looking up a human-readable label from an ICD-O/ontology code already
# present in the data (e.g. specimen_anatomic_location, relapse site).
reverse_mohccn_map <- function(map) {
  if (is.null(map)) return(NULL)
  setNames(names(map), map)
}

# Look up values in a MOHCCN map; returns NA_character_ where there is no match or empty code.
apply_mohccn_map <- function(values, map) {
  if (is.null(map)) return(rep(NA_character_, length(values)))
  result <- map[as.character(values)]
  result[is.na(result) | result == ""] <- NA_character_
  names(result) <- NULL
  result
}

# Apply a reversed primary-site-style map: ICD-O code → human-readable label.
# Extracts the parent code (C + 2 digits) before lookup, handling both
# decimal format ("C22.0" → "C22") and concatenated MOHCCN format ("C220" → "C22").
apply_ps_label <- function(codes, reversed_map) {
  if (is.null(reversed_map)) return(rep(NA_character_, length(codes)))
  prefixes <- sub("^(C\\d{2}).*", "\\1", trimws(as.character(codes)))
  result   <- reversed_map[prefixes]
  result[is.na(result) | result == ""] <- NA_character_
  unname(result)
}

# ── Sample linking ────────────────────────────────────────────────────────────
# A registration is a usable cBioPortal sample when it is a Total DNA solid-tissue
# tumour whose sample ID carries the MoHQ analyte/designation suffix (-1DT/-2FRT);
# other suffixes are registry bookkeeping rows, not sequenced samples. Germline
# normals (buffy coat and any other tissue source) are patient-level context only
# and never become rows in data_clinical_sample.txt.
is_cbio_sample <- function(reg) {
  reg$sample_type               == "Total DNA" &
  reg$tumour_normal_designation == "Tumour" &
  reg$specimen_tissue_source    == "Solid tissue" &
  grepl("-\\d+[A-Z]*[DR]T$", reg$submitter_sample_id)
}

# Read the genomic linking file (subject_id TAB sample_id) produced by the
# genomic workflow, used to filter clinical output to genomic subjects and
# adopt the genomic sample_id convention (see main.nf "both" mode wiring).
read_genomic_subjects <- function(filepath) {
  read.table(filepath, header=TRUE, sep="\t", stringsAsFactors=FALSE)
}

# Bind data frames with different columns; missing columns filled with NA.
# Used to union the per-EVENT_TYPE timeline parts into one data_timeline.txt.
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
