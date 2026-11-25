# Shared helpers for the clinical scripts; staged as a path() input so it is part of the task hash

# day_interval from an ARGO JSON date string; diagnosis is day 0 (negative = before)
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

# MOHCCN map CSV → named vector (value → code): first 3 columns, empty rows skipped, both sides trimmed
read_mohccn_map <- function(filepath) {
  if (is.null(filepath) || !file.exists(filepath)) return(NULL)
  raw <- read.csv(filepath, header=FALSE, skip=1, stringsAsFactors=FALSE, fill=TRUE)
  raw <- raw[, 1:3, drop=FALSE]
  colnames(raw) <- c("full_str", "value", "code")
  raw <- raw[nchar(trimws(raw$value)) > 0, ]
  setNames(trimws(raw$code), trimws(raw$value))
}

# Reverse a read_mohccn_map() map into code → label
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

# ICD-O codes → labels at the C+2-digit parent ("C22.0"/"C220" → "C22"); code lists join with "; ", unmatched dropped
apply_ps_label <- function(codes, reversed_map) {
  if (is.null(reversed_map)) return(rep(NA_character_, length(codes)))
  vapply(as.character(codes), function(value) {
    if (is.na(value)) return(NA_character_)
    prefixes <- sub("^(C\\d{2}).*", "\\1", trimws(strsplit(value, ",", fixed = TRUE)[[1]]))
    labels   <- unname(reversed_map[prefixes])
    labels   <- labels[!is.na(labels) & labels != ""]
    if (length(labels) == 0) NA_character_ else paste(labels, collapse = "; ")
  }, character(1), USE.NAMES = FALSE)
}

# ── Sample linking: only Total DNA solid-tissue tumours with an -<n>[A-Z]*[DR]T suffix become sample rows ──
is_cbio_sample <- function(reg) {
  reg$sample_type               == "Total DNA" &
  reg$tumour_normal_designation == "Tumour" &
  reg$specimen_tissue_source    == "Solid tissue" &
  grepl("-\\d+[A-Z]*[DR]T$", reg$submitter_sample_id)
}

# Genomic linking file (subject_id TAB sample_id): filters to genomic subjects and gives their sample IDs
read_genomic_subjects <- function(filepath) {
  read.table(filepath, header=TRUE, sep="\t", stringsAsFactors=FALSE)
}

# Bind data frames with different columns, filling missing ones with NA
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
