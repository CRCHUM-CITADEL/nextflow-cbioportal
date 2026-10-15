#!/usr/bin/env Rscript

# Convert Oarfish gene expression output to cBioPortal expression format.

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
})

option_list <- list(
    make_option(c("-i", "--input"),       type = "character", help = "Expression TSV"),
    make_option(c("-m", "--gtf_gz"), 	  type = "character", help = "GTF.gz file for mapping"),
    make_option(c("-s", "--sample_id"),   type = "character", default = "Sample1", help = "Sample identifier [default=%default]"),
    make_option(c("-o", "--output"),      type = "character", default = "expression_data.txt", help = "Output file [default=%default]"),
    make_option("--min_tpm",             type = "double",    default = 0,
                help = "Minimum TPM threshold; 0 = no filter [default=%default]")
)

opt <- parse_args(OptionParser(option_list = option_list))


cat("Configuration:\n")
cat("  Input:           ", opt$input, "\n")
cat("  GTF file:        ", opt$gtf_gz, "\n")
cat("  Sample ID:       ", opt$sample_id, "\n")
cat("  Output:          ", opt$output, "\n")
cat("  Min TPM:         ", opt$min_tpm, "\n\n")

strip_version <- function(x) sub("\\.\\d+$", "", x)

extract_attr <- function(attr, key) {
  pattern <- paste0('.*', key, ' "([^"]+)".*')
  ifelse(grepl(paste0(key, ' "'), attr), sub(pattern, "\\1", attr), NA_character_)
}

# ── 1. Read expression ─────────────────────────────────────────────────

cat("Reading Isofox expression file...\n")
expr <- fread(opt$input, header = TRUE)
cat("  Records:", nrow(expr), "\n")

# calculate TPM
rpk <- expr$num_reads / expr$len
expr$TPM <- rpk / sum(rpk) * 1e6

expr_dt <- data.table(
    ensembl_id  = expr[["tname"]],
    tpm         = as.numeric(expr[["TPM"]])
)

# ── 2. Map Ensembl → Entrez ID ────────────────────────────────────────────────
con <- gzfile(opt$gtf_gz, "r")
pieces <- list()

repeat {
  lines <- readLines(con, n = 500000)
  if (length(lines) == 0) break
  lines <- lines[grepl("\ttranscript\t", lines, fixed = TRUE)]
  if (length(lines) == 0) next
  attr <- vapply(strsplit(lines, "\t", fixed = TRUE), `[`, "", 9)
  pieces[[length(pieces) + 1]] <- data.frame(
    tx_id     = strip_version(extract_attr(attr, "transcript_id")),
    gene_id   = strip_version(extract_attr(attr, "gene_id")),
    gene_name = extract_attr(attr, "gene_name"),
    stringsAsFactors = FALSE
  )
}

close(con)
tx2gene <- unique(do.call(rbind, pieces))
message("Read ", nrow(tx2gene), " transcript-to-gene mappings from the GTF.")

expr_dt$ensembl_id <- sub("\\.[0-9]+$", "", expr_dt$ensembl_id)

expr_dt_gene <- merge(expr_dt, tx2gene, by.x = "ensembl_id", by.y = "tx_id", all.x = TRUE, sort = FALSE)

# ── 3. Format and deduplicate ─────────────────────────────────────────────────

formatted <- expr_dt_gene[, .(
    Hugo_Symbol      = as.character(ifelse(is.na(gene_name), "", gene_name)),
    Entrez_Gene_Id   = as.character(ifelse(is.na(gene_id), "", gene_id)), 
    tpm_col          = tpm
)]
setnames(formatted, "tpm_col", opt$sample_id)

# Resolve duplicate gene symbols: keep highest TPM
cat("Resolving duplicate gene symbols...\n")
n_dups <- sum(duplicated(formatted$Hugo_Symbol))
cat("  Duplicate symbols:", n_dups, "\n")
formatted <- formatted[order(-get(opt$sample_id))]
formatted <- unique(formatted, by = "Hugo_Symbol")

# ── 4. Apply minimum TPM threshold ───────────────────────────────────────────

if (opt$min_tpm > 0) {
    before <- nrow(formatted)
    formatted <- formatted[get(opt$sample_id) >= opt$min_tpm]
    cat("  Filtered", before - nrow(formatted), "genes below TPM threshold", opt$min_tpm, "\n")
}

# ── 5. Write output ───────────────────────────────────────────────────────────

cat("Writing expression file:", opt$output, "\n")
write.table(formatted, opt$output, sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Total genes in output:", nrow(formatted), "\n")
cat("Done!\n")

