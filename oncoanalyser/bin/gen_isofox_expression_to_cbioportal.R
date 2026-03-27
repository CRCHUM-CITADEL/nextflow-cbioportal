#!/usr/bin/env Rscript

# Convert Isofox gene expression output to cBioPortal expression format.
#
# Isofox exp.tsv already contains gene names (GeneName) and Ensembl IDs (GeneId).
# This script:
#   1. Reads Isofox exp.tsv.
#   2. Maps Ensembl GeneId → Entrez ID using the annotation file.
#   3. Deduplicates gene symbols (keeps highest TPM).
#   4. Outputs cBioPortal format: Hugo_Symbol, Entrez_Gene_Id, <sample_id>.

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
})

option_list <- list(
    make_option(c("-i", "--input"),       type = "character", help = "Isofox expression TSV (<sample>.isofox.exp.tsv)"),
    make_option(c("-g", "--gene_map"),    type = "character", help = "Ensembl annotations TSV (cols: ensembl_id, entrez_ncbi_id, gene_symbol, ...)"),
    make_option(c("-s", "--sample_id"),   type = "character", default = "Sample1", help = "Sample identifier [default=%default]"),
    make_option(c("-o", "--output"),      type = "character", default = "expression_data.txt", help = "Output file [default=%default]"),
    make_option("--use_adjusted_tpm",     action = "store_true", default = FALSE,
                help = "Use AdjustedTPM instead of TPM (excludes highly abundant non-coding genes) [default=FALSE]"),
    make_option("--min_tpm",             type = "double",    default = 0,
                help = "Minimum TPM threshold; 0 = no filter [default=%default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

for (arg in c("input", "gene_map")) {
    if (is.null(opt[[arg]])) stop(paste("Missing required argument: --", arg, sep = ""))
}

cat("Configuration:\n")
cat("  Input:           ", opt$input, "\n")
cat("  Gene map:        ", opt$gene_map, "\n")
cat("  Sample ID:       ", opt$sample_id, "\n")
cat("  Output:          ", opt$output, "\n")
cat("  Use AdjustedTPM:", opt$use_adjusted_tpm, "\n")
cat("  Min TPM:         ", opt$min_tpm, "\n\n")

# ── 1. Read Isofox expression ─────────────────────────────────────────────────

cat("Reading Isofox expression file...\n")
expr <- fread(opt$input, header = TRUE)
cat("  Records:", nrow(expr), "\n")

cols <- names(expr)

# Auto-detect column names (Isofox column naming is stable but version can vary)
gene_id_col   <- if ("GeneId"      %in% cols) "GeneId"      else "geneId"
gene_name_col <- if ("GeneName"    %in% cols) "GeneName"    else "geneName"

tpm_col <- if (opt$use_adjusted_tpm && "AdjustedTPM" %in% cols) {
    cat("  Using AdjustedTPM\n")
    "AdjustedTPM"
} else if ("TPM" %in% cols) {
    "TPM"
} else {
    # Fallback: look for any column containing 'TPM' (case-insensitive)
    hit <- grep("tpm", cols, ignore.case = TRUE, value = TRUE)[1]
    if (is.na(hit)) stop("Cannot find a TPM column in Isofox expression file")
    cat("  Using column:", hit, "\n")
    hit
}

expr_dt <- data.table(
    ensembl_id  = expr[[gene_id_col]],
    Hugo_Symbol = expr[[gene_name_col]],
    tpm         = as.numeric(expr[[tpm_col]])
)

# ── 2. Map Ensembl → Entrez ID ────────────────────────────────────────────────

cat("Reading gene annotation file...\n")
gene_map <- fread(opt$gene_map, header = TRUE)
# Expected columns (Ensembl BioMart format):
#   ensembl_id, entrez_ncbi_id, gene_symbol, chr, start, end, strand, description
entrez_map <- unique(gene_map[, .(ensembl_id, entrez_ncbi_id)])
setkey(entrez_map, ensembl_id)

# Trim version numbers from Ensembl IDs (e.g. ENSG00000001234.5 → ENSG00000001234)
expr_dt[, ensembl_id_trim := sub("\\.[0-9]+$", "", ensembl_id)]

expr_ann <- merge(expr_dt, entrez_map, by.x = "ensembl_id_trim", by.y = "ensembl_id", all.x = TRUE)
cat("  Mapped", sum(!is.na(expr_ann$entrez_ncbi_id)), "of", nrow(expr_ann), "genes to Entrez IDs\n")

# ── 3. Format and deduplicate ─────────────────────────────────────────────────

formatted <- expr_ann[, .(
    Hugo_Symbol    = Hugo_Symbol,
    Entrez_Gene_Id = as.character(ifelse(is.na(entrez_ncbi_id), "", entrez_ncbi_id)),
    tpm_col        = tpm
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
