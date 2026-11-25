#!/usr/bin/env Rscript

# PURPLE CNV → cBioPortal SEG (log2(CN / baseline)) and DISCRETE_LONG (CN scaled to diploid, then binned; bins in the fcase() below)

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
})

option_list <- list(
    make_option("--purple_cnv_somatic",  type = "character", help = "PURPLE somatic CNV segment file (<sample>.purple.cnv.somatic.tsv)"),
    make_option("--purple_cnv_gene",     type = "character", help = "PURPLE gene-level CNV file (<sample>.purple.cnv.gene.tsv)"),
    make_option("--purple_purity",       type = "character", default = NULL, help = "PURPLE purity file (<sample>.purple.purity.tsv), for the gender-based chrX/chrY baseline. Optional: falls back to a diploid baseline when absent"),
    make_option("--sample_id",           type = "character", help = "Tumor sample ID (used as column name in output)"),
    make_option("--ensembl_annotations", type = "character", help = "Ensembl annotations TSV (cols: ensembl_id, entrez_ncbi_id, gene_symbol, ...)"),
    make_option("--output_seg",          type = "character", help = "Output SEG file path"),
    make_option("--output_long",         type = "character", help = "Output discrete CNA DISCRETE_LONG file path")
)

opt <- parse_args(OptionParser(option_list = option_list))

for (arg in c("purple_cnv_somatic", "purple_cnv_gene", "sample_id",
              "ensembl_annotations", "output_seg", "output_long")) {
    if (is.null(opt[[arg]])) stop(paste("Missing required argument: --", arg, sep = ""))
}

# ── 0. Sample sex for the chrX/chrY baseline: PURPLE's gender call, else NA (diploid everywhere) ──

sample_sex <- NA_character_
if (!is.null(opt$purple_purity) && file.exists(opt$purple_purity)) {
    purity <- fread(opt$purple_purity, header = TRUE)
    if ("gender" %in% names(purity) && nrow(purity) > 0) {
        sample_sex <- toupper(trimws(purity$gender[1]))
        if (!sample_sex %in% c("MALE", "FEMALE", "MALE_KLINEFELTER")) sample_sex <- NA_character_
    }
}
cat("Resolved sample sex:", ifelse(is.na(sample_sex), "unknown (diploid baseline)", sample_sex), "\n")

# Expected copy number for a chromosome given sex; NA (chrY in a female) means the row is dropped
expected_cn <- function(chromosome, sex) {
    chrom <- sub("^chr", "", chromosome)
    fcase(
        chrom == "X" & identical(sex, "MALE"),  1,
        chrom == "X",                            2,
        chrom == "Y" & identical(sex, "FEMALE"), NA_real_,
        chrom == "Y" & is.na(sex),                2,
        chrom == "Y",                            1,
        default =                                2
    )
}

# ── 1. Gene-level copy numbers (DISCRETE_LONG) ─────────────────────────────

cat("Reading PURPLE gene-level CNV:", opt$purple_cnv_gene, "\n")
purple_gene <- fread(opt$purple_cnv_gene, header = TRUE)
print(head(purple_gene))
cat("Reading Ensembl annotations:", opt$ensembl_annotations, "\n")
annotations <- fread(opt$ensembl_annotations, header = TRUE)
# Build gene_symbol → entrez_ncbi_id lookup (de-duplicated)
gene_entrez <- unique(annotations[!is.na(gene_symbol), .(gene_symbol, entrez_ncbi_id)])
setkey(gene_entrez, gene_symbol)

# Normalise column names across PURPLE versions
gene_col   <- if ("gene"          %in% names(purple_gene)) "gene"          else "Gene"
min_cn_col <- if ("minCopyNumber" %in% names(purple_gene)) "minCopyNumber" else "minCN"

purple_gene[, Hugo_Symbol := get(gene_col)]
purple_gene[, min_cn      := as.numeric(get(min_cn_col))]
purple_gene[, baseline    := expected_cn(chromosome, sample_sex)]

# chrY genes in a female sample were never measured -- drop, don't score.
purple_gene <- purple_gene[!is.na(baseline)]

purple_gene[, cn_equiv := min_cn * 2 / baseline]

# An unparseable copy number is dropped, not scored (fcase() would default it to 2) nor written as NA (aborts the import)
n_unscored <- sum(is.na(purple_gene$cn_equiv))
if (n_unscored > 0) {
    warning(sprintf("Dropping %d gene(s) with an unparseable %s",
                    n_unscored, min_cn_col), call. = FALSE)
    purple_gene <- purple_gene[!is.na(cn_equiv)]
}

purple_gene[, cn_value := fcase(
    baseline < 2 & min_cn < 0.5, -2L,   # haploid: the only copy is gone
    cn_equiv <  0.5,             -2L,
    cn_equiv <  1.5,             -1L,
    cn_equiv <= 3.5,              0L,
    cn_equiv <= 6.0,              1L,
    default  =                    2L
)]

# Merge with Entrez IDs
purple_annotated <- merge(purple_gene, gene_entrez, by.x = "Hugo_Symbol", by.y = "gene_symbol", all.x = TRUE)

discrete_long <- data.table(
    Hugo_Symbol    = purple_annotated$Hugo_Symbol,
    Entrez_Gene_Id = purple_annotated$entrez_ncbi_id,
    Sample_Id      = opt$sample_id,
    Value          = purple_annotated$cn_value
)

# Remove rows with NA Entrez_Gene_Id, then deduplicate keeping first occurrence
discrete_long <- discrete_long[!is.na(Entrez_Gene_Id)]
discrete_long <- discrete_long[!duplicated(discrete_long, by = c("Sample_Id", "Entrez_Gene_Id"))]

cat("Writing discrete CNA long file:", opt$output_long, "\n")
write.table(discrete_long, opt$output_long, sep = "\t", quote = FALSE, row.names = FALSE)

# ── 2. Segment-level copy numbers (SEG) ────────────────────────────────────

cat("Reading PURPLE somatic CNV segments:", opt$purple_cnv_somatic, "\n")
purple_somatic <- fread(opt$purple_cnv_somatic, header = TRUE)

# Normalise chromosome names: remove 'chr' prefix for cBioPortal compatibility
purple_somatic[, chrom_clean := sub("^chr", "", chromosome)]
purple_somatic[, baseline    := expected_cn(chromosome, sample_sex)]

# chrY segments in a female sample were never measured -- drop, don't score.
purple_somatic <- purple_somatic[!is.na(baseline)]

seg_data <- data.table(
    ID        = opt$sample_id,
    chrom     = purple_somatic$chrom_clean,
    loc.start = purple_somatic$start,
    loc.end   = purple_somatic$end,
    num.mark  = purple_somatic$bafCount,
    seg.mean  = round(log2(pmax(purple_somatic$copyNumber, 0.001) / purple_somatic$baseline), 4)
)

cat("Writing SEG file:", opt$output_seg, "\n")
write.table(seg_data, opt$output_seg, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Done! Genes:", nrow(discrete_long), "| Segments:", nrow(seg_data), "\n")
