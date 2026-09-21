#!/usr/bin/env Rscript

# Convert PURPLE CNV output to cBioPortal SEG and discrete CNA DISCRETE_LONG formats.
# PURPLE produces purity/ploidy-adjusted absolute copy numbers.
#
# Copy numbers are scored against the copy number the sample is expected to carry,
# which is not 2 on the sex chromosomes (see expected_cn() below). Copy number is
# first converted to its diploid equivalent, min_cn * 2 / baseline, and then binned:
#   < 0.5  → -2  (homozygous deletion)
#   < 1.5  → -1  (hemizygous deletion)
#   <= 3.5 →  0  (neutral / diploid)
#   <= 6.0 →  1  (gain)
#   > 6.0  →  2  (high-level amplification)
# On a haploid baseline an absolute copy number below 0.5 is always -2: there is
# only one copy to lose, so the hemizygous bin does not apply.
#
# SEG output uses log2(copyNumber / baseline) as the segment mean value.

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

# ── 0. Sample sex, for the chrX/chrY baseline ──────────────────────────────
# PURPLE's own gender call: MALE, FEMALE or MALE_KLINEFELTER. Anything else
# (including no purity file at all) falls back to NA, which expected_cn()
# treats as diploid everywhere -- today's behaviour.

sample_sex <- NA_character_
if (!is.null(opt$purple_purity) && file.exists(opt$purple_purity)) {
    purity <- fread(opt$purple_purity, header = TRUE)
    if ("gender" %in% names(purity) && nrow(purity) > 0) {
        sample_sex <- toupper(trimws(purity$gender[1]))
        if (!sample_sex %in% c("MALE", "FEMALE", "MALE_KLINEFELTER")) sample_sex <- NA_character_
    }
}
cat("Resolved sample sex:", ifelse(is.na(sample_sex), "unknown (diploid baseline)", sample_sex), "\n")

# Expected copy number baseline for a chromosome given the sample's sex.
# NA means "not measured for this sex" (chrY in a female) -- the caller drops
# those rows rather than scoring them against a baseline that doesn't exist.
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
