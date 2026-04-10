#!/usr/bin/env Rscript

# Convert PURPLE CNV output to cBioPortal SEG and discrete CNA DISCRETE_LONG formats.
# PURPLE produces purity/ploidy-adjusted absolute copy numbers.
#
# Discrete CNA thresholds on minCopyNumber (gene-level minimum):
#   < 0.5  → -2  (homozygous deletion)
#   < 1.5  → -1  (hemizygous deletion)
#   <= 3.5 →  0  (neutral / diploid)
#   <= 6.0 →  1  (gain)
#   > 6.0  →  2  (high-level amplification)
#
# SEG output uses log2(copyNumber / 2) as the segment mean value.

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
})

option_list <- list(
    make_option("--purple_cnv_somatic",  type = "character", help = "PURPLE somatic CNV segment file (<sample>.purple.cnv.somatic.tsv)"),
    make_option("--purple_cnv_gene",     type = "character", help = "PURPLE gene-level CNV file (<sample>.purple.cnv.gene.tsv)"),
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

purple_gene[, cn_value := fcase(
    min_cn < 0.5,  -2L,
    min_cn < 1.5,  -1L,
    min_cn <= 3.5,  0L,
    min_cn <= 6.0,  1L,
    default =       2L
)]

# Merge with Entrez IDs
purple_annotated <- merge(purple_gene, gene_entrez, by.x = "Hugo_Symbol", by.y = "gene_symbol", all.x = TRUE)

discrete_long <- data.table(
    Hugo_Symbol    = purple_annotated$Hugo_Symbol,
    Entrez_Gene_Id = ifelse(is.na(purple_annotated$entrez_ncbi_id), "", as.character(purple_annotated$entrez_ncbi_id)),
    Sample_Id      = opt$sample_id,
    Value          = purple_annotated$cn_value
)

# One entry per gene per sample (keep first occurrence)
discrete_long <- unique(discrete_long, by = c("Hugo_Symbol", "Sample_Id"))

cat("Writing discrete CNA long file:", opt$output_long, "\n")
write.table(discrete_long, opt$output_long, sep = "\t", quote = FALSE, row.names = FALSE)

# ── 2. Segment-level copy numbers (SEG) ────────────────────────────────────

cat("Reading PURPLE somatic CNV segments:", opt$purple_cnv_somatic, "\n")
purple_somatic <- fread(opt$purple_cnv_somatic, header = TRUE)

# Normalise chromosome names: remove 'chr' prefix for cBioPortal compatibility
purple_somatic[, chrom_clean := sub("^chr", "", chromosome)]

seg_data <- data.table(
    ID        = opt$sample_id,
    chrom     = purple_somatic$chrom_clean,
    loc.start = purple_somatic$start,
    loc.end   = purple_somatic$end,
    num.mark  = purple_somatic$bafCount,
    seg.mean  = round(log2(pmax(purple_somatic$copyNumber, 0.001) / 2), 4)
)

cat("Writing SEG file:", opt$output_seg, "\n")
write.table(seg_data, opt$output_seg, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Done! Genes:", nrow(discrete_long), "| Segments:", nrow(seg_data), "\n")
