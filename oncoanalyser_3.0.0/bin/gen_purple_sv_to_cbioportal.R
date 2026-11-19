#!/usr/bin/env Rscript

# Convert PURPLE's somatic SV VCF (<tumor>.purple.sv.vcf.gz) to cBioPortal data_sv.txt.
#
# PURPLE passes ESVEE's records through (filtered and unfiltered), adds inferred single breakends
# (FILTER=INFERRED) and tumour-in-normal recovered SVs (INFO TINC_RECOVERED), and annotates each
# breakend with PURPLE_AF / PURPLE_CN / PURPLE_CN_CHANGE (lists: [this breakend, mate]) and PURPLE_JCN.
# This script:
#   1. Keeps PASS BND records that have a mate (drops INFERRED, PON-filtered and single breakends).
#   2. Reads the tumour's VF / SF / DF from the sample column named --tumor_id.
#   3. Pairs BNDs by MATEID, annotates each site with the overlapping gene (both sites must hit one).
#   4. Classifies DELETION/DUPLICATION/INVERSION/TRANSLOCATION from the BND orientations.
#   5. Writes cBioPortal SV format plus HMF.PURPLE_* namespace columns (SV_Status = SOMATIC).

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
})

option_list <- list(
    make_option(c("-i", "--input"),       type = "character", help = "PURPLE somatic SV VCF (.purple.sv.vcf.gz or .vcf)"),
    make_option(c("-s", "--sample"),      type = "character", help = "Sample ID"),
    make_option(c("-t", "--tumor_id"),    type = "character", help = "Tumour sample column name in the VCF (e.g. <subject>-T)"),
    make_option("--ensembl_annotations",  type = "character", help = "Ensembl annotations TSV (cols: ensembl_id, entrez_ncbi_id, gene_symbol, chr, start, end, ...)"),
    make_option(c("-o", "--output"),      type = "character", help = "Output data_sv.txt path")
)

opt <- parse_args(OptionParser(option_list = option_list))

for (arg in c("input", "sample", "tumor_id", "ensembl_annotations", "output")) {
    if (is.null(opt[[arg]])) stop(paste("Missing required argument: --", arg, sep = ""))
}

# Column order of data_sv.txt. HMF.* columns are always written so every SV file has the same width
# (gen_isofox_fusion_to_cbioportal.R writes the same set as NA).
SV_COLUMNS <- c(
    "Sample_Id", "SV_Status",
    "Site1_Hugo_Symbol", "Site1_Ensembl_Transcript_Id", "Site1_Region_Number", "Site1_Region",
    "Site2_Hugo_Symbol", "Site2_Ensembl_Transcript_Id", "Site2_Region_Number", "Site2_Region",
    "Site2_Effect_On_Frame", "NCBI_Build", "Class", "DNA_Support", "RNA_Support",
    "Tumor_Variant_Count", "Connection_Type", "Breakpoint_Type", "Event_Info", "Annotation",
    "Site1_Chromosome", "Site1_Position", "Site2_Chromosome", "Site2_Position",
    "Tumor_Split_Read_Count", "Tumor_Paired_End_Read_Count",
    "HMF.PURPLE_JCN",
    "HMF.PURPLE_AF_SITE1", "HMF.PURPLE_AF_SITE2",
    "HMF.PURPLE_CN_SITE1", "HMF.PURPLE_CN_SITE2",
    "HMF.PURPLE_CN_CHANGE_SITE1", "HMF.PURPLE_CN_CHANGE_SITE2"
)

# ── BND helpers ───────────────────────────────────────────────────────────────

# REF before the bracket ("N[chr:pos[", "N]chr:pos]") → + strand; after it → - strand.
# perl = TRUE: in R's default POSIX regex "\\[" inside a bracket is a literal backslash, so it never matched.
get_bnd_strand <- function(alt) {
    if (grepl("^[A-Za-z.]+[\\[\\]]", alt, perl = TRUE)) "+" else "-"
}

# Same chromosome: (+,-) DELETION, (-,+) DUPLICATION, (+,+)/(-,-) INVERSION; else TRANSLOCATION.
classify_sv <- function(chr1, chr2, strand1, strand2) {
    if (chr1 != chr2) return("TRANSLOCATION")
    if (strand1 == "+" && strand2 == "-") return("DELETION")
    if (strand1 == "-" && strand2 == "+") return("DUPLICATION")
    return("INVERSION")
}

# + strand → the gene contributes its 5' end; - strand → its 3' end.
get_connection_type <- function(strand1, strand2) {
    s1 <- if (strand1 == "+") "5" else "3"
    s2 <- if (strand2 == "+") "5" else "3"
    paste0(s1, "to", s2)
}

# ── Helpers ──────────────────────────────────────────────────────────────────

parse_info <- function(info_str) {
    parts <- strsplit(info_str, ";", fixed = TRUE)[[1]]
    kv <- lapply(parts, function(p) {
        eq <- regexpr("=", p, fixed = TRUE)
        if (eq == -1) return(list(key = p, val = TRUE))
        list(key = substr(p, 1, eq - 1), val = substr(p, eq + 1, nchar(p)))
    })
    vals <- lapply(kv, `[[`, "val")
    names(vals) <- sapply(kv, `[[`, "key")
    vals
}

parse_format <- function(fmt_str, smp_str) {
    keys <- strsplit(fmt_str, ":", fixed = TRUE)[[1]]
    vals <- strsplit(smp_str, ":", fixed = TRUE)[[1]]
    if (length(vals) < length(keys)) vals <- c(vals, rep(NA_character_, length(keys) - length(vals)))
    setNames(as.list(vals[seq_along(keys)]), keys)
}

# Integer FORMAT field, NA when absent or ".".
format_int <- function(fmt_data, key) {
    v <- fmt_data[[key]]
    if (is.null(v) || is.na(v) || v == ".") return(NA_integer_)
    suppressWarnings(as.integer(v))
}

# First element of a PURPLE INFO list: PURPLE writes [this breakend, mate] on each record.
info_first <- function(info, key) {
    v <- info[[key]]
    if (is.null(v) || isTRUE(v)) return(NA_real_)
    suppressWarnings(as.numeric(strsplit(v, ",", fixed = TRUE)[[1]][1]))
}

# ── Read VCF ─────────────────────────────────────────────────────────────────

cat("Reading PURPLE SV VCF:", opt$input, "\n")
con <- if (grepl("\\.gz$", opt$input)) gzfile(opt$input, "rt") else file(opt$input, "rt")
vcf_lines <- readLines(con)
close(con)

# Tumour column by name: ESVEE's sample column order is not documented.
header_line <- grep("^#CHROM", vcf_lines, value = TRUE)
if (length(header_line) == 0) stop("VCF header (#CHROM) not found")
col_names <- strsplit(header_line, "\t", fixed = TRUE)[[1]]
sample_cols <- if (length(col_names) > 9) col_names[10:length(col_names)] else character(0)
tumor_col_idx <- match(opt$tumor_id, col_names)
if (is.na(tumor_col_idx)) {
    if (length(sample_cols) != 1) {
        stop(sprintf("Tumour column '%s' not found in VCF sample columns: %s",
                     opt$tumor_id, paste(sample_cols, collapse = ", ")))
    }
    tumor_col_idx <- 10L
    cat("Tumour column", opt$tumor_id, "not found; using the only sample column:", sample_cols, "\n")
}
cat("Tumour column:", col_names[tumor_col_idx], "\n")

# ── Parse BND records ─────────────────────────────────────────────────────────

data_lines <- vcf_lines[!grepl("^#", vcf_lines)]
cat("Total variant records:", length(data_lines), "\n")

bnd_list <- list()

for (line in data_lines) {
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(fields) < tumor_col_idx) next
    if (fields[7] != "PASS") next

    info <- parse_info(fields[8])
    if (!identical(info[["SVTYPE"]], "BND")) next
    mateid <- info[["MATEID"]]
    if (is.null(mateid) || isTRUE(mateid)) next

    fmt_data <- parse_format(fields[9], fields[tumor_col_idx])

    bnd_list[[fields[3]]] <- list(
        id             = fields[3],
        chrom          = fields[1],
        pos            = suppressWarnings(as.integer(fields[2])),
        strand         = get_bnd_strand(fields[5]),
        mateid         = as.character(mateid),
        total_frags    = format_int(fmt_data, "VF"),
        split_frags    = format_int(fmt_data, "SF"),
        disc_frags     = format_int(fmt_data, "DF"),
        purple_jcn     = info_first(info, "PURPLE_JCN"),
        purple_af      = info_first(info, "PURPLE_AF"),
        purple_cn      = info_first(info, "PURPLE_CN"),
        purple_cn_chg  = info_first(info, "PURPLE_CN_CHANGE")
    )
}

cat("PASS BND records:", length(bnd_list), "\n")

# ── Load gene annotations for interval overlap ────────────────────────────────

cat("Reading Ensembl annotations:", opt$ensembl_annotations, "\n")
ann <- fread(opt$ensembl_annotations, header = TRUE)
chr_col <- if ("chr" %in% names(ann)) "chr" else names(ann)[4]
ann[, chrom_key := sub("^chr", "", get(chr_col))]
ann_genes <- ann[, .(chrom_key, start, end, gene_symbol)]
setkey(ann_genes, chrom_key, start, end)

find_gene <- function(chrom, pos) {
    q <- data.table(chrom_key = sub("^chr", "", chrom), start = pos, end = pos)
    hits <- foverlaps(q, ann_genes, by.x = c("chrom_key", "start", "end"),
                      by.y = c("chrom_key", "start", "end"), type = "within", nomatch = NULL)
    if (nrow(hits) > 0) hits$gene_symbol[1] else NA_character_
}

# ── Pair BND records and build SV table ──────────────────────────────────────

seen_pairs <- new.env(hash = TRUE)
sv_rows <- list()

for (id in names(bnd_list)) {
    rec  <- bnd_list[[id]]
    mate <- bnd_list[[rec$mateid]]
    if (is.null(mate)) next

    pair_key <- paste(sort(c(id, rec$mateid)), collapse = "_")
    if (!is.null(seen_pairs[[pair_key]])) next
    seen_pairs[[pair_key]] <- TRUE

    gene1 <- tryCatch(find_gene(rec$chrom, rec$pos),  error = function(e) NA_character_)
    gene2 <- tryCatch(find_gene(mate$chrom, mate$pos), error = function(e) NA_character_)

    # cBioPortal needs a Hugo symbol at both sites
    if (is.na(gene1) || is.na(gene2) || gene1 == "" || gene2 == "") next

    chr1 <- sub("^chr", "", rec$chrom)
    chr2 <- sub("^chr", "", mate$chrom)

    sv_rows[[length(sv_rows) + 1]] <- data.table(
        Sample_Id                   = opt$sample,
        SV_Status                   = "SOMATIC",
        Site1_Hugo_Symbol           = gene1,
        Site1_Ensembl_Transcript_Id = NA_character_,
        Site1_Region_Number         = NA_character_,
        Site1_Region                = NA_character_,
        Site2_Hugo_Symbol           = gene2,
        Site2_Ensembl_Transcript_Id = NA_character_,
        Site2_Region_Number         = NA_character_,
        Site2_Region                = NA_character_,
        Site2_Effect_On_Frame       = NA_character_,
        NCBI_Build                  = "GRCh38",
        Class                       = classify_sv(chr1, chr2, rec$strand, mate$strand),
        DNA_Support                 = "Yes",
        RNA_Support                 = "No",
        Tumor_Variant_Count         = rec$total_frags,
        Connection_Type             = get_connection_type(rec$strand, mate$strand),
        Breakpoint_Type             = "PRECISE",
        Event_Info                  = paste0("DNA WGS Fusion: ", gene1, "--", gene2),
        Annotation                  = NA_character_,
        Site1_Chromosome            = chr1,
        Site1_Position              = rec$pos,
        Site2_Chromosome            = chr2,
        Site2_Position              = mate$pos,
        Tumor_Split_Read_Count      = rec$split_frags,
        Tumor_Paired_End_Read_Count = rec$disc_frags,
        HMF.PURPLE_JCN              = rec$purple_jcn,
        HMF.PURPLE_AF_SITE1         = rec$purple_af,
        HMF.PURPLE_AF_SITE2         = mate$purple_af,
        HMF.PURPLE_CN_SITE1         = rec$purple_cn,
        HMF.PURPLE_CN_SITE2         = mate$purple_cn,
        HMF.PURPLE_CN_CHANGE_SITE1  = rec$purple_cn_chg,
        HMF.PURPLE_CN_CHANGE_SITE2  = mate$purple_cn_chg
    )
}

if (length(sv_rows) > 0) {
    result <- rbindlist(sv_rows, fill = TRUE)
} else {
    cat("No PASS BND pairs with genes at both sites. Writing header only.\n")
    result <- as.data.table(setNames(replicate(length(SV_COLUMNS), character(), simplify = FALSE), SV_COLUMNS))
}
setcolorder(result, SV_COLUMNS)

# Same breakpoint pair called more than once: keep the call with the most variant fragments.
result <- result[order(-nafill(as.integer(Tumor_Variant_Count), fill = 0L))]
result <- unique(result, by = c("Site1_Chromosome", "Site1_Position",
                                "Site2_Chromosome", "Site2_Position"))

cat("Writing SV output:", opt$output, "\n")
cat("Total SV records:", nrow(result), "\n")
write.table(result, opt$output, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
cat("Done!\n")
