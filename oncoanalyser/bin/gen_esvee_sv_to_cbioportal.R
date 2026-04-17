#!/usr/bin/env Rscript

# Convert ESVEE somatic structural variant VCF to cBioPortal data_sv.txt format.
#
# ESVEE represents each SV as a pair of BND (breakend) records linked by MATEID.
# This script:
#   1. Reads the ESVEE somatic VCF (may be gzipped).
#   2. Filters to PASS BND records.
#   3. Pairs BNDs by MATEID to recover the two-end breakpoint.
#   4. Annotates each breakpoint with the overlapping gene from the Ensembl
#      annotation file (interval overlap using data.table::foverlaps).
#   5. Writes cBioPortal SV format.

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
})

option_list <- list(
    make_option(c("-i", "--input"),       type = "character", help = "Input ESVEE VCF (.esvee.unfiltered.vcf.gz or .vcf)"),
    make_option(c("-s", "--sample"),      type = "character", help = "Sample ID"),
    make_option("--ensembl_annotations",  type = "character", help = "Ensembl annotations TSV (cols: ensembl_id, entrez_ncbi_id, gene_symbol, chr, start, end, ...)"),
    make_option("--sv_status",            type = "character", default = "SOMATIC", help = "SV_Status value: SOMATIC or GERMLINE (default: SOMATIC)"),
    make_option(c("-o", "--output"),      type = "character", help = "Output data_sv.txt path")
)

opt <- parse_args(OptionParser(option_list = option_list))

for (arg in c("input", "sample", "ensembl_annotations", "output")) {
    if (is.null(opt[[arg]])) stop(paste("Missing required argument: --", arg, sep = ""))
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

extract_read_counts <- function(fmt_data) {
    # ESVEE FORMAT fields vary; try in priority order
    split_reads <- NA_integer_
    paired_reads <- NA_integer_

    # BVF = breakend variant fragments (split reads), VF = variant fragments (total)
    if ("BVF" %in% names(fmt_data) && !is.na(fmt_data$BVF)) split_reads  <- suppressWarnings(as.integer(fmt_data$BVF))
    if ("VF"  %in% names(fmt_data) && !is.na(fmt_data$VF))  paired_reads <- suppressWarnings(as.integer(fmt_data$VF))

    # Fall back to SR/PR (split reads / paired reads)
    if (is.na(split_reads)  && "SR" %in% names(fmt_data)) split_reads  <- suppressWarnings(as.integer(fmt_data$SR))
    if (is.na(paired_reads) && "PR" %in% names(fmt_data)) paired_reads <- suppressWarnings(as.integer(fmt_data$PR))

    # Last resort: second allele from AD (ref, alt)
    if (is.na(split_reads) && "AD" %in% names(fmt_data)) {
        ad <- strsplit(as.character(fmt_data$AD), ",", fixed = TRUE)[[1]]
        if (length(ad) >= 2) split_reads <- suppressWarnings(as.integer(ad[2]))
    }

    list(split = split_reads, paired = paired_reads)
}

# ── Read VCF ─────────────────────────────────────────────────────────────────

cat("Reading ESVEE VCF:", opt$input, "\n")
if (grepl("\\.gz$", opt$input)) {
    con <- gzfile(opt$input, "rt")
    vcf_lines <- readLines(con)
    close(con)
} else {
    vcf_lines <- readLines(opt$input)
}

# Identify tumor sample column (last column in #CHROM line)
header_line <- grep("^#CHROM", vcf_lines, value = TRUE)
if (length(header_line) == 0) stop("VCF header (#CHROM) not found")
col_names <- strsplit(header_line, "\t", fixed = TRUE)[[1]]
tumor_col_idx <- length(col_names)  # last sample column = tumor

# ── Parse BND records ─────────────────────────────────────────────────────────

data_lines <- vcf_lines[!grepl("^#", vcf_lines)]
cat("Total variant records:", length(data_lines), "\n")

bnd_list <- list()

for (line in data_lines) {
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(fields) < 9) next
    if (fields[7] != "PASS") next            # FILTER must be PASS

    info <- parse_info(fields[8])
    if (!identical(info[["SVTYPE"]], "BND")) next   # BNDs only

    mateid <- info[["MATEID"]]
    if (is.null(mateid) || is.na(mateid)) next

    fmt_data <- if (length(fields) >= tumor_col_idx) {
        parse_format(fields[9], fields[tumor_col_idx])
    } else {
        parse_format(fields[9], fields[10])
    }

    rc <- extract_read_counts(fmt_data)

    bnd_list[[fields[3]]] <- list(
        id          = fields[3],
        chrom       = fields[1],
        pos         = suppressWarnings(as.integer(fields[2])),
        alt         = fields[5],
        mateid      = as.character(mateid),
        split_reads = rc$split,
        paired_reads= rc$paired
    )
}

cat("PASS BND records:", length(bnd_list), "\n")

# ── Load gene annotations for interval overlap ────────────────────────────────

cat("Reading Ensembl annotations:", opt$ensembl_annotations, "\n")
ann <- fread(opt$ensembl_annotations, header = TRUE)
# Normalise chromosome column
chr_col <- if ("chr" %in% names(ann)) "chr" else names(ann)[4]
ann[, chrom_key := sub("^chr", "", get(chr_col))]
ann_genes <- ann[, .(chrom_key, start, end, gene_symbol)]
setkey(ann_genes, chrom_key, start, end)

find_gene <- function(chrom, pos) {
    chrom_clean <- sub("^chr", "", chrom)
    q <- data.table(chrom_key = chrom_clean, start = pos, end = pos)
    hits <- foverlaps(q, ann_genes, by.x = c("chrom_key","start","end"),
                      by.y = c("chrom_key","start","end"), type = "within", nomatch = NULL)
    if (nrow(hits) > 0) hits$gene_symbol[1] else NA_character_
}

# ── Pair BND records and build SV table ──────────────────────────────────────

seen_pairs <- character(0)
sv_rows <- list()

for (id in names(bnd_list)) {
    rec  <- bnd_list[[id]]
    mate <- bnd_list[[rec$mateid]]
    if (is.null(mate)) next

    pair_key <- paste(sort(c(id, rec$mateid)), collapse = "_")
    if (pair_key %in% seen_pairs) next
    seen_pairs <- c(seen_pairs, pair_key)

    gene1 <- tryCatch(find_gene(rec$chrom, rec$pos),  error = function(e) NA_character_)
    gene2 <- tryCatch(find_gene(mate$chrom, mate$pos), error = function(e) NA_character_)

    # Skip rows without Hugo gene symbols - both sites must have gene annotation
    if (is.na(gene1) || is.na(gene2) || gene1 == "" || gene2 == "") next

    chr1 <- sub("^chr", "", rec$chrom)
    chr2 <- sub("^chr", "", mate$chrom)
    sv_class <- if (chr1 == chr2) "INVERSION" else "TRANSLOCATION"

    event_info <- paste0("DNA WGS Fusion: ", gene1, "--", gene2)

    sv_rows <- c(sv_rows, list(data.table(
        Sample_Id                  = opt$sample,
        SV_Status                  = opt$sv_status,
        Site1_Hugo_Symbol          = gene1,
        Site1_Ensembl_Transcript_Id= NA_character_,
        Site1_Region_Number        = NA_character_,
        Site1_Region               = NA_character_,
        Site2_Hugo_Symbol          = gene2,
        Site2_Ensembl_Transcript_Id= NA_character_,
        Site2_Region_Number        = NA_character_,
        Site2_Region               = NA_character_,
        Site2_Effect_On_Frame      = NA_character_,
        NCBI_Build                 = "GRCh38",
        Class                      = sv_class,
        DNA_Support                = "Yes",
        RNA_Support                = "No",
        Tumor_Variant_Count        = rec$split_reads,
        Connection_Type            = "5to3",
        Breakpoint_Type            = "PRECISE",
        Event_Info                 = event_info,
        Annotation                 = NA_character_,
        Site1_Chromosome           = chr1,
        Site1_Position             = rec$pos,
        Site2_Chromosome           = chr2,
        Site2_Position             = mate$pos,
        Tumor_Split_Read_Count     = rec$split_reads,
        Tumor_Paired_End_Read_Count= rec$paired_reads
    )))
}

if (length(sv_rows) > 0) {
    result <- rbindlist(sv_rows, fill = TRUE)
} else {
    cat("No PASS BND pairs found. Writing empty output.\n")
    result <- data.table(
        Sample_Id = character(), SV_Status = character(),
        Site1_Hugo_Symbol = character(), Site1_Ensembl_Transcript_Id = character(),
        Site1_Region_Number = character(), Site1_Region = character(),
        Site2_Hugo_Symbol = character(), Site2_Ensembl_Transcript_Id = character(),
        Site2_Region_Number = character(), Site2_Region = character(),
        Site2_Effect_On_Frame = character(), NCBI_Build = character(),
        Class = character(), DNA_Support = character(), RNA_Support = character(),
        Tumor_Variant_Count = integer(), Connection_Type = character(),
        Breakpoint_Type = character(), Event_Info = character(), Annotation = character(),
        Site1_Chromosome = character(), Site1_Position = integer(),
        Site2_Chromosome = character(), Site2_Position = integer(),
        Tumor_Split_Read_Count = integer(), Tumor_Paired_End_Read_Count = integer()
    )
}

cat("Writing SV output:", opt$output, "\n")
cat("Total SV records:", nrow(result), "\n")
write.table(result, opt$output, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
cat("Done!\n")
