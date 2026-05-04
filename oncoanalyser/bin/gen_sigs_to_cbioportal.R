#!/usr/bin/env Rscript
# Convert HMFTOOLS SIGS allocation output to per-sample cBioPortal GENERIC_ASSAY format.
# Input: {subject}-T.sig.allocation.tsv (columns: signature, allocation, percent)
#        signature column uses Sig1/Sig2 notation.
# Output: {sample}.data_mutational_signatures_contribution_SBS.txt
#         (ENTITY_STABLE_ID, NAME, DESCRIPTION, {sample})
# DESCRIPTION is "{etiology} ({main effect})" from signatures_etiology.tsv.

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
})

option_list <- list(
    make_option(c("-i", "--input"),    type = "character", help = "Path to .sig.allocation.tsv file [required]"),
    make_option(c("-s", "--sample"),   type = "character", help = "Sample ID for output column header [required]"),
    make_option(c("-e", "--etiology"), type = "character", help = "Path to signatures_etiology.tsv [required]"),
    make_option(c("-o", "--output"),   type = "character", help = "Output file path [required]")
)

opt <- parse_args(OptionParser(option_list = option_list))

for (arg in c("input", "sample", "etiology", "output")) {
    if (is.null(opt[[arg]])) stop(paste0("Missing required argument: --", arg))
}

# Short broad category labels for the NAME column, keyed by Sig notation
categories <- c(
    Sig1  = "Aging",
    Sig2  = "APOBEC",
    Sig3  = "HR deficiency",
    Sig4  = "Tobacco",
    Sig5  = "Clock-like",
    Sig6  = "MMR deficiency",
    Sig7  = "UV",
    Sig8  = "Unknown",
    Sig9  = "AID activity",
    Sig10 = "POLE deficiency",
    Sig11 = "Chemotherapy",
    Sig12 = "Unknown",
    Sig13 = "APOBEC",
    Sig14 = "MMR deficiency",
    Sig15 = "MMR deficiency",
    Sig16 = "Unknown",
    Sig17 = "Unknown",
    Sig18 = "ROS",
    Sig19 = "Unknown",
    Sig20 = "MMR deficiency",
    Sig21 = "MMR deficiency",
    Sig22 = "Carcinogen",
    Sig23 = "Unknown",
    Sig24 = "Carcinogen",
    Sig25 = "Chemotherapy",
    Sig26 = "MMR deficiency",
    Sig27 = "Unknown",
    Sig28 = "Unknown",
    Sig29 = "Tobacco",
    Sig30 = "BER deficiency"
)

# Main biological effects for parenthetical annotation, keyed by Sig notation
main_effects <- c(
    Sig1  = "aging",
    Sig2  = "breast, bladder, cervical cancer",
    Sig3  = "BRCA1/2 deficiency, breast, ovarian cancer",
    Sig4  = "lung cancer",
    Sig5  = "ubiquitous, clock-like",
    Sig6  = "colorectal, endometrial cancer",
    Sig7  = "skin cancer (melanoma)",
    Sig8  = "breast, colorectal cancer",
    Sig9  = "B-cell lymphoma, CLL",
    Sig10 = "colorectal, endometrial cancer",
    Sig11 = "glioma",
    Sig12 = "liver cancer",
    Sig13 = "breast, bladder, cervical cancer",
    Sig14 = "colorectal, endometrial cancer",
    Sig15 = "stomach, colorectal cancer",
    Sig16 = "liver cancer",
    Sig17 = "esophageal, stomach cancer",
    Sig18 = "neuroblastoma, breast cancer",
    Sig19 = "pilocytic astrocytoma",
    Sig20 = "breast cancer",
    Sig21 = "colorectal cancer",
    Sig22 = "urothelial cancer",
    Sig23 = "many cancer types",
    Sig24 = "liver cancer",
    Sig25 = "lymphoma",
    Sig26 = "breast, colorectal cancer",
    Sig27 = "kidney cancer",
    Sig28 = "many cancer types",
    Sig29 = "oral cancer",
    Sig30 = "breast cancer"
)

cat("Reading SIGS allocation file:", opt$input, "\n")
sigs <- fread(opt$input, sep = "\t", header = TRUE)

expected_cols <- c("signature", "allocation", "percent")
missing_cols <- setdiff(expected_cols, colnames(sigs))
if (length(missing_cols) > 0) {
    stop("Missing columns in input: ", paste(missing_cols, collapse = ", "))
}

cat("Reading etiology file:", opt$etiology, "\n")
etiology <- fread(opt$etiology, sep = "\t", header = TRUE)
setkey(etiology, signature)

cat("Found", nrow(sigs), "allocated signatures\n")

# Input uses Sig notation (Sig1); etiology file and categories/main_effects are keyed the same way
sig_key <- sigs$signature

matched_etiology <- etiology[sig_key, etiology]
matched_effect   <- main_effects[sig_key]

unmapped <- sigs$signature[is.na(matched_etiology)]
if (length(unmapped) > 0) {
    warning("No etiology mapping found for: ", paste(unmapped, collapse = ", "),
            " — using signature ID as DESCRIPTION")
}

descriptions <- ifelse(
    is.na(matched_etiology),
    paste0("Mutational signature ", sigs$signature),
    ifelse(
        !is.na(matched_effect),
        paste0(matched_etiology, " (", matched_effect, ")"),
        matched_etiology
    )
)

sig_num      <- gsub("Sig", "", sigs$signature)
matched_cat  <- categories[sig_key]
sig_names    <- ifelse(
    !is.na(matched_cat),
    paste0("SBS", sig_num, " (", matched_cat, ")"),
    paste0("SBS", sig_num)
)

out <- data.table(
    ENTITY_STABLE_ID = paste0("mutational_signatures_contribution_SBS", sig_num),
    NAME             = sig_names,
    DESCRIPTION      = descriptions
)
out[[opt$sample]] <- sigs$percent

write.table(out, opt$output, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
cat("Written", nrow(out), "rows to", opt$output, "\n")
