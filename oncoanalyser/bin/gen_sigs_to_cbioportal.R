#!/usr/bin/env Rscript
# Convert HMFTOOLS SIGS allocation output to per-sample cBioPortal GENERIC_ASSAY format.
# Input: {subject}-T.sig.allocation.tsv (columns: signature, allocation, percent)
# Output: {sample}.data_sigs.txt (ENTITY_STABLE_ID, NAME, DESCRIPTION, URL, {sample})
# ENTITY_STABLE_ID / NAME / DESCRIPTION / URL are mapped from COSMIC v3.3 SBS signatures.

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
})

option_list <- list(
    make_option(c("-i", "--input"),  type = "character", help = "Path to .sig.allocation.tsv file [required]"),
    make_option(c("-s", "--sample"), type = "character", help = "Sample ID for output column header [required]"),
    make_option(c("-o", "--output"), type = "character", help = "Output file path [required]")
)

opt <- parse_args(OptionParser(option_list = option_list))

for (arg in c("input", "sample", "output")) {
    if (is.null(opt[[arg]])) stop(paste0("Missing required argument: --", arg))
}

# COSMIC v3.3 SBS signature metadata -----------------------------------------
cosmic_sbs <- data.table(
    id = c(
        "SBS1", "SBS2", "SBS3", "SBS4", "SBS5",
        "SBS6", "SBS7a", "SBS7b", "SBS7c", "SBS7d",
        "SBS8", "SBS9", "SBS10a", "SBS10b", "SBS10c", "SBS10d",
        "SBS11", "SBS12", "SBS13", "SBS14",
        "SBS15", "SBS16", "SBS17a", "SBS17b",
        "SBS18", "SBS19", "SBS20", "SBS21",
        "SBS22", "SBS22a", "SBS22b",
        "SBS23", "SBS24", "SBS25", "SBS26",
        "SBS27", "SBS28", "SBS29", "SBS30",
        "SBS31", "SBS32", "SBS33", "SBS34",
        "SBS35", "SBS36", "SBS37", "SBS38",
        "SBS39", "SBS40", "SBS40a", "SBS40b", "SBS40c",
        "SBS41", "SBS42", "SBS43", "SBS44",
        "SBS45", "SBS46", "SBS47", "SBS48",
        "SBS49", "SBS50", "SBS51", "SBS52",
        "SBS53", "SBS54", "SBS55", "SBS56",
        "SBS57", "SBS58",
        "SBS84", "SBS85", "SBS86", "SBS87",
        "SBS88", "SBS89", "SBS90", "SBS91",
        "SBS92", "SBS93", "SBS94"
    ),
    name = c(
        "SBS1", "SBS2", "SBS3", "SBS4", "SBS5",
        "SBS6", "SBS7a", "SBS7b", "SBS7c", "SBS7d",
        "SBS8", "SBS9", "SBS10a", "SBS10b", "SBS10c", "SBS10d",
        "SBS11", "SBS12", "SBS13", "SBS14",
        "SBS15", "SBS16", "SBS17a", "SBS17b",
        "SBS18", "SBS19", "SBS20", "SBS21",
        "SBS22", "SBS22a", "SBS22b",
        "SBS23", "SBS24", "SBS25", "SBS26",
        "SBS27", "SBS28", "SBS29", "SBS30",
        "SBS31", "SBS32", "SBS33", "SBS34",
        "SBS35", "SBS36", "SBS37", "SBS38",
        "SBS39", "SBS40", "SBS40a", "SBS40b", "SBS40c",
        "SBS41", "SBS42", "SBS43", "SBS44",
        "SBS45", "SBS46", "SBS47", "SBS48",
        "SBS49", "SBS50", "SBS51", "SBS52",
        "SBS53", "SBS54", "SBS55", "SBS56",
        "SBS57", "SBS58",
        "SBS84", "SBS85", "SBS86", "SBS87",
        "SBS88", "SBS89", "SBS90", "SBS91",
        "SBS92", "SBS93", "SBS94"
    ),
    description = c(
        "Spontaneous deamination of 5-methylcytosine; clock-like signature associated with age",
        "Activity of APOBEC family of cytidine deaminases",
        "Defective homologous recombination-based DNA repair (e.g., BRCA1/BRCA2 loss)",
        "Tobacco smoking (benzo[a]pyrene and other carcinogens)",
        "Unknown aetiology; clock-like signature associated with replication",
        "Defective DNA mismatch repair (MMR)",
        "Ultraviolet light exposure (pyrimidine dimers)",
        "Ultraviolet light exposure (pyrimidine dimers)",
        "Ultraviolet light exposure (pyrimidine dimers)",
        "Ultraviolet light exposure (pyrimidine dimers)",
        "Unknown aetiology; possibly late replication errors",
        "Polymerase eta activity in somatic hypermutation",
        "Defective POLE proofreading (polymerase epsilon exonuclease domain mutations)",
        "Defective POLE proofreading (polymerase epsilon exonuclease domain mutations)",
        "Defective POLD1 proofreading (polymerase delta exonuclease domain mutations)",
        "Defective POLD1 proofreading (polymerase delta exonuclease domain mutations)",
        "Alkylating agent chemotherapy (e.g., temozolomide)",
        "Unknown aetiology",
        "Activity of APOBEC family of cytidine deaminases",
        "Concurrent POLE exonuclease domain mutations and mismatch repair deficiency",
        "Defective DNA mismatch repair",
        "Unknown aetiology; associated with alcohol consumption in liver",
        "Unknown aetiology",
        "Unknown aetiology",
        "Damage from reactive oxygen species (ROS)",
        "Unknown aetiology",
        "Concurrent POLD1 exonuclease domain mutations and mismatch repair deficiency",
        "Defective DNA mismatch repair",
        "Aristolochic acid exposure",
        "Aristolochic acid exposure",
        "Aristolochic acid exposure",
        "Unknown aetiology",
        "Aflatoxin exposure",
        "Unknown aetiology; associated with prior chemotherapy",
        "Defective DNA mismatch repair",
        "Possible sequencing artifact",
        "Unknown aetiology",
        "Tobacco chewing",
        "Defective DNA base excision repair due to NTHL1 mutations",
        "Prior treatment with platinum drugs",
        "Prior treatment with azathioprine",
        "Unknown aetiology",
        "Unknown aetiology",
        "Prior treatment with platinum drugs",
        "Defective DNA base excision repair due to MUTYH mutations",
        "Unknown aetiology",
        "Indirect DNA damage from UV light",
        "Unknown aetiology",
        "Unknown aetiology; clock-like signature",
        "Unknown aetiology; clock-like signature",
        "Unknown aetiology; clock-like signature",
        "Unknown aetiology; clock-like signature",
        "Unknown aetiology",
        "Occupational exposure to haloalkanes",
        "Possible sequencing artifact",
        "Defective DNA mismatch repair",
        "Possible sequencing artifact",
        "Possible sequencing artifact",
        "Possible sequencing artifact",
        "Possible sequencing artifact",
        "Possible sequencing artifact",
        "Possible sequencing artifact",
        "Possible sequencing artifact",
        "Possible sequencing artifact",
        "Possible sequencing artifact",
        "Possible sequencing artifact",
        "Possible sequencing artifact",
        "Possible sequencing artifact",
        "Possible sequencing artifact",
        "Activity of activation-induced cytidine deaminase (AID)",
        "Activity of AID (somatic hypermutation at Ig loci)",
        "Unknown chemotherapy treatment",
        "Thiopurine chemotherapy treatment",
        "Colibactin exposure (E. coli genotoxin)",
        "Duocarmycin treatment",
        "Duocarmycin treatment",
        "Unknown aetiology",
        "Tobacco smoking",
        "Unknown aetiology",
        "Unknown aetiology"
    ),
    url = paste0("https://cancer.sanger.ac.uk/signatures/sbs/", tolower(c(
        "sbs1", "sbs2", "sbs3", "sbs4", "sbs5",
        "sbs6", "sbs7a", "sbs7b", "sbs7c", "sbs7d",
        "sbs8", "sbs9", "sbs10a", "sbs10b", "sbs10c", "sbs10d",
        "sbs11", "sbs12", "sbs13", "sbs14",
        "sbs15", "sbs16", "sbs17a", "sbs17b",
        "sbs18", "sbs19", "sbs20", "sbs21",
        "sbs22", "sbs22a", "sbs22b",
        "sbs23", "sbs24", "sbs25", "sbs26",
        "sbs27", "sbs28", "sbs29", "sbs30",
        "sbs31", "sbs32", "sbs33", "sbs34",
        "sbs35", "sbs36", "sbs37", "sbs38",
        "sbs39", "sbs40", "sbs40a", "sbs40b", "sbs40c",
        "sbs41", "sbs42", "sbs43", "sbs44",
        "sbs45", "sbs46", "sbs47", "sbs48",
        "sbs49", "sbs50", "sbs51", "sbs52",
        "sbs53", "sbs54", "sbs55", "sbs56",
        "sbs57", "sbs58",
        "sbs84", "sbs85", "sbs86", "sbs87",
        "sbs88", "sbs89", "sbs90", "sbs91",
        "sbs92", "sbs93", "sbs94"
    )))
)
setkey(cosmic_sbs, id)

# ---------------------------------------------------------------------------

cat("Reading SIGS allocation file:", opt$input, "\n")
sigs <- fread(opt$input, sep = "\t", header = TRUE)

expected_cols <- c("signature", "allocation", "percent")
missing_cols <- setdiff(expected_cols, colnames(sigs))
if (length(missing_cols) > 0) {
    stop("Missing columns in input: ", paste(missing_cols, collapse = ", "))
}

cat("Found", nrow(sigs), "allocated signatures\n")

# Map to COSMIC metadata; fall back to the raw ID for unknown signatures
matched <- cosmic_sbs[sigs$signature]
fallback_name <- ifelse(is.na(matched$name), sigs$signature, matched$name)
fallback_desc <- ifelse(
    is.na(matched$description),
    paste0("Mutational signature ", sigs$signature),
    matched$description
)
fallback_url <- ifelse(is.na(matched$url), NA_character_, matched$url)

unmapped <- sigs$signature[is.na(matched$name)]
if (length(unmapped) > 0) {
    warning("No COSMIC mapping found for: ", paste(unmapped, collapse = ", "),
            " — using raw signature ID as NAME/DESCRIPTION")
}

out <- data.table(
    ENTITY_STABLE_ID = sigs$signature,
    NAME             = fallback_name,
    DESCRIPTION      = fallback_desc,
    URL              = fallback_url
)
out[[opt$sample]] <- sigs$percent

write.table(out, opt$output, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
cat("Written", nrow(out), "rows to", opt$output, "\n")
