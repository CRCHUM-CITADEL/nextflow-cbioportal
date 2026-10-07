#!/usr/bin/env Rscript

# Convert Isofox RNA fusion output to cBioPortal data_sv.txt format.
#
# Expects the oncoanalyser 3.0 Isofox pass-fusion table
# (<sample>.isf.pass_fusions.tsv), whose columns are:
#   Name KnownType ChromosomeUp ChromosomeDown PositionUp PositionDown
#   OrientationUp OrientationDown JunctionTypeUp JunctionTypeDown
#   TranscriptUp TranscriptDown ExonUp ExonDown SvType
#   SplitFrags RealignedFrags DiscordantFrags DepthUp DepthDown
#   MaxAnchorLengthUp MaxAnchorLengthDown CohortFrequency
#
# Both gene symbols are packed into Name as "<up>_<down>"; either half may be
# empty when Isofox could not annotate that breakend. Such rows are dropped —
# data_sv.txt requires a Hugo symbol at both sites.
#
# KnownType, JunctionType*, SvType, Depth*, MaxAnchorLength* and CohortFrequency
# have no cBioPortal SV counterpart and are deliberately not emitted.

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
})

option_list <- list(
    make_option(c("-i", "--input"),  type = "character", help = "Isofox fusion TSV file (<sample>.isf.pass_fusions.tsv)"),
    make_option(c("-s", "--sample"), type = "character", help = "Tumor RNA sample ID"),
    make_option(c("-o", "--output"), type = "character", help = "Output data_sv.txt path")
)

opt <- parse_args(OptionParser(option_list = option_list))

for (arg in c("input", "sample", "output")) {
    if (is.null(opt[[arg]])) stop(paste("Missing required argument: --", arg, sep = ""))
}

cat("Reading Isofox fusion file:", opt$input, "\n")
fusions <- fread(opt$input, header = TRUE)
cat("Found", nrow(fusions), "fusion records\n")

if (nrow(fusions) == 0) {
    cat("No fusions found. Writing empty output.\n")

    empty <- data.table(
        Sample_Id = character(), SV_Status = character(),
        Site1_Hugo_Symbol = character(), Site1_Ensembl_Transcript_Id = character(),
        Site1_Region_Number = integer(), Site1_Region = character(),
        Site2_Hugo_Symbol = character(), Site2_Ensembl_Transcript_Id = character(),
        Site2_Region_Number = integer(), Site2_Region = character(),
        Site2_Effect_On_Frame = character(), NCBI_Build = character(),
        Class = character(), DNA_Support = character(), RNA_Support = character(),
        Tumor_Variant_Count = integer(), Connection_Type = character(),
        Breakpoint_Type = character(), Event_Info = character(), Annotation = character(),
        Site1_Chromosome = character(), Site1_Position = integer(),
        Site2_Chromosome = character(), Site2_Position = integer(),
        Tumor_Split_Read_Count = integer(), Tumor_Paired_End_Read_Count = integer()
    )
    write.table(empty, opt$output, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    cat("Done!\n")
    quit(status = 0)
}

required_cols <- c("Name", "ChromosomeUp", "ChromosomeDown", "PositionUp", "PositionDown",
                   "TranscriptUp", "TranscriptDown", "ExonUp", "ExonDown",
                   "SplitFrags", "RealignedFrags", "DiscordantFrags")
missing_cols <- setdiff(required_cols, names(fusions))
if (length(missing_cols) > 0) {
    stop(paste("Missing expected Isofox 3.0 columns:", paste(missing_cols, collapse = ", "),
               "\nAvailable:", paste(names(fusions), collapse = ", ")))
}

# Name is "<up gene>_<down gene>" — split on the first underscore.
# A Name without any underscore is malformed: keep it as the 5' gene only so the
# row is dropped by the Hugo-symbol filter rather than becoming a self-fusion.
fusion_name <- as.character(fusions$Name)
gene_up     <- sub("_.*$", "", fusion_name)
gene_down   <- ifelse(grepl("_", fusion_name, fixed = TRUE),
                      sub("^[^_]*_", "", fusion_name), "")

# Blank strings mean "not annotated" — carry them through as NA
blank_to_na <- function(x) {
    x <- as.character(x)
    x[is.na(x) | x == ""] <- NA_character_
    x
}

# Exon 0 is Isofox's "no exon" placeholder
exon_or_na <- function(x) {
    x <- suppressWarnings(as.integer(x))
    x[!is.na(x) & x == 0] <- NA_integer_
    x
}

exon_up   <- exon_or_na(fusions$ExonUp)
exon_down <- exon_or_na(fusions$ExonDown)

# cBioPortal's Site1/2_Region vocabulary is {5_Prime_UTR, 3_Prime_UTR, Promoter, Exon, Intron}.
#
# Isofox fills Transcript/Exon from the breakend's TransExonRef list
# (FusionReadData.java) — a transcript and exon rank are recorded only when the
# breakend falls on an exon of that transcript, and are left empty/0 otherwise.
# The fusion caller assigns a gene by overlap with the gene body and has no
# upstream/promoter allowance, so a named gene with no exon reference means the
# breakend sits inside that gene but outside every exon: intronic.
#
#   transcript + exon rank  -> Exon
#   gene but no exon ref    -> Intron
#   no gene                 -> NA (the row is dropped by the Hugo-symbol filter below)
#
# 5_Prime_UTR / 3_Prime_UTR / Promoter are not derivable from pass_fusions.tsv:
# it carries no coding-type or CDS boundaries (RnaFusionFile.java lists all 23
# columns), and UTR bases are exonic, so they surface here as Exon.
classify_region <- function(gene, transcript, exon) {
    ifelse(!is.na(exon) & !is.na(transcript), "Exon",
           ifelse(!is.na(gene) & gene != "", "Intron", NA_character_))
}

transcript_up   <- blank_to_na(fusions$TranscriptUp)
transcript_down <- blank_to_na(fusions$TranscriptDown)

int_or_zero <- function(x) {
    x <- suppressWarnings(as.integer(x))
    x[is.na(x)] <- 0L
    x
}

result <- data.table(
    Sample_Id                   = opt$sample,
    SV_Status                   = "SOMATIC",
    Site1_Hugo_Symbol           = gene_up,
    Site1_Ensembl_Transcript_Id = transcript_up,
    Site1_Region_Number         = exon_up,
    Site1_Region                = classify_region(gene_up, transcript_up, exon_up),
    Site2_Hugo_Symbol           = gene_down,
    Site2_Ensembl_Transcript_Id = transcript_down,
    Site2_Region_Number         = exon_down,
    Site2_Region                = classify_region(gene_down, transcript_down, exon_down),
    Site2_Effect_On_Frame       = NA_character_,
    NCBI_Build                  = "GRCh38",
    Class                       = "FUSION",
    DNA_Support                 = "No",
    RNA_Support                 = "Yes",
    # Isofox 3.0 dropped TotalFragments — sum the three supporting-fragment counts
    Tumor_Variant_Count         = int_or_zero(fusions$SplitFrags) +
                                  int_or_zero(fusions$RealignedFrags) +
                                  int_or_zero(fusions$DiscordantFrags),
    Connection_Type             = "5to3",
    Breakpoint_Type             = "PRECISE",
    Event_Info                  = paste0("RNA-seq Fusion: ", gene_up, "--", gene_down),
    Annotation                  = paste0(gene_up, " - ", gene_down, " fusion"),
    Site1_Chromosome            = sub("^chr", "", as.character(fusions$ChromosomeUp)),
    Site1_Position              = fusions$PositionUp,
    Site2_Chromosome            = sub("^chr", "", as.character(fusions$ChromosomeDown)),
    Site2_Position              = fusions$PositionDown,
    Tumor_Split_Read_Count      = suppressWarnings(as.integer(fusions$SplitFrags)),
    Tumor_Paired_End_Read_Count = suppressWarnings(as.integer(fusions$DiscordantFrags))
)

# Deduplicate on genomic coordinates
result <- unique(result, by = c("Site1_Chromosome", "Site1_Position",
                                "Site2_Chromosome", "Site2_Position"))

# Filter out rows with empty or NA Hugo symbols at either site
result <- result[
    !is.na(Site1_Hugo_Symbol) & Site1_Hugo_Symbol != "" &
    !is.na(Site2_Hugo_Symbol) & Site2_Hugo_Symbol != ""
]

cat("Writing fusion SV output:", opt$output, "\n")
cat("Total records:", nrow(result), "\n")
write.table(result, opt$output, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
cat("Done!\n")
