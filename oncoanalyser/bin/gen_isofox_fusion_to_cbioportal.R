#!/usr/bin/env Rscript

# Convert Isofox RNA fusion output to cBioPortal data_sv.txt format.
#
# Isofox fusion.tsv column naming can vary slightly by version.
# This script auto-detects the 5'/3' gene, chromosome, and position columns.
# Typical column sets:
#   - GeneName / OtherGeneName  (Isofox ≥v1.5)
#   - GeneNameUp / GeneNameDown (older)
#   - 5pGeneName / 3pGeneName   (some builds)

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
})

option_list <- list(
    make_option(c("-i", "--input"),  type = "character", help = "Isofox fusion TSV file (<sample>.isofox.fusion.tsv)"),
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

# Helper: pick first matching column name from a priority list
pick_col <- function(cols_available, candidates) {
    for (c in candidates) {
        if (c %in% cols_available) return(c)
    }
    stop(paste("None of the expected columns found:", paste(candidates, collapse = ", "),
               "\nAvailable:", paste(cols_available, collapse = ", ")))
}

if (nrow(fusions) == 0) {
    cat("No fusions found. Writing empty output.\n")

    empty <- data.table(
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
    write.table(empty, opt$output, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    cat("Done!\n")
    quit(status = 0)
}

cols <- names(fusions)

# Auto-detect 5' gene info columns
gene1_col <- pick_col(cols, c("GeneName", "GeneNameUp", "5pGeneName", "gene_name_up", "GeneA"))
chr1_col  <- pick_col(cols, c("Chromosome", "ChrUp", "5pChromosome", "chr_up", "ChrA"))
pos1_col  <- pick_col(cols, c("Position", "PosUp", "JunctionPositionUp", "pos_up", "PosA", "JunctionPosition"))

# Auto-detect 3' gene info columns
gene2_col <- pick_col(cols, c("OtherGeneName", "GeneNameDown", "3pGeneName", "gene_name_down", "GeneB"))
chr2_col  <- pick_col(cols, c("OtherChromosome", "ChrDown", "3pChromosome", "chr_down", "ChrB"))
pos2_col  <- pick_col(cols, c("OtherPosition", "PosDown", "JunctionPositionDown", "pos_down", "PosB"))

# Auto-detect read count columns
split_col   <- Filter(function(c) c %in% cols,
                      c("SplitFragments", "JunctionReadCount", "SpliceReadCount",
                        "SplitFragmentCount", "junction_read_count"))[1]
discord_col <- Filter(function(c) c %in% cols,
                      c("DiscordantPairs", "DiscordantFragments", "DiscordantFragmentCount",
                        "discordant_fragment_count", "DiscordantReads"))[1]
total_col   <- Filter(function(c) c %in% cols,
                      c("TotalFragments", "total_fragments"))[1]

cat("Using columns: 5'gene=", gene1_col, ", 5'chr=", chr1_col, ", 5'pos=", pos1_col, "\n")
cat("              3'gene=", gene2_col, ", 3'chr=", chr2_col, ", 3'pos=", pos2_col, "\n")

result <- data.table(
    Sample_Id                   = opt$sample,
    SV_Status                   = "SOMATIC",
    Site1_Hugo_Symbol           = fusions[[gene1_col]],
    Site1_Ensembl_Transcript_Id = NA_character_,
    Site1_Region_Number         = NA_character_,
    Site1_Region                = NA_character_,
    Site2_Hugo_Symbol           = fusions[[gene2_col]],
    Site2_Ensembl_Transcript_Id = NA_character_,
    Site2_Region_Number         = NA_character_,
    Site2_Region                = NA_character_,
    Site2_Effect_On_Frame       = NA_character_,
    NCBI_Build                  = "GRCh38",
    Class                       = "FUSION",
    DNA_Support                 = "No",
    RNA_Support                 = "Yes",
    Tumor_Variant_Count         = if (length(total_col) > 0 && !is.na(total_col))
                                        fusions[[total_col]]
                                    else if (length(split_col) > 0 && !is.na(split_col))
                                        fusions[[split_col]]
                                    else NA_integer_,
    Connection_Type             = "5to3",
    Breakpoint_Type             = "PRECISE",
    Event_Info                  = paste0("RNA-seq Fusion: ", fusions[[gene1_col]], "--", fusions[[gene2_col]]),
    Annotation                  = paste0(fusions[[gene1_col]], " - ", fusions[[gene2_col]], " fusion"),
    Site1_Chromosome            = sub("^chr", "", as.character(fusions[[chr1_col]])),
    Site1_Position              = fusions[[pos1_col]],
    Site2_Chromosome            = sub("^chr", "", as.character(fusions[[chr2_col]])),
    Site2_Position              = fusions[[pos2_col]],
    Tumor_Split_Read_Count      = if (length(split_col) > 0 && !is.na(split_col))
                                        fusions[[split_col]] else NA_integer_,
    Tumor_Paired_End_Read_Count = if (length(discord_col) > 0 && !is.na(discord_col))
                                        fusions[[discord_col]] else NA_integer_
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
