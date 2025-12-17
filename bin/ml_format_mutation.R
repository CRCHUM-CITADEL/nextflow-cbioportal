#!/usr/bin/env Rscript

# Single MAF file processor:
# → Takes a single pre-combined MAF file
# → Splits it into somatic and germline mutations
# → Outputs two separate files

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Error: Provide a MAF file.\nUsage: Rscript format_mutation.R <maf_file>")
}

library(tidyverse)

input_file <- args[1]

if (!file.exists(input_file)) {
  stop(paste("Error: File", input_file, "does not exist."))
}

cat("Processing MAF file:", input_file, "\n\n")

# --- Amino acid conversion ---
convert_aa_code <- function(aa_change) {
  if (is.na(aa_change)) return(NA_character_)
  
  aa_map <- c(
    "Ala"="A","Arg"="R","Asn"="N","Asp"="D",
    "Cys"="C","Gln"="Q","Glu"="E","Gly"="G",
    "His"="H","Ile"="I","Leu"="L","Lys"="K",
    "Met"="M","Phe"="F","Pro"="P","Ser"="S",
    "Thr"="T","Trp"="W","Tyr"="Y","Val"="V",
    "Ter"="*","Xxx"="X"
  )
  
  result <- aa_change
  for (three in names(aa_map)) {
    result <- gsub(three, aa_map[[three]], result, ignore.case = TRUE)
  }
  return(result)
}

# --- Read and process the MAF file ---
maf <- read_tsv(input_file, comment = "#", show_col_types = FALSE)

# Transform
transformed <- maf %>%
  mutate(
    sample = as.character(Tumor_Sample_Barcode),
    gene = as.character(Hugo_Symbol),
    chrom = as.character(Chromosome),
    start = as.integer(Start_Position),
    end = as.integer(End_Position),
    ref = as.character(Reference_Allele),
    alt = as.character(Tumor_Seq_Allele2),
    Amino_Acid_Change = as.character(sapply(HGVSp_Short, convert_aa_code)),
    effect = as.character(Variant_Classification),
    callers = as.character(NA_character_),
    dna_vaf = as.double(ifelse(t_depth > 0, t_alt_count / t_depth, NA_real_)),
  ) %>%
  select(sample, gene, chrom, start, end, ref, alt,
         Tumor_Sample_Barcode, Amino_Acid_Change, effect,
         callers, dna_vaf, Mutation_Status)

# Split into somatic and germline
somatic <- transformed %>% 
  filter(Mutation_Status != "Germline") %>%
  select(-Mutation_Status)

germline <- transformed %>% 
  filter(Mutation_Status == "Germline") %>%
  select(-Mutation_Status)

cat("  → Somatic mutations:", nrow(somatic), "\n")
cat("  → Germline mutations:", nrow(germline), "\n\n")

# --- Write final files ---
write_tsv(somatic, file.path(getwd(), "all_somatic_mutations.tsv"))
write_tsv(germline, file.path(getwd(), "all_germline_mutations.tsv"))

cat("✔ DONE!\n")
cat("Final files created:\n")
cat("  -", file.path(getwd(), "all_somatic_mutations.tsv"), "\n")
cat("  -", file.path(getwd(), "all_germline_mutations.tsv"), "\n")
cat("Somatic rows:", nrow(somatic), "\n")
cat("Germline rows:", nrow(germline), "\n")