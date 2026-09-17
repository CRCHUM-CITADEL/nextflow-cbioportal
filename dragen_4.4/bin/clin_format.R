#!/usr/bin/env Rscript

library(optparse)

option_list <- list(
  make_option(c("-p", "--patient"), type="character", default=NULL,
              help="patient CSV [REQUIRED]", metavar="FILE"),
  make_option(c("-d", "--diagnosis"), type="character", default=NULL,
              help="diagnosis CSV [REQUIRED]", metavar="FILE"),
  make_option(c("-t", "--treatment"), type="character", default=NULL,
              help="treatment CSV [OPTIONAL]", metavar="FILE"),
  make_option(c("-s", "--surgeries"), type="character", default=NULL,
              help="surgeries CSV [OPTIONAL]", metavar="FILE"),
  make_option(c("-e", "--systreat"), type="character", default=NULL,
              help="systemic treatment CSV [OPTIONAL]", metavar="FILE"),
  make_option(c("-i", "--specimen"), type="character", default=NULL,
              help="specimen CSV [REQUIRED]", metavar="FILE"),
  make_option(c("-r", "--radiotherapy"), type="character", default=NULL,
              help="radio therapy CSV [OPTIONAL]", metavar="FILE"),
  make_option(c("-l", "--linking_file"), type="character", default=NULL,
              help="file to link tumor subjects to respective sample id [OPTIONAL]", metavar="FILE"),
  make_option(c("-o", "--output"), type="character", default="data_clinical_sample.txt",
              help="Output file path [default= %default]", metavar="FILE"),
  make_option(c("-m", "--mode"), type="character", default="sample",
              help="Between 'sample' or 'patient' mode.")
)

opt_parser <- OptionParser(option_list=option_list,
                          usage="Usage: %prog -p patient.csv -d diagnosis.csv -i specimen.csv [options]",
                          description="Merge clinical information into final sample file")
opt <- parse_args(opt_parser)

# Validation checks for required arguments
if (is.null(opt$patient)) {
  stop("Error: --patient argument is required")
}
if (is.null(opt$diagnosis)) {
  stop("Error: --diagnosis argument is required")
}
if (is.null(opt$specimen)) {
  stop("Error: --specimen argument is required")
}

# Check that provided files exist
input_files <- c(opt$patient, opt$diagnosis, opt$treatment, opt$surgeries,
                 opt$systreat, opt$specimen, opt$radiotherapy)

for (file in input_files[!sapply(input_files, is.null)]) {
  if (!file.exists(file)) {
    stop(paste("Error: File does not exist:", file))
  }
}
# Print parameters
cat("=== Parameters ===\n")
cat(sprintf("Patient:        %s\n", ifelse(is.null(opt$patient), "NULL", opt$patient)))
cat(sprintf("Diagnosis:      %s\n", ifelse(is.null(opt$diagnosis), "NULL", opt$diagnosis)))
cat(sprintf("Treatment:      %s\n", ifelse(is.null(opt$treatment), "NULL", opt$treatment)))
cat(sprintf("Surgeries:      %s\n", ifelse(is.null(opt$surgeries), "NULL", opt$surgeries)))
cat(sprintf("Systemic Treat: %s\n", ifelse(is.null(opt$systreat), "NULL", opt$systreat)))
cat(sprintf("Specimen:       %s\n", ifelse(is.null(opt$specimen), "NULL", opt$specimen)))
cat(sprintf("Radiotherapy:   %s\n", ifelse(is.null(opt$radiotherapy), "NULL", opt$radiotherapy)))
cat(sprintf("Linking file:   %s\n", ifelse(is.null(opt$linking_file), "NULL", opt$linking_file)))
cat(sprintf("Output:         %s\n", opt$output))
cat(sprintf("Mode:           %s\n", opt$mode))
cat("==================\n\n")

# Read in files
cat("Reading in files: \n")
sampleList=read.csv(opt$linking_file, header=T, sep='\t')
patient_data=read.csv(opt$patient, header=T, sep=",")
diagnostic_data=read.csv(opt$diagnosis, header=T, sep=",")
specimen_data=read.csv(opt$specimen, header=T, sep=",")

colnames(sampleList)=c("patient", "sample")

# Patient data formatting

colnames(patient_data)[1]="patient"
patient_data=patient_data[,c("patient","sex_at_birth","date_of_birth","is_deceased")]
patient_data$sex_at_birth=gsub("Female","F",patient_data$sex_at_birth)
patient_data$sex_at_birth=gsub("Male","M",patient_data$sex_at_birth)
patient_data$is_deceased=gsub("YES","1:DECEASED",patient_data$is_deceased)
patient_data$is_deceased=gsub("NO","0:LIVING",patient_data$is_deceased)

# Diagnosis data formatting
colnames(diagnostic_data)[1]="patient"
diagnostic_data=diagnostic_data[,c("patient","date_of_diagnosis","cancer_type_code","clinical_stage_group","submitter_primary_diagnosis_id", "date_of_data_extraction")]
diagnostic_data=diagnostic_data[order(diagnostic_data$patient, diagnostic_data$date_of_diagnosis),]
diagnostic_data=diagnostic_data[!duplicated(diagnostic_data$patient), ]

# Specimen dat formatting
colnames(specimen_data)[1]="patient"
specimen_data=specimen_data[,c("patient","tumour_histological_type","submitter_primary_diagnosis_id")]

#merge
m=merge(sampleList, patient_data, by='patient', all.x=TRUE)

print(head(m))

m=merge(m,diagnostic_data, by='patient',all.x=TRUE)
m=merge(m,specimen_data, by='patient',all.x=TRUE)
m=m[(m$submitter_primary_diagnosis_id.x==m$submitter_primary_diagnosis_id.y) | is.na(m$submitter_primary_diagnosis_id.x), ]
m=subset(m, select = -c(submitter_primary_diagnosis_id.x,submitter_primary_diagnosis_id.y))


# Age at diagnosis, overall survival in months since initial diagnosis
m$age_at_diagnosis_years=as.numeric((as.Date(m$date_of_diagnosis)-as.Date(m$date_of_birth))/365)
m$age_at_diagnosis_years=round(m$age_at_diagnosis_years,3)
m$overall_survival_months[m$is_deceased == "1:DECEASED" | is.na(m$is_deceased)]=NA

living_idx <- m$is_deceased == "0:LIVING" & !is.na(m$is_deceased)
m$overall_survival_months[living_idx] <- as.numeric((as.Date(m$date_of_data_extraction[living_idx]) - as.Date(m$date_of_diagnosis[living_idx])) / 30)
m$overall_survival_months <- round(m$overall_survival_months, 3)

# Fill in missing information for now
m$disease_free_survival_status=NA
m$disease_free_months=NA
m$sample_type="Primary"
m$cancer_type_details="NA"
m$tumour_site="NA"

# Fill empty values with NA
m[m==""]=NA

print(paste("writing with mode", opt$mode))
print(str(m))

# Keep final columns
if(opt$mode=="patient") {
	writeLines(c("#Patient Identifier	Sex	Diagnosis Age	Tumor Site	Tumor Grade	Tumor Histological Type	Overall Survival (Months)	Overall Survival Status	Disease Free Status	Disease free Months",
				"#Identifier to uniquely specify a patient.	Sex.	Tumor site.	Tumor Grade.	Tumor histological type.	Age at which a condition or disease was first diagnosed.	Overall survival in months since initial diagonosis.	Overall patient survival status.	 Disease free status since initial treatment.	Disease free (months) since initial treatment.",
				"#STRING	STRING	NUMBER	STRING	STRING	STRING	NUMBER	STRING	STRING	NUMBER",
				"#1	1	1	1	1	1	1	1	1	1",
				"PATIENT_ID	SEX	AGE	TUMOR_SITE	TUMOR_GRADE	TUMOR_HISTOLOGICAL_TYPE	OS_MONTHS	OS_STATUS	DFS_STATUS	DFS_MONTHS"), con = opt$output)
	write.table(m[,c("patient","sex_at_birth","age_at_diagnosis_years","tumour_site","clinical_stage_group","tumour_histological_type","overall_survival_months","is_deceased","disease_free_survival_status","disease_free_months")], file = opt$output, sep = "\t", row.names = FALSE, col.names=F, quote = FALSE, append = TRUE)
} else if(opt$mode=="sample") {
	writeLines(c("#Patient Identifier	Sample Identifier	Sample Type	Cancer Type	Cancer Code",
        "#Identifier to uniquely specify a patient.	A unique sample identifier.	The type of sample (i.e., normal, primary, met, recurrence).	Cancer Type Details.	Cancer Type Code.",
				"#STRING	STRING	STRING	STRING	STRING",
				"#1	1	1	1	1",
	   		"PATIENT_ID	SAMPLE_ID	SAMPLE_TYPE	CANCER_TYPE_DETAILS	CANCER_TYPE_CODE"), con = opt$output)
	# here we asssume patient_id and sample_id are the same. might not be the case
	write.table(m[,c("patient", "sample", "sample_type","cancer_type_details","cancer_type_code")], file = opt$output, sep = "\t", row.names = FALSE, col.names=F, quote = FALSE, append = TRUE)
}
