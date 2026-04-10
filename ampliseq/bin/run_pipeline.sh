BASE_DIR="/data/optilab/patients/Cohorte_10_patients_new"
SCRIPTS_DIR="/shared/cbioportal/formatting"
DATA_DIR="/data/optilab/cbioportal"

FORMAT_SCRIPT="${SCRIPTS_DIR}/format_tsv.py"
FORMAT_CNA="${SCRIPTS_DIR}/format_cna.py"
CLINICAL_PATIENT_SCRIPT="${SCRIPTS_DIR}/clinical_patients_format.py"
CLINICAL_SAMPLE_SCRIPT="${SCRIPTS_DIR}/clinical_sample_format.py"
FORMAT_MUTATIONS_SCRIPT="${SCRIPTS_DIR}/format_mutations.py"
FORMAT_SV_SCRIPT="${SCRIPTS_DIR}/format_sv.py"
FORMAT_CNA_DEANON_SCRIPT="${SCRIPTS_DIR}/format_cna_deanon.py"
FORMAT_META_SCRIPT="${SCRIPTS_DIR}/format_meta.py"


APPTAINER_SIF="/shared/cbioportal/formatting/vcf2maf_ensembl-vep_v1.6.22_2.sif"
LINKING_FILE="/shared/cbioportal/data_mount/data_copy/linking_file.txt"

OUT_DIR="$(pwd)/output"
WORK_DIR="$(pwd)/work"
mkdir -p "$OUT_DIR" "$WORK_DIR"

SAMPLE_IDS=()

echo "Test"

for SAMPLE_DIR in "$BASE_DIR"/*/; do
    SAMPLE=$(basename "$SAMPLE_DIR")
    echo "Processing: $SAMPLE"

    # Resolve VCF early to derive SAMPLE_ID
    VCF_GZ=$(ls "$SAMPLE_DIR"*-basespace-pisces.final.vcf.gz 2>/dev/null | head -1)
    [[ -z "$VCF_GZ" ]] && { echo "  SKIP: no VCF found"; continue; }
    SAMPLE_ID=$(basename "$VCF_GZ" | sed 's/-basespace.*//')
    SAMPLE_IDS+=("$SAMPLE_ID")

    # 1. Format TSV export — writes data_sv.txt to output directory
    TSV=$(ls "$SAMPLE_DIR"analysis_*_export.tsv 2>/dev/null | head -1)
    [[ -z "$TSV" ]] && { echo "  SKIP: no TSV found"; continue; }
    (cd "$OUT_DIR" && python3.12 "$FORMAT_SCRIPT" "$TSV" "$SAMPLE_ID")
    # 1b. Format CNA - appends to data_cna.txt in output directory
    (cd "$OUT_DIR" && python3.12 "$FORMAT_CNA" "$TSV" "$SAMPLE_ID")
    # 2. Decompress VCF to work directory and run Apptainer
    VCF="${WORK_DIR}/${SAMPLE_ID}.final.vcf"
    gzip -dc "$VCF_GZ" > "$VCF"

    APPTAINER_OUT="${WORK_DIR}/${SAMPLE_ID}-basespace-pisces.final.maf"
    apptainer exec \
        --bind ${DATA_DIR}:/home/jbellavance/ \
        "$APPTAINER_SIF" vcf2maf.pl \
        --tumor-id "$SAMPLE_ID" --cache-version 113 \
        --ncbi-build GRCh37 --vep-path /opt/conda/bin \
        --vep-data ~/ --ref-fasta ~/hg19.fa \
        --input-vcf "$VCF" --output-maf "$APPTAINER_OUT"

    # Strip first header line from MAF (keep only the second line as the column header)
    sed -i '1d' "$APPTAINER_OUT"

    # 3. AWK post-processing — append filtered rows to final output
    FINAL_OUT="${OUT_DIR}/data_mutations.txt"
    awk -v skip_header="$([ -s "$FINAL_OUT" ] && echo 1 || echo 0)" \
        'NR==FNR {regions["chr"$1":"$2"-"$3]=1; next} FNR==1 {if (!skip_header) print; next} $5":"$6"-"$7 in regions' \
        "$TSV" \
        "$APPTAINER_OUT" >> "$FINAL_OUT"

    echo "  Done: $FINAL_OUT (appended $SAMPLE_ID)"
done

# Deanonymise clinical and data files
echo "Deanonymising clinical and data files..."
(cd "$OUT_DIR" && python3.12 "$CLINICAL_PATIENT_SCRIPT" ${BASE_DIR}/Patient_file_cohort_test.txt)
(cd "$OUT_DIR" && python3.12 "$CLINICAL_SAMPLE_SCRIPT"  ${BASE_DIR}/Sample_file_cohort_test.txt)
( cd "$OUT_DIR" && python3.12 "$FORMAT_MUTATIONS_SCRIPT" data_mutations.txt "$LINKING_FILE")
( cd "$OUT_DIR" && python3.12 "$FORMAT_SV_SCRIPT"        data_sv.txt        "$LINKING_FILE")
( cd "$OUT_DIR" && python3.12 "$FORMAT_CNA_DEANON_SCRIPT" data_cna.txt 	 "$LINKING_FILE") 

# and in the deanonymisation block at the bottom:
# Write case_lists files using deanonymised sample IDs from linking file
mkdir -p "${OUT_DIR}/case_lists"
IDS_TAB=$(awk 'NR>1 {printf "%s\t", $2} END {print ""}' "$LINKING_FILE" | sed 's/\t$//')

cat > "${OUT_DIR}/case_lists/cases_sequenced.txt" << CEOF
cancer_study_identifier: optilab_test
stable_id: optilab_test_sequenced
case_list_name: all mutations
case_list_description: all mutations of optilab test
case_list_ids: ${IDS_TAB}
CEOF

cat > "${OUT_DIR}/case_lists/cases_sv.txt" << CEOF
cancer_study_identifier: optilab_test
stable_id: optilab_test_sv
case_list_name: all sv
case_list_description: all sv of optilab test
case_list_ids: ${IDS_TAB}
CEOF

cat > "${OUT_DIR}/case_lists/cases_cna.txt" << CEOF
cancer_study_identifier: optilab_test
stable_id: optilab_test_cna
case_list_name: all cna
case_list_description: all cna of optilab test
case_list_ids: ${IDS_TAB}
CEOF

echo "Written: ${OUT_DIR}/case_lists/cases_sequenced.txt, cases_cna.txt and cases_sv.txt (${#SAMPLE_IDS[@]} samples)"

# Generate cBioPortal meta files
python3.12 "$FORMAT_META_SCRIPT" "optilab_study" "$OUT_DIR"
echo "Written: meta files to ${OUT_DIR}"
