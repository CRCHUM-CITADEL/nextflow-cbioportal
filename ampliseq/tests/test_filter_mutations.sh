#!/usr/bin/env bash
# tests/test_filter_mutations.sh
#
# Tests FILTER_MUTATIONS and PASSTHROUGH_MUTATIONS logic using the awk/cp commands
# from modules/local/filter_mutations/main.nf and modules/local/passthrough_mutations/main.nf
#
# Two scenarios:
#   1. All MAF mutations match TSV coordinates → filtered == passthrough (same count)
#   2. MAF has extra mutations NOT in TSV     → passthrough > filtered (non-TSV variants preserved)
#
# Usage: bash tests/test_filter_mutations.sh
# Run from project root.

set -euo pipefail

TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

PASS=0
FAIL=0

check() {
    local desc="$1" expected="$2" actual="$3"
    if [ "$expected" = "$actual" ]; then
        echo "PASS: $desc"
        PASS=$((PASS + 1))
    else
        echo "FAIL: $desc (expected=$expected, actual=$actual)"
        FAIL=$((FAIL + 1))
    fi
}

# Mirrors FILTER_MUTATIONS awk from modules/local/filter_mutations/main.nf
filter_maf() {
    local maf="$1" tsv="$2" out="$3"
    awk '
        NR==FNR { regions["chr" $1 ":" $2 "-" $3] = 1; next }
        FNR==1  { print; next }
        $5 ":" $6 "-" $7 in regions
    ' "$tsv" "$maf" > "$out"
}

# Mirrors PASSTHROUGH_MUTATIONS cp from modules/local/passthrough_mutations/main.nf
passthrough_maf() {
    local maf="$1" out="$2"
    cp "$maf" "$out"
}

MAF_HEADER="Hugo_Symbol	Entrez_Gene_Id	Center	NCBI_Build	Chromosome	Start_Position	End_Position	Strand	Variant_Classification	Variant_Type	Reference_Allele	Tumor_Seq_Allele1	Tumor_Seq_Allele2	dbSNP_RS	dbSNP_Val_Status	Tumor_Sample_Barcode"

# ─────────────────────────────────────────────────────────────────────────────
# Scenario 1: all MAF mutations have exact coordinate matches in TSV
#   Expected: filter == passthrough (same line count)
# ─────────────────────────────────────────────────────────────────────────────
echo "--- Scenario 1: all mutations in TSV (filter == passthrough) ---"

TSV_001="assets/samples/SAMPLE_001/analysis_001_export.tsv"
TSV_002="assets/samples/SAMPLE_002/analysis_002_export.tsv"

printf '%s\n' "$MAF_HEADER" > "$TMPDIR/s001_all_match.maf"
printf 'BRAF\t673\ttest\thg19\tchr1\t150000\t150000\t+\tMissense_Mutation\tSNP\tA\tT\tT\t.\t.\tSAMPLE_001\n' >> "$TMPDIR/s001_all_match.maf"
printf 'TP53\t7157\ttest\thg19\tchr17\t7571720\t7571720\t+\tMissense_Mutation\tSNP\tG\tA\tA\t.\t.\tSAMPLE_001\n' >> "$TMPDIR/s001_all_match.maf"

filter_maf     "$TMPDIR/s001_all_match.maf" "$TSV_001" "$TMPDIR/s001_filtered.txt"
passthrough_maf "$TMPDIR/s001_all_match.maf"            "$TMPDIR/s001_passthrough.txt"

filtered_n=$(wc -l < "$TMPDIR/s001_filtered.txt")
passthrough_n=$(wc -l < "$TMPDIR/s001_passthrough.txt")

check "SAMPLE_001: filter keeps all TSV-matching mutations" "3" "$filtered_n"
check "SAMPLE_001: passthrough keeps all mutations" "3" "$passthrough_n"
check "SAMPLE_001: filter == passthrough count" "$passthrough_n" "$filtered_n"

printf '%s\n' "$MAF_HEADER" > "$TMPDIR/s002_all_match.maf"
printf 'PIK3CA\t5290\ttest\thg19\tchr10\t89685251\t89685251\t+\tMissense_Mutation\tSNP\tT\tG\tG\t.\t.\tSAMPLE_002\n' >> "$TMPDIR/s002_all_match.maf"
printf 'KRAS\t3845\ttest\thg19\tchr12\t25380275\t25380275\t+\tMissense_Mutation\tSNP\tC\tA\tA\t.\t.\tSAMPLE_002\n' >> "$TMPDIR/s002_all_match.maf"

filter_maf     "$TMPDIR/s002_all_match.maf" "$TSV_002" "$TMPDIR/s002_filtered.txt"
passthrough_maf "$TMPDIR/s002_all_match.maf"            "$TMPDIR/s002_passthrough.txt"

filtered_n=$(wc -l < "$TMPDIR/s002_filtered.txt")
passthrough_n=$(wc -l < "$TMPDIR/s002_passthrough.txt")

check "SAMPLE_002: filter keeps all TSV-matching mutations" "3" "$filtered_n"
check "SAMPLE_002: passthrough keeps all mutations" "3" "$passthrough_n"
check "SAMPLE_002: filter == passthrough count" "$passthrough_n" "$filtered_n"

# ─────────────────────────────────────────────────────────────────────────────
# Scenario 2: MAF has mutations NOT in TSV (the real-world bug case)
#   Passthrough must preserve them; filter must drop them.
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "--- Scenario 2: mutations NOT in TSV (passthrough > filter) ---"

# chr5:1234567 and chr9:9876543 are not in any TSV
printf '%s\n' "$MAF_HEADER" > "$TMPDIR/s001_extra.maf"
printf 'BRAF\t673\ttest\thg19\tchr1\t150000\t150000\t+\tMissense_Mutation\tSNP\tA\tT\tT\t.\t.\tSAMPLE_001\n' >> "$TMPDIR/s001_extra.maf"
printf 'TP53\t7157\ttest\thg19\tchr17\t7571720\t7571720\t+\tMissense_Mutation\tSNP\tG\tA\tA\t.\t.\tSAMPLE_001\n' >> "$TMPDIR/s001_extra.maf"
printf 'UNKNOWN\t0\ttest\thg19\tchr5\t1234567\t1234567\t+\tMissense_Mutation\tSNP\tC\tT\tT\t.\t.\tSAMPLE_001\n'  >> "$TMPDIR/s001_extra.maf"
printf 'UNKNOWN\t0\ttest\thg19\tchr9\t9876543\t9876543\t+\tMissense_Mutation\tSNP\tA\tG\tG\t.\t.\tSAMPLE_001\n'  >> "$TMPDIR/s001_extra.maf"

filter_maf     "$TMPDIR/s001_extra.maf" "$TSV_001" "$TMPDIR/s001_extra_filtered.txt"
passthrough_maf "$TMPDIR/s001_extra.maf"            "$TMPDIR/s001_extra_passthrough.txt"

filtered_n=$(wc -l < "$TMPDIR/s001_extra_filtered.txt")
passthrough_n=$(wc -l < "$TMPDIR/s001_extra_passthrough.txt")

check "non-TSV variants: filter drops them (only 2 TSV matches + header)" "3" "$filtered_n"
check "non-TSV variants: passthrough keeps all 4 mutations + header" "5" "$passthrough_n"
check "non-TSV variants: passthrough > filter" "true" "$([ "$passthrough_n" -gt "$filtered_n" ] && echo true || echo false)"

# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "Results: $PASS passed, $FAIL failed"
[ "$FAIL" -eq 0 ]
