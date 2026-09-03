#!/usr/bin/env bash
# create_test_data.sh — Generate minimal test fixtures for oncoanalyser nf-test tests.
#
# Requirements: bgzip, tabix (from htslib)
# Usage: bash bin/create_test_data.sh  (run from oncoanalyser/ directory)
#
# Creates:
#   assets/test_data/oncoanalyser_output/TEST/T{1,2}/  — pipeline input fixtures
#   assets/test_data/precomputed/                       — pre-computed aggregate/ML fixtures

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="$(cd "${SCRIPT_DIR}/.." && pwd)/assets/test_data"

echo "=== Creating oncoanalyser test fixtures ==="
echo "Base directory: ${BASE}"

# Check dependencies
for cmd in bgzip tabix; do
    if ! command -v "$cmd" &>/dev/null; then
        echo "ERROR: '$cmd' not found. Install htslib (conda install -c bioconda htslib)."
        exit 1
    fi
done

# ──────────────────────────────────────────────────────────────────────────────
# Helper: write + bgzip + tabix a VCF
# ──────────────────────────────────────────────────────────────────────────────
make_vcf() {
    local path="$1"
    local content="$2"
    mkdir -p "$(dirname "$path")"
    printf '%s' "$content" > "${path%.gz}"
    bgzip -f "${path%.gz}"
    tabix -p vcf "$path"
}

# ──────────────────────────────────────────────────────────────────────────────
# VCF templates
# Using chr21 positions that fall within genes in the annotation file:
#   PRDM15: chr21:41798225-41879482
#   PFKL:   chr21:44300051-44327376
#   HUNK:   chr21:31873020-32044633
#   ERG (not in annotation but common chr21 gene): chr21:38380000
# ──────────────────────────────────────────────────────────────────────────────

# Shared VCF header for somatic/germline SNV files (2 sample columns: normal T1-N, tumor T1)
snv_header_t1='##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=LOW_QUAL,Description="Low quality">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##contig=<ID=chr21,length=46709983>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	T1-N	T1
'

snv_header_t2='##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=LOW_QUAL,Description="Low quality">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##contig=<ID=chr21,length=46709983>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	T2-N	T2
'

# BND VCF header for ESVEE SVs
bnd_header_t1='##fileformat=VCFv4.1
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">
##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakend">
##FORMAT=<ID=BVF,Number=1,Type=Integer,Description="Breakend variant fragments">
##FORMAT=<ID=VF,Number=1,Type=Integer,Description="Variant fragments">
##contig=<ID=chr21,length=46709983>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	T1
'

bnd_header_t2='##fileformat=VCFv4.1
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">
##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakend">
##FORMAT=<ID=BVF,Number=1,Type=Integer,Description="Breakend variant fragments">
##FORMAT=<ID=VF,Number=1,Type=Integer,Description="Variant fragments">
##contig=<ID=chr21,length=46709983>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	T2
'

echo ""
echo "--- Creating T1 fixtures ---"

# pave somatic VCF (T1)
make_vcf "${BASE}/oncoanalyser_output/TEST/T1/pave/somatic/T1-T.pave.somatic.vcf.gz" \
"${snv_header_t1}chr21	41800000	.	C	T	60	PASS	DP=100	GT:AD:DP	0/0:50,0:50	0/1:40,30:70
"

# pave germline VCF (T1) — must have sample column T1 so FILTER_GERMLINE_DNA can exclude it
make_vcf "${BASE}/oncoanalyser_output/TEST/T1/pave/germline/T1-T.pave.germline.vcf.gz" \
"${snv_header_t1}chr21	44300100	.	A	G	80	PASS	DP=120	GT:AD:DP	0/1:55,10:65	0/0:60,0:60
"

# sage append VCF (RNA-append somatic, T1)
make_vcf "${BASE}/oncoanalyser_output/TEST/T1/sage_append/T1-T/T1-T.sage.append.vcf.gz" \
"${snv_header_t1}chr21	41800000	.	C	T	60	PASS	DP=80	GT:AD:DP	0/0:40,0:40	0/1:30,25:55
"

# esvee BND VCF (T1) — one PASS BND pair
make_vcf "${BASE}/oncoanalyser_output/TEST/T1/esvee/T1-T.esvee.unfiltered.vcf.gz" \
"${bnd_header_t1}chr21	41800000	BND_T1_1	A	A[chr21:44300000[	.	PASS	SVTYPE=BND;MATEID=BND_T1_2	BVF:VF	3:5
chr21	44300000	BND_T1_2	G	]chr21:41800000]G	.	PASS	SVTYPE=BND;MATEID=BND_T1_1	BVF:VF	3:5
"

# PURPLE somatic CNV segments (T1)
mkdir -p "${BASE}/oncoanalyser_output/TEST/T1/purple"
cat > "${BASE}/oncoanalyser_output/TEST/T1/purple/T1-T.purple.cnv.somatic.tsv" << 'EOF'
chromosome	start	end	copyNumber	bafCount	observedBAF	baf	depthWindowCount	gcContent	method	minStart	maxStart
chr21	10000000	20000000	2.1	150	0.52	0.51	200	0.42	PCNI	9990000	10010000
chr21	20000000	30000000	3.8	200	0.48	0.49	300	0.43	PCNI	19990000	20010000
chr21	30000000	41000000	1.2	80	0.55	0.54	150	0.41	PCNI	29990000	30010000
EOF

# PURPLE gene-level CNV (T1) — needs gene + minCopyNumber + extras (script uses cut -f1-19)
cat > "${BASE}/oncoanalyser_output/TEST/T1/purple/T1-T.purple.cnv.gene.tsv" << 'EOF'
gene	chromosome	start	end	minCopyNumber	maxCopyNumber	unused1	unused2	unused3	unused4	unused5	unused6	unused7	unused8	unused9	unused10	unused11	unused12	unused13
PRDM15	chr21	41798225	41879482	3.5	4.1	0	0	0	0	0	0	0	0	0	0	0	0	0
PFKL	chr21	44300051	44327376	1.2	1.8	0	0	0	0	0	0	0	0	0	0	0	0	0
HUNK	chr21	31873020	32044633	0.3	0.5	0	0	0	0	0	0	0	0	0	0	0	0	0
EOF

# Isofox gene data TSV (T1)
mkdir -p "${BASE}/oncoanalyser_output/TEST/T1/isofox"
cat > "${BASE}/oncoanalyser_output/TEST/T1/isofox/T1-T.isf.gene_data.tsv" << 'EOF'
GeneId	GeneName	rawCounts	fragmentsPerKbLength	medianGCContent	TPM	AdjustedTPM	spliceJunctionCount
ENSG00000141956	PRDM15	1200	120.5	0.55	145.82	140.10	85
ENSG00000141959	PFKL	800	80.2	0.52	23.45	22.80	40
ENSG00000142149	HUNK	300	30.1	0.48	8.12	7.90	15
ENSG00000142156	COL6A1	2000	200.0	0.56	312.50	308.20	120
EOF

# Isofox pass fusions (T1) — PRDM15--PFKL, a TMPRSS2--ERG deletion, and three
# rows whose Name has an empty half (unannotated breakend) that must be filtered
cat > "${BASE}/oncoanalyser_output/TEST/T1/isofox/T1-T.isf.pass_fusions.tsv" << 'EOF'
Name	KnownType	ChromosomeUp	ChromosomeDown	PositionUp	PositionDown	OrientationUp	OrientationDown	JunctionTypeUp	JunctionTypeDown	TranscriptUp	TranscriptDown	ExonUp	ExonDown	SvType	SplitFrags	RealignedFrags	DiscordantFrags	DepthUp	DepthDown	MaxAnchorLengthUp	MaxAnchorLengthDown	CohortFrequency
PRDM15_PFKL	KNOWN_PAIR	chr21	chr21	41800000	44300000	1	-1	KNOWN	KNOWN	ENST00000398642	ENST00000349048	2	3	DEL	8	0	3	24	26	75	80	0
TMPRSS2_ERG	NONE	chr21	chr21	41500000	38400000	-1	1	CANONICAL	KNOWN	ENST00000332149	ENST00000442448	1	4	DEL	4	1	2	18	15	90	85	0.01
_ORPHAN3P	NONE	chr21	chr21	42000000	42500000	1	-1	CANONICAL	KNOWN		ENST00000399329	0	2	DEL	4	0	0	12	11	70	100	0
ORPHAN5P_	NONE	chr21	chr21	43000000	43500000	-1	1	KNOWN	CANONICAL	ENST00000399330		3	0	DEL	2	0	1	15	9	110	40	0.05
_	NONE	chr21	chr21	43800000	43900000	1	-1	CANONICAL	CANONICAL			0	0	DEL	3	0	0	10	10	120	120	0
EOF

echo "--- Creating T2 fixtures ---"

# pave somatic VCF (T2)
make_vcf "${BASE}/oncoanalyser_output/TEST/T2/pave/somatic/T2-T.pave.somatic.vcf.gz" \
"${snv_header_t2}chr21	41800500	.	G	A	55	PASS	DP=90	GT:AD:DP	0/0:45,0:45	0/1:35,28:63
"

# pave germline VCF (T2)
make_vcf "${BASE}/oncoanalyser_output/TEST/T2/pave/germline/T2-T.pave.germline.vcf.gz" \
"${snv_header_t2}chr21	44300200	.	T	C	75	PASS	DP=110	GT:AD:DP	0/1:50,12:62	0/0:55,0:55
"

# sage append VCF (T2)
make_vcf "${BASE}/oncoanalyser_output/TEST/T2/sage_append/T2-T/T2-T.sage.append.vcf.gz" \
"${snv_header_t2}chr21	41800500	.	G	A	55	PASS	DP=70	GT:AD:DP	0/0:35,0:35	0/1:28,22:50
"

# esvee BND VCF (T2)
make_vcf "${BASE}/oncoanalyser_output/TEST/T2/esvee/T2-T.esvee.unfiltered.vcf.gz" \
"${bnd_header_t2}chr21	41800500	BND_T2_1	C	C[chr21:44300200[	.	PASS	SVTYPE=BND;MATEID=BND_T2_2	BVF:VF	4:6
chr21	44300200	BND_T2_2	A	]chr21:41800500]A	.	PASS	SVTYPE=BND;MATEID=BND_T2_1	BVF:VF	4:6
"

# PURPLE somatic CNV segments (T2)
mkdir -p "${BASE}/oncoanalyser_output/TEST/T2/purple"
cat > "${BASE}/oncoanalyser_output/TEST/T2/purple/T2-T.purple.cnv.somatic.tsv" << 'EOF'
chromosome	start	end	copyNumber	bafCount	observedBAF	baf	depthWindowCount	gcContent	method	minStart	maxStart
chr21	10000000	22000000	1.8	130	0.53	0.52	180	0.42	PCNI	9990000	10010000
chr21	22000000	35000000	4.2	220	0.47	0.48	320	0.43	PCNI	21990000	22010000
chr21	35000000	42000000	2.0	100	0.50	0.50	160	0.41	PCNI	34990000	35010000
EOF

cat > "${BASE}/oncoanalyser_output/TEST/T2/purple/T2-T.purple.cnv.gene.tsv" << 'EOF'
gene	chromosome	start	end	minCopyNumber	maxCopyNumber	unused1	unused2	unused3	unused4	unused5	unused6	unused7	unused8	unused9	unused10	unused11	unused12	unused13
PRDM15	chr21	41798225	41879482	1.8	2.2	0	0	0	0	0	0	0	0	0	0	0	0	0
PFKL	chr21	44300051	44327376	4.5	5.0	0	0	0	0	0	0	0	0	0	0	0	0	0
HUNK	chr21	31873020	32044633	0.4	0.6	0	0	0	0	0	0	0	0	0	0	0	0	0
EOF

mkdir -p "${BASE}/oncoanalyser_output/TEST/T2/isofox"
cat > "${BASE}/oncoanalyser_output/TEST/T2/isofox/T2-T.isf.gene_data.tsv" << 'EOF'
GeneId	GeneName	rawCounts	fragmentsPerKbLength	medianGCContent	TPM	AdjustedTPM	spliceJunctionCount
ENSG00000141956	PRDM15	900	90.5	0.54	110.20	106.50	70
ENSG00000141959	PFKL	1500	150.0	0.53	45.80	44.20	65
ENSG00000142149	HUNK	200	20.0	0.47	5.30	5.10	10
ENSG00000142156	COL6A1	1800	180.0	0.55	280.40	276.80	105
EOF

# Isofox pass fusions (T2) — two rows on identical coordinates, which dedupe to one
cat > "${BASE}/oncoanalyser_output/TEST/T2/isofox/T2-T.isf.pass_fusions.tsv" << 'EOF'
Name	KnownType	ChromosomeUp	ChromosomeDown	PositionUp	PositionDown	OrientationUp	OrientationDown	JunctionTypeUp	JunctionTypeDown	TranscriptUp	TranscriptDown	ExonUp	ExonDown	SvType	SplitFrags	RealignedFrags	DiscordantFrags	DepthUp	DepthDown	MaxAnchorLengthUp	MaxAnchorLengthDown	CohortFrequency
PRDM15_PFKL	KNOWN_PAIR	chr21	chr21	41800500	44300200	1	-1	KNOWN	KNOWN	ENST00000398642	ENST00000349048	2	3	DEL	6	0	2	20	22	75	80	0
PRDM15_PFKL	KNOWN_PAIR	chr21	chr21	41800500	44300200	1	-1	CANONICAL	KNOWN	ENST00000398643	ENST00000349048	2	3	DEL	1	0	0	20	22	60	65	0
EOF

# ──────────────────────────────────────────────────────────────────────────────
# Precomputed fixture files for aggregate_output and ML subworkflow tests
# ──────────────────────────────────────────────────────────────────────────────
echo ""
echo "--- Creating precomputed fixtures ---"
mkdir -p "${BASE}/precomputed"

# SEG files (output of PURPLE_CNV_TO_CBIOPORTAL)
cat > "${BASE}/precomputed/T1_data_cna_hg38.seg" << 'EOF'
ID	chrom	loc.start	loc.end	num.mark	seg.mean
T1	21	10000000	20000000	150	0.0704
T1	21	20000000	30000000	200	0.9236
T1	21	30000000	41000000	80	-0.7370
EOF

cat > "${BASE}/precomputed/T2_data_cna_hg38.seg" << 'EOF'
ID	chrom	loc.start	loc.end	num.mark	seg.mean
T2	21	10000000	22000000	130	-0.1520
T2	21	22000000	35000000	220	1.0700
T2	21	35000000	42000000	100	0.0000
EOF

# CNA long files (output of PURPLE_CNV_TO_CBIOPORTAL — discrete CNA)
cat > "${BASE}/precomputed/T1_data_cna_long.txt" << 'EOF'
Hugo_Symbol	Entrez_Gene_Id	Sample_Id	Value
PRDM15	100001	T1	1
PFKL	100002	T1	-1
HUNK	100003	T1	-2
EOF

cat > "${BASE}/precomputed/T2_data_cna_long.txt" << 'EOF'
Hugo_Symbol	Entrez_Gene_Id	Sample_Id	Value
PRDM15	100001	T2	0
PFKL	100002	T2	2
HUNK	100003	T2	-2
EOF

# Group-merged CNA long (used by FORMAT_ML_CNV — must be named exactly data_cna_long.txt)
cat > "${BASE}/precomputed/data_cna_long.txt" << 'EOF'
Hugo_Symbol	Entrez_Gene_Id	Sample_Id	Value
PRDM15	100001	T1	1
PFKL	100002	T1	-1
HUNK	100003	T1	-2
PRDM15	100001	T2	0
PFKL	100002	T2	2
HUNK	100003	T2	-2
EOF

# SV files (output of ESVEE_SV_TO_CBIOPORTAL)
cat > "${BASE}/precomputed/T1.data_sv.txt" << 'EOF'
Sample_Id	SV_Status	Site1_Hugo_Symbol	Site1_Ensembl_Transcript_Id	Site1_Region_Number	Site1_Region	Site2_Hugo_Symbol	Site2_Ensembl_Transcript_Id	Site2_Region_Number	Site2_Region	Site2_Effect_On_Frame	NCBI_Build	Class	DNA_Support	RNA_Support	Tumor_Variant_Count	Connection_Type	Breakpoint_Type	Event_Info	Annotation	Site1_Chromosome	Site1_Position	Site2_Chromosome	Site2_Position	Tumor_Split_Read_Count	Tumor_Paired_End_Read_Count
T1	SOMATIC	PRDM15		.	.	PFKL		.	.		GRCh38	INVERSION	Yes	No	3	5to3	PRECISE	DNA SV: PRDM15--PFKL		21	41800000	21	44300000	3	5
EOF

cat > "${BASE}/precomputed/T2.data_sv.txt" << 'EOF'
Sample_Id	SV_Status	Site1_Hugo_Symbol	Site1_Ensembl_Transcript_Id	Site1_Region_Number	Site1_Region	Site2_Hugo_Symbol	Site2_Ensembl_Transcript_Id	Site2_Region_Number	Site2_Region	Site2_Effect_On_Frame	NCBI_Build	Class	DNA_Support	RNA_Support	Tumor_Variant_Count	Connection_Type	Breakpoint_Type	Event_Info	Annotation	Site1_Chromosome	Site1_Position	Site2_Chromosome	Site2_Position	Tumor_Split_Read_Count	Tumor_Paired_End_Read_Count
T2	SOMATIC	PRDM15		.	.	PFKL		.	.		GRCh38	INVERSION	Yes	No	4	5to3	PRECISE	DNA SV: PRDM15--PFKL		21	41800500	21	44300200	4	6
EOF

# Group-merged SV (used by FORMAT_PROCESS_ML_SV — named data_sv.txt)
cat > "${BASE}/precomputed/data_sv.txt" << 'EOF'
Sample_Id	SV_Status	Site1_Hugo_Symbol	Site1_Ensembl_Transcript_Id	Site1_Region_Number	Site1_Region	Site2_Hugo_Symbol	Site2_Ensembl_Transcript_Id	Site2_Region_Number	Site2_Region	Site2_Effect_On_Frame	NCBI_Build	Class	DNA_Support	RNA_Support	Tumor_Variant_Count	Connection_Type	Breakpoint_Type	Event_Info	Annotation	Site1_Chromosome	Site1_Position	Site2_Chromosome	Site2_Position	Tumor_Split_Read_Count	Tumor_Paired_End_Read_Count
T1	SOMATIC	PRDM15		.	.	PFKL		.	.		GRCh38	INVERSION	Yes	No	3	5to3	PRECISE	DNA SV: PRDM15--PFKL		21	41800000	21	44300000	3	5
T2	SOMATIC	PRDM15		.	.	PFKL		.	.		GRCh38	INVERSION	Yes	No	4	5to3	PRECISE	DNA SV: PRDM15--PFKL		21	41800500	21	44300200	4	6
EOF

# TPM files (output of ISOFOX_EXPRESSION_TO_CBIOPORTAL — per-sample TPM)
cat > "${BASE}/precomputed/T1.tpm.tsv" << 'EOF'
Hugo_Symbol	Entrez_Gene_Id	T1
PRDM15	100001	145.82
PFKL	100002	23.45
HUNK	100003	8.12
COL6A1	100004	312.50
EOF

cat > "${BASE}/precomputed/T2.tpm.tsv" << 'EOF'
Hugo_Symbol	Entrez_Gene_Id	T2
PRDM15	100001	110.20
PFKL	100002	45.80
HUNK	100003	5.30
COL6A1	100004	280.40
EOF

# Group-merged expression (used by FORMAT_ML_EXPRESSION — must be named data_expression.txt)
cat > "${BASE}/precomputed/data_expression.txt" << 'EOF'
Hugo_Symbol	Entrez_Gene_Id	T1	T2
PRDM15	100001	145.82	110.20
PFKL	100002	23.45	45.80
HUNK	100003	8.12	5.30
COL6A1	100004	312.50	280.40
EOF

# MAF files (output of CONVERT_CPSR_TO_MAF — somatic+germline combined)
# MAF format: line 1 = version comment, line 2 = column header, line 3+ = data
# aggregate_output uses skip: 2, so lines 1+2 are skipped when merging
MAF_HEADER="Hugo_Symbol	Entrez_Gene_Id	Center	NCBI_Build	Chromosome	Start_Position	End_Position	Strand	Variant_Classification	Variant_Type	Reference_Allele	Tumor_Seq_Allele1	Tumor_Seq_Allele2	dbSNP_RS	dbSNP_Val_Status	Tumor_Sample_Barcode	Matched_Norm_Sample_Barcode	Match_Norm_Seq_Allele1	Match_Norm_Seq_Allele2	Tumor_Validation_Allele1	Tumor_Validation_Allele2	Match_Norm_Validation_Allele1	Match_Norm_Validation_Allele2	Verification_Status	Validation_Status	Mutation_Status	Sequencing_Phase	Sequence_Source	Validation_Method	Score	BAM_File	Sequencer	Tumor_Sample_UUID	Matched_Norm_Sample_UUID	HGVSc	HGVSp	HGVSp_Short	Transcript_ID	Exon_Number	t_depth	t_ref_count	t_alt_count	n_depth	n_ref_count	n_alt_count	all_effects	Allele	Gene	Feature	Feature_type	Consequence	cDNA_position	CDS_position	Protein_position	Amino_acids	Codons	Existing_variation	ALLELE_NUM	DISTANCE	TRANSCRIPT_STRAND	SYMBOL	SYMBOL_SOURCE	HGNC_ID	BIOTYPE	CANONICAL	CCDS	ENSP	SWISSPROT	TREMBL	UNIPARC	RefSeq	SIFT	PolyPhen	EXON	INTRON	DOMAINS	AF	AFR_AF	AMR_AF	ASIAN_AF	ASN_AF	EAS_AF	EUR_AF	SAS_AF	AA_AF	EA_AF	gnomAD_AF	gnomAD_AFR_AF	gnomAD_AMR_AF	gnomAD_ASJ_AF	gnomAD_EAS_AF	gnomAD_FIN_AF	gnomAD_NFE_AF	gnomAD_OTH_AF	gnomAD_SAS_AF	MAX_AF	MAX_AF_POPS	CLIN_SIG	SOMATIC	PUBMED	MOTIF_NAME	MOTIF_POS	HIGH_INF_POS	MOTIF_SCORE_CHANGE	IMPACT	PICK	VARIANT_CLASS	TSL	HGVS_OFFSET	PHENO	MINIMISED	ExAC_AF	ExAC_AF_AFR	ExAC_AF_AMR	ExAC_AF_EAS	ExAC_AF_FIN	ExAC_AF_NFE	ExAC_AF_OTH	ExAC_AF_SAS	GENE_PHENO	FILTER	CONTEXT	src_vcf_id	tumor_bam_uuid	normal_bam_uuid	case_id	GDC_Validation_Status	GDC_Valid_Somatic	vcf_region	vcf_info	vcf_format	vcf_tumor_gt	Mutation_Status_extra"

printf '#version 2.4\n' > "${BASE}/precomputed/T1.somatic_rna_germline.maf"
printf '%s\n' "$MAF_HEADER" >> "${BASE}/precomputed/T1.somatic_rna_germline.maf"
printf 'PRDM15\t100001\t.\tGRCh38\tchr21\t41800000\t41800000\t+\tMissense_Mutation\tSNP\tC\tC\tT\t.\t.\tT1\tT1-N\tC\tC\t.\t.\t.\t.\t.\t.\tSOMATC\t.\t.\t.\t.\t.\t.\t.\t.\tc.123C>T\tp.Arg41Cys\tp.R41C\tENST00000399231\t2/10\t70\t40\t30\t65\t50\t0\t.\tT\tENSG00000141956\tENST00000399231\tTranscript\tMissense_Variant\t123/1200\t123/1200\t41/400\tR/C\tcGt/tGt\t.\t1\t.\t-1\tPRDM15\tHGNC\tHGNC:13999\tprotein_coding\tYES\t.\t.\t.\t.\t.\t.\t.\t.\t.\t2/10\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\tMODERATE\t1\tSNV\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\tPASS\tACTCG\t.\t.\t.\t.\t.\t.\tchr21:41800000-41800000\t.\t.\t.\tSOMATC\n' >> "${BASE}/precomputed/T1.somatic_rna_germline.maf"

printf '#version 2.4\n' > "${BASE}/precomputed/T2.somatic_rna_germline.maf"
printf '%s\n' "$MAF_HEADER" >> "${BASE}/precomputed/T2.somatic_rna_germline.maf"
printf 'PFKL\t100002\t.\tGRCh38\tchr21\t44300200\t44300200\t+\tMissense_Mutation\tSNP\tT\tT\tC\t.\t.\tT2\tT2-N\tT\tT\t.\t.\t.\t.\t.\t.\tSOMATC\t.\t.\t.\t.\t.\t.\t.\t.\tc.456T>C\tp.Phe152Leu\tp.F152L\tENST00000399500\t5/15\t63\t35\t28\t62\t50\t0\t.\tC\tENSG00000141959\tENST00000399500\tTranscript\tMissense_Variant\t456/1500\t456/1500\t152/500\tF/L\ttTt/tCt\t.\t1\t.\t1\tPFKL\tHGNC\tHGNC:8876\tprotein_coding\tYES\t.\t.\t.\t.\t.\t.\t.\t.\t.\t5/15\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\tMODERATE\t1\tSNV\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\tPASS\tGATCG\t.\t.\t.\t.\t.\t.\tchr21:44300200-44300200\t.\t.\t.\tSOMATC\n' >> "${BASE}/precomputed/T2.somatic_rna_germline.maf"

# Group-merged mutations (used by FORMAT_ML_MUTATION — named data_mutations_dna_rna_germline.txt)
# skip:2 merges: line 1 (#version) and line 2 (header) are skipped in all but first file
{
    cat "${BASE}/precomputed/T1.somatic_rna_germline.maf"
    tail -n +3 "${BASE}/precomputed/T2.somatic_rna_germline.maf"
} > "${BASE}/precomputed/data_mutations_dna_rna_germline.txt"

echo ""
echo "=== Done! ==="
echo ""
echo "Created:"
find "${BASE}/oncoanalyser_output" -type f | sort | sed "s|${BASE}/||"
echo ""
find "${BASE}/precomputed" -type f | sort | sed "s|${BASE}/||"
echo ""
echo "Next step: commit these files and run:"
echo "  nf-test test tests/subworkflows/genomic_cnv.nf.test --profile test,apptainer"
