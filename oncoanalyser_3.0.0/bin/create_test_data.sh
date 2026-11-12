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

# PURPLE somatic VCF header: the SNV header plus the INFO fields PURPLE adds
purple_somatic_header() {
    printf '%s\n' \
'##fileformat=VCFv4.2' \
'##FILTER=<ID=PASS,Description="All filters passed">' \
'##FILTER=<ID=LOW_QUAL,Description="Low quality">' \
'##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">' \
'##INFO=<ID=PURPLE_AF,Number=1,Type=Float,Description="Purity adjusted allele frequency">' \
'##INFO=<ID=PURPLE_CN,Number=1,Type=Float,Description="Purity adjusted copy number surrounding variant location">' \
'##INFO=<ID=PURPLE_VCN,Number=1,Type=Float,Description="Purity adjusted variant copy number">' \
'##INFO=<ID=PURPLE_MACN,Number=1,Type=Float,Description="Purity adjusted minor allele copy number">' \
'##INFO=<ID=BIALLELIC,Number=0,Type=Flag,Description="Variant is biallelic">' \
'##INFO=<ID=SUBCL,Number=1,Type=Float,Description="Subclonal likelihood">' \
'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' \
'##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">' \
'##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">' \
'##contig=<ID=chr21,length=46709983>' \
"#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	$1-N	$1"
}

# PURPLE germline VCF header: PURPLE calls GT where it can (else ./.) and adds LOW_VAF / LOW_TUMOR_VCN
purple_germline_header() {
    printf '%s\n' \
'##fileformat=VCFv4.2' \
'##FILTER=<ID=PASS,Description="All filters passed">' \
'##FILTER=<ID=LOW_QUAL,Description="Low quality">' \
'##FILTER=<ID=LOW_VAF,Description="Germline variant has very low allelic frequency">' \
'##FILTER=<ID=LOW_TUMOR_VCN,Description="Germline variant has very low tumor variant copy number">' \
'##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">' \
'##INFO=<ID=PURPLE_CN,Number=1,Type=Float,Description="Purity adjusted copy number surrounding variant location">' \
'##INFO=<ID=PURPLE_VCN,Number=1,Type=Float,Description="Purity adjusted variant copy number">' \
'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' \
'##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">' \
'##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allelic frequency calculated from read context counts">' \
'##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">' \
'##contig=<ID=chr21,length=46709983>' \
"#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	$1-N	$1"
}

# PURPLE somatic SV VCF header: ESVEE's records + PURPLE's INFO, normal then tumour column
purple_sv_header() {
printf '%s\n' \
'##fileformat=VCFv4.2' \
'##FILTER=<ID=PASS,Description="All filters passed">' \
'##FILTER=<ID=PON,Description="Found in panel of normals">' \
'##FILTER=<ID=INFERRED,Description="Breakend inferred from copy number transition">' \
'##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">' \
'##INFO=<ID=MATEID,Number=1,Type=String,Description="Id of other breakend in junction">' \
'##INFO=<ID=INFERRED,Number=0,Type=Flag,Description="Breakend inferred from copy number transition">' \
'##INFO=<ID=TINC_RECOVERED,Number=0,Type=Flag,Description="Recovered from germline by tumor-in-normal contamination">' \
'##INFO=<ID=PURPLE_AF,Number=.,Type=Float,Description="Purity adjusted allele frequency at each breakend">' \
'##INFO=<ID=PURPLE_CN,Number=.,Type=Float,Description="Purity adjusted copy number at each breakend">' \
'##INFO=<ID=PURPLE_CN_CHANGE,Number=.,Type=Float,Description="Purity adjusted change in copy number at each breakend">' \
'##INFO=<ID=PURPLE_JCN,Number=1,Type=Float,Description="Purity adjusted junction copy number">' \
'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' \
'##FORMAT=<ID=DF,Number=1,Type=Integer,Description="Count of discordant fragments with a read either side of the breakend">' \
'##FORMAT=<ID=SF,Number=1,Type=Integer,Description="Count of fragments supporting the breakend with a read overlapping the breakend">' \
'##FORMAT=<ID=VF,Number=1,Type=Integer,Description="Total variant fragments supporting the breakend">' \
'##contig=<ID=chr21,length=46709983>' \
"#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	$1-N	$1-T"
}

echo ""
echo "--- Creating T1 fixtures ---"

# PURPLE somatic VCF (T1) — PAVE's records plus PURPLE's INFO. Adjacent SNV pairs (300/301,
# 500/501, 1000/1001) are DBSs for SigProfiler DBS; 800 is LOW_QUAL.
make_vcf "${BASE}/oncoanalyser_output/TEST/T1/purple/T1-T.purple.somatic.vcf.gz" \
"$(purple_somatic_header T1)
chr21	10000100	.	C	T	60	PASS	DP=100;PURPLE_AF=0.40;PURPLE_CN=2;PURPLE_VCN=0.80;PURPLE_MACN=1;BIALLELIC;SUBCL=0	GT:AD:DP	0/0:50,0:50	0/1:40,30:70
chr21	10000200	.	G	A	55	PASS	DP=90;PURPLE_AF=0.41;PURPLE_CN=2;PURPLE_VCN=0.81;PURPLE_MACN=1;SUBCL=0	GT:AD:DP	0/0:45,0:45	0/1:35,25:60
chr21	10000300	.	A	G	58	PASS	DP=95;PURPLE_AF=0.42;PURPLE_CN=2;PURPLE_VCN=0.82;PURPLE_MACN=1;SUBCL=0	GT:AD:DP	0/0:48,0:48	0/1:38,28:66
chr21	10000301	.	C	T	62	PASS	DP=88;PURPLE_AF=0.43;PURPLE_CN=2;PURPLE_VCN=0.83;PURPLE_MACN=1;SUBCL=0	GT:AD:DP	0/0:44,0:44	0/1:36,26:62
chr21	10000500	.	T	C	50	PASS	DP=80;PURPLE_AF=0.44;PURPLE_CN=2;PURPLE_VCN=0.84;PURPLE_MACN=1;SUBCL=0	GT:AD:DP	0/0:40,0:40	0/1:32,22:54
chr21	10000501	.	G	A	48	PASS	DP=85;PURPLE_AF=0.45;PURPLE_CN=2;PURPLE_VCN=0.85;PURPLE_MACN=1;SUBCL=0	GT:AD:DP	0/0:42,0:42	0/1:34,24:58
chr21	10000700	.	C	A	65	PASS	DP=110;PURPLE_AF=0.46;PURPLE_CN=2;PURPLE_VCN=0.86;PURPLE_MACN=1;SUBCL=0	GT:AD:DP	0/0:55,0:55	0/1:44,32:76
chr21	10000800	.	G	T	45	LOW_QUAL	DP=50;PURPLE_AF=0.47;PURPLE_CN=2;PURPLE_VCN=0.87;PURPLE_MACN=1;SUBCL=0	GT:AD:DP	0/0:25,0:25	0/1:20,10:30
chr21	10001000	.	T	A	52	PASS	DP=92;PURPLE_AF=0.48;PURPLE_CN=2;PURPLE_VCN=0.88;PURPLE_MACN=1;SUBCL=0	GT:AD:DP	0/0:46,0:46	0/1:36,26:62
chr21	10001001	.	C	G	54	PASS	DP=88;PURPLE_AF=0.49;PURPLE_CN=2;PURPLE_VCN=0.89;PURPLE_MACN=1;SUBCL=0	GT:AD:DP	0/0:44,0:44	0/1:35,25:60
"

# PURPLE germline VCF (T1) — must have sample column T1 so FILTER_GERMLINE_DNA can exclude it.
# 100-500 are ./.; the normal's DP/AF cover each genotype rule and its boundaries:
#   100 hom, 150 hom at AF 0.9, 200 het, 250 het at AF 0.2, 300 ref (AF < 0.2),
#   350 stays ./. (DP 9), 400 ref at low depth, 500 LOW_QUAL (dropped).
# 600 PURPLE called 0/1 (kept, though AF 0.95 would infer 1/1), 700 LOW_TUMOR_VCN (kept, inferred 0/1),
# 800 LOW_VAF and 900 LOW_TUMOR_VCN;LOW_VAF (both dropped).
make_vcf "${BASE}/oncoanalyser_output/TEST/T1/purple/T1-T.purple.germline.vcf.gz" \
"$(purple_germline_header T1)
chr21	44300100	.	A	G	80	PASS	DP=80	GT:AD:AF:DP	./.:2,38:0.95:40	./.:2,38:0.95:40
chr21	44300150	.	C	T	80	PASS	DP=60	GT:AD:AF:DP	./.:3,27:0.9:30	./.:3,27:0.9:30
chr21	44300200	.	G	A	80	PASS	DP=120	GT:AD:AF:DP	./.:30,30:0.5:60	./.:30,30:0.5:60
chr21	44300250	.	T	C	80	PASS	DP=100	GT:AD:AF:DP	./.:40,10:0.2:50	./.:40,10:0.2:50
chr21	44300300	.	A	C	80	PASS	DP=160	GT:AD:AF:DP	./.:68,12:0.15:80	./.:60,20:0.25:80
chr21	44300350	.	G	T	80	PASS	DP=18	GT:AD:AF:DP	./.:4,5:0.5:9	./.:4,5:0.5:9
chr21	44300400	.	C	G	80	PASS	DP=10	GT:AD:AF:DP	./.:4,1:0.1:5	./.:4,1:0.1:5
chr21	44300500	.	T	A	40	LOW_QUAL	DP=100	GT:AD:AF:DP	./.:25,25:0.5:50	./.:25,25:0.5:50
chr21	44300600	.	A	T	80	PASS	DP=80;PURPLE_CN=2;PURPLE_VCN=1	GT:AD:AF:DP	0/1:2,38:0.95:40	./.:2,38:0.95:40
chr21	44300700	.	G	C	80	LOW_TUMOR_VCN	DP=100;PURPLE_CN=1;PURPLE_VCN=0.1	GT:AD:AF:DP	./.:25,25:0.5:50	./.:48,2:0.04:50
chr21	44300800	.	C	A	80	LOW_VAF	DP=100	GT:AD:AF:DP	./.:25,25:0.5:50	./.:25,25:0.5:50
chr21	44300900	.	T	G	80	LOW_TUMOR_VCN;LOW_VAF	DP=100	GT:AD:AF:DP	./.:25,25:0.5:50	./.:25,25:0.5:50
"

# sage append VCF (RNA-append somatic, T1)
make_vcf "${BASE}/oncoanalyser_output/TEST/T1/sage_append/T1-T/T1-T.sage.append.vcf.gz" \
"${snv_header_t1}chr21	41800000	.	C	T	60	PASS	DP=80	GT:AD:DP	0/0:40,0:40	0/1:30,25:55
"

# PURPLE SV VCF (T1). PURPLE_AF/CN/CN_CHANGE are [this breakend, mate], so the mate's list is reversed.
# Kept: the PRDM15-PFKL PASS pair and the TINC_RECOVERED SOD1-HUNK pair. Dropped: IFNAR1 INFERRED breakend,
# SIK1 single breakend (no MATEID), COL6A1-COL6A2 PON pair. The normal's counts differ from the tumour's.
make_vcf "${BASE}/oncoanalyser_output/TEST/T1/purple/T1-T.purple.sv.vcf.gz" \
"$(purple_sv_header T1)
chr21	31660000	SV_T1_3o	T	T]chr21:31900000]	800	PASS	SVTYPE=BND;MATEID=SV_T1_3h;TINC_RECOVERED;PURPLE_AF=0.12,0.15;PURPLE_CN=2.1,2.3;PURPLE_CN_CHANGE=0.3,0.35;PURPLE_JCN=0.32	GT:DF:SF:VF	./.:1:1:2	./.:2:3:5
chr21	31900000	SV_T1_3h	C	C]chr21:31660000]	800	PASS	SVTYPE=BND;MATEID=SV_T1_3o;TINC_RECOVERED;PURPLE_AF=0.15,0.12;PURPLE_CN=2.3,2.1;PURPLE_CN_CHANGE=0.35,0.3;PURPLE_JCN=0.32	GT:DF:SF:VF	./.:1:1:2	./.:2:3:5
chr21	33340000	purple_1	G	.G	.	INFERRED	SVTYPE=BND;INFERRED;PURPLE_AF=0.4;PURPLE_CN=3;PURPLE_CN_CHANGE=1	GT:DF:SF:VF	./.:0:0:0	./.:0:0:0
chr21	41800000	SV_T1_1o	A	A[chr21:44300100[	1500	PASS	SVTYPE=BND;MATEID=SV_T1_1h;PURPLE_AF=0.31,0.27;PURPLE_CN=3.1,2.4;PURPLE_CN_CHANGE=0.9,0.7;PURPLE_JCN=0.95	GT:DF:SF:VF	./.:0:1:1	./.:5:3:8
chr21	43420000	SV_T1_4	C	C.	300	PASS	SVTYPE=SGL;PURPLE_AF=0.2;PURPLE_CN=2;PURPLE_CN_CHANGE=0.4;PURPLE_JCN=0.4	GT:DF:SF:VF	./.:0:0:0	./.:0:4:4
chr21	44300100	SV_T1_1h	G	]chr21:41800000]G	1500	PASS	SVTYPE=BND;MATEID=SV_T1_1o;PURPLE_AF=0.27,0.31;PURPLE_CN=2.4,3.1;PURPLE_CN_CHANGE=0.7,0.9;PURPLE_JCN=0.95	GT:DF:SF:VF	./.:0:1:1	./.:5:3:8
chr21	45990000	SV_T1_2o	A	A[chr21:46100000[	400	PON	SVTYPE=BND;MATEID=SV_T1_2h;PURPLE_AF=0.1,0.1;PURPLE_CN=2,2;PURPLE_CN_CHANGE=0.2,0.2;PURPLE_JCN=0.2	GT:DF:SF:VF	./.:0:0:0	./.:2:2:4
chr21	46100000	SV_T1_2h	G	]chr21:45990000]G	400	PON	SVTYPE=BND;MATEID=SV_T1_2o;PURPLE_AF=0.1,0.1;PURPLE_CN=2,2;PURPLE_CN_CHANGE=0.2,0.2;PURPLE_JCN=0.2	GT:DF:SF:VF	./.:0:0:0	./.:2:2:4
"

# PURPLE somatic CNV segments (T1)
mkdir -p "${BASE}/oncoanalyser_output/TEST/T1/purple"
cat > "${BASE}/oncoanalyser_output/TEST/T1/purple/T1-T.purple.cnv.somatic.tsv" << 'EOF'
chromosome	start	end	copyNumber	bafCount	observedBAF	baf	depthWindowCount	gcContent	method	minStart	maxStart
chr21	10000000	20000000	2.1	150	0.52	0.51	200	0.42	PCNI	9990000	10010000
chr21	20000000	30000000	3.8	200	0.48	0.49	300	0.43	PCNI	19990000	20010000
chr21	30000000	41000000	1.2	80	0.55	0.54	150	0.41	PCNI	29990000	30010000
chrX	49000000	51000000	1.0	120	0.50	0.50	170	0.42	PCNI	48990000	49010000
chrY	19000000	21000000	1.0	60	0.50	0.50	90	0.40	PCNI	18990000	19010000
EOF

# PURPLE gene-level CNV (T1) — needs gene + minCopyNumber + extras (script uses cut -f1-19)
cat > "${BASE}/oncoanalyser_output/TEST/T1/purple/T1-T.purple.cnv.gene.tsv" << 'EOF'
gene	chromosome	start	end	minCopyNumber	maxCopyNumber	unused1	unused2	unused3	unused4	unused5	unused6	unused7	unused8	unused9	unused10	unused11	unused12	unused13
PRDM15	chr21	41798225	41879482	3.5	4.1	0	0	0	0	0	0	0	0	0	0	0	0	0
PFKL	chr21	44300051	44327376	1.2	1.8	0	0	0	0	0	0	0	0	0	0	0	0	0
HUNK	chr21	31873020	32044633	0.3	0.5	0	0	0	0	0	0	0	0	0	0	0	0	0
XGENE_NEUTRAL	chrX	50000000	50010000	1.0	1.0	0	0	0	0	0	0	0	0	0	0	0	0	0
XGENE_DEL	chrX	60000000	60010000	0.2	0.2	0	0	0	0	0	0	0	0	0	0	0	0	0
YGENE_NEUTRAL	chrY	20000000	20010000	1.0	1.0	0	0	0	0	0	0	0	0	0	0	0	0	0
EOF

# PURPLE purity (T1) — gender for the chrX/chrY baseline, TMB/TML for the sample attributes
cat > "${BASE}/oncoanalyser_output/TEST/T1/purple/T1-T.purple.purity.tsv" << 'EOF'
purity	ploidy	gender	status	msIndelsPerMb	msStatus	tml	tmlStatus	tmbPerMb	tmbStatus
0.75	2.0	MALE	PASS	0.12	MSS	312	HIGH	12.5	HIGH
EOF

# CHORD prediction (T1 only, so T2 exercises the HRD = NA path)
mkdir -p "${BASE}/oncoanalyser_output/TEST/T1/chord"
# printf, not a heredoc: the empty remarks columns are trailing tabs, which pre-commit strips.
printf 'sample\tBRCA1\tBRCA2\thrd\thr_status\thrd_type\tremarks_hr_status\tremarks_hrd_type\nT1-T\t0.02\t0.81\t0.83\tHR_deficient\tBRCA2_type\t\t\n' \
    > "${BASE}/oncoanalyser_output/TEST/T1/chord/T1-T.chord.prediction.tsv"

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
# rows whose Name has an empty half (unannotated breakend) that must be filtered.
# PFKL--COL6A1's 5' breakend has no exonic transcript ref -> Site1_Region = Intron
cat > "${BASE}/oncoanalyser_output/TEST/T1/isofox/T1-T.isf.pass_fusions.tsv" << 'EOF'
Name	KnownType	ChromosomeUp	ChromosomeDown	PositionUp	PositionDown	OrientationUp	OrientationDown	JunctionTypeUp	JunctionTypeDown	TranscriptUp	TranscriptDown	ExonUp	ExonDown	SvType	SplitFrags	RealignedFrags	DiscordantFrags	DepthUp	DepthDown	MaxAnchorLengthUp	MaxAnchorLengthDown	CohortFrequency
PRDM15_PFKL	KNOWN_PAIR	chr21	chr21	41800000	44300000	1	-1	KNOWN	KNOWN	ENST00000398642	ENST00000349048	2	3	DEL	8	0	3	24	26	75	80	0
TMPRSS2_ERG	NONE	chr21	chr21	41500000	38400000	-1	1	CANONICAL	KNOWN	ENST00000332149	ENST00000442448	1	4	DEL	4	1	2	18	15	90	85	0.01
PFKL_COL6A1	NONE	chr21	chr21	44310000	45400000	1	-1	UNKNOWN	KNOWN		ENST00000361866	0	3	DEL	3	0	1	14	16	65	70	0
_ORPHAN3P	NONE	chr21	chr21	42000000	42500000	1	-1	CANONICAL	KNOWN		ENST00000399329	0	2	DEL	4	0	0	12	11	70	100	0
ORPHAN5P_	NONE	chr21	chr21	43000000	43500000	-1	1	KNOWN	CANONICAL	ENST00000399330		3	0	DEL	2	0	1	15	9	110	40	0.05
_	NONE	chr21	chr21	43800000	43900000	1	-1	CANONICAL	CANONICAL			0	0	DEL	3	0	0	10	10	120	120	0
EOF

echo "--- Creating T2 fixtures ---"

# PURPLE somatic VCF (T2) — DBS pairs at 200/201 and 400/401
make_vcf "${BASE}/oncoanalyser_output/TEST/T2/purple/T2-T.purple.somatic.vcf.gz" \
"$(purple_somatic_header T2)
chr21	20000100	.	G	A	55	PASS	DP=90;PURPLE_AF=0.40;PURPLE_CN=2;PURPLE_VCN=0.80;PURPLE_MACN=1;BIALLELIC;SUBCL=0	GT:AD:DP	0/0:45,0:45	0/1:35,28:63
chr21	20000200	.	C	T	60	PASS	DP=100;PURPLE_AF=0.41;PURPLE_CN=2;PURPLE_VCN=0.81;PURPLE_MACN=1;SUBCL=0	GT:AD:DP	0/0:50,0:50	0/1:40,30:70
chr21	20000201	.	C	A	58	PASS	DP=95;PURPLE_AF=0.42;PURPLE_CN=2;PURPLE_VCN=0.82;PURPLE_MACN=1;SUBCL=0	GT:AD:DP	0/0:48,0:48	0/1:38,28:66
chr21	20000400	.	T	G	52	PASS	DP=85;PURPLE_AF=0.43;PURPLE_CN=2;PURPLE_VCN=0.83;PURPLE_MACN=1;SUBCL=0	GT:AD:DP	0/0:42,0:42	0/1:34,24:58
chr21	20000401	.	A	C	50	PASS	DP=80;PURPLE_AF=0.44;PURPLE_CN=2;PURPLE_VCN=0.84;PURPLE_MACN=1;SUBCL=0	GT:AD:DP	0/0:40,0:40	0/1:32,22:54
chr21	20000600	.	A	T	48	PASS	DP=78;PURPLE_AF=0.45;PURPLE_CN=2;PURPLE_VCN=0.85;PURPLE_MACN=1;SUBCL=0	GT:AD:DP	0/0:39,0:39	0/1:30,20:50
chr21	20000800	.	C	G	62	PASS	DP=105;PURPLE_AF=0.46;PURPLE_CN=2;PURPLE_VCN=0.86;PURPLE_MACN=1;SUBCL=0	GT:AD:DP	0/0:52,0:52	0/1:42,30:72
"

# PURPLE germline VCF (T2)
make_vcf "${BASE}/oncoanalyser_output/TEST/T2/purple/T2-T.purple.germline.vcf.gz" \
"$(purple_germline_header T2)
chr21	44300200	.	T	C	75	PASS	DP=110	GT:AD:AF:DP	./.:31,31:0.5:62	./.:24,24:0.5:48
"

# sage append VCF (T2)
make_vcf "${BASE}/oncoanalyser_output/TEST/T2/sage_append/T2-T/T2-T.sage.append.vcf.gz" \
"${snv_header_t2}chr21	41800500	.	G	A	55	PASS	DP=70	GT:AD:DP	0/0:35,0:35	0/1:28,22:50
"

# PURPLE SV VCF (T2) — one PASS PRDM15-PFKL pair
make_vcf "${BASE}/oncoanalyser_output/TEST/T2/purple/T2-T.purple.sv.vcf.gz" \
"$(purple_sv_header T2)
chr21	41800500	SV_T2_1o	C	C[chr21:44300200[	1200	PASS	SVTYPE=BND;MATEID=SV_T2_1h;PURPLE_AF=0.42,0.38;PURPLE_CN=2.8,2.2;PURPLE_CN_CHANGE=1.1,0.9;PURPLE_JCN=1.05	GT:DF:SF:VF	./.:0:0:0	./.:6:4:10
chr21	44300200	SV_T2_1h	A	]chr21:41800500]A	1200	PASS	SVTYPE=BND;MATEID=SV_T2_1o;PURPLE_AF=0.38,0.42;PURPLE_CN=2.2,2.8;PURPLE_CN_CHANGE=0.9,1.1;PURPLE_JCN=1.05	GT:DF:SF:VF	./.:0:0:0	./.:6:4:10
"

# PURPLE SV VCF (T3) — two PRDM15-DNMT3L PASS pairs at the same coordinates (VF 8 and 21) for the dedupe test
make_vcf "${BASE}/oncoanalyser_output/TEST/T3/purple/T3-T.purple.sv.vcf.gz" \
"$(purple_sv_header T3)
chr21	41850000	SV_T3_1a	A	A[chr21:44250000[	900	PASS	SVTYPE=BND;MATEID=SV_T3_2a;PURPLE_AF=0.2,0.2;PURPLE_CN=2,2;PURPLE_CN_CHANGE=0.4,0.4;PURPLE_JCN=0.4	GT:DF:SF:VF	./.:0:0:0	./.:5:3:8
chr21	41850000	SV_T3_1b	A	A[chr21:44250000[	1900	PASS	SVTYPE=BND;MATEID=SV_T3_2b;PURPLE_AF=0.5,0.5;PURPLE_CN=2,2;PURPLE_CN_CHANGE=1,1;PURPLE_JCN=1	GT:DF:SF:VF	./.:0:0:0	./.:12:9:21
chr21	44250000	SV_T3_2a	G	]chr21:41850000]G	900	PASS	SVTYPE=BND;MATEID=SV_T3_1a;PURPLE_AF=0.2,0.2;PURPLE_CN=2,2;PURPLE_CN_CHANGE=0.4,0.4;PURPLE_JCN=0.4	GT:DF:SF:VF	./.:0:0:0	./.:5:3:8
chr21	44250000	SV_T3_2b	G	]chr21:41850000]G	1900	PASS	SVTYPE=BND;MATEID=SV_T3_1b;PURPLE_AF=0.5,0.5;PURPLE_CN=2,2;PURPLE_CN_CHANGE=1,1;PURPLE_JCN=1	GT:DF:SF:VF	./.:0:0:0	./.:12:9:21
"

# PURPLE somatic CNV segments (T2)
mkdir -p "${BASE}/oncoanalyser_output/TEST/T2/purple"
cat > "${BASE}/oncoanalyser_output/TEST/T2/purple/T2-T.purple.cnv.somatic.tsv" << 'EOF'
chromosome	start	end	copyNumber	bafCount	observedBAF	baf	depthWindowCount	gcContent	method	minStart	maxStart
chr21	10000000	22000000	1.8	130	0.53	0.52	180	0.42	PCNI	9990000	10010000
chr21	22000000	35000000	4.2	220	0.47	0.48	320	0.43	PCNI	21990000	22010000
chr21	35000000	42000000	2.0	100	0.50	0.50	160	0.41	PCNI	34990000	35010000
chrX	49000000	51000000	2.0	130	0.51	0.50	180	0.42	PCNI	48990000	49010000
chrY	19000000	21000000	0.0	0	0.50	0.50	0	0.40	PCNI	18990000	19010000
EOF

cat > "${BASE}/oncoanalyser_output/TEST/T2/purple/T2-T.purple.cnv.gene.tsv" << 'EOF'
gene	chromosome	start	end	minCopyNumber	maxCopyNumber	unused1	unused2	unused3	unused4	unused5	unused6	unused7	unused8	unused9	unused10	unused11	unused12	unused13
PRDM15	chr21	41798225	41879482	1.8	2.2	0	0	0	0	0	0	0	0	0	0	0	0	0
PFKL	chr21	44300051	44327376	4.5	5.0	0	0	0	0	0	0	0	0	0	0	0	0	0
HUNK	chr21	31873020	32044633	0.4	0.6	0	0	0	0	0	0	0	0	0	0	0	0	0
XGENE_FEMALE	chrX	50000000	50010000	2.0	2.0	0	0	0	0	0	0	0	0	0	0	0	0	0
YGENE_DROP	chrY	20000000	20010000	0.0	0.0	0	0	0	0	0	0	0	0	0	0	0	0	0
EOF

# PURPLE purity (T2) — no CHORD file for T2
cat > "${BASE}/oncoanalyser_output/TEST/T2/purple/T2-T.purple.purity.tsv" << 'EOF'
purity	ploidy	gender	status	msIndelsPerMb	msStatus	tml	tmlStatus	tmbPerMb	tmbStatus
0.68	2.0	FEMALE	PASS	0.08	MSS	45	LOW	2.1	LOW
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
PRDM15	100001	T1-1DT	1
PFKL	100002	T1-1DT	-1
HUNK	100003	T1-1DT	-2
PRDM15	100001	T2-T	0
PFKL	100002	T2-T	2
HUNK	100003	T2-T	-2
EOF

# SV files (output of PURPLE_SV_TO_CBIOPORTAL + Isofox, merged)
cat > "${BASE}/precomputed/T1.data_sv.txt" << 'EOF'
Sample_Id	SV_Status	Site1_Hugo_Symbol	Site1_Ensembl_Transcript_Id	Site1_Region_Number	Site1_Region	Site2_Hugo_Symbol	Site2_Ensembl_Transcript_Id	Site2_Region_Number	Site2_Region	Site2_Effect_On_Frame	NCBI_Build	Class	DNA_Support	RNA_Support	Tumor_Variant_Count	Connection_Type	Breakpoint_Type	Event_Info	Annotation	Site1_Chromosome	Site1_Position	Site2_Chromosome	Site2_Position	Tumor_Split_Read_Count	Tumor_Paired_End_Read_Count	HMF.PURPLE_JCN	HMF.PURPLE_AF_SITE1	HMF.PURPLE_AF_SITE2	HMF.PURPLE_CN_SITE1	HMF.PURPLE_CN_SITE2	HMF.PURPLE_CN_CHANGE_SITE1	HMF.PURPLE_CN_CHANGE_SITE2
T1	SOMATIC	PRDM15		.	.	PFKL		.	.		GRCh38	INVERSION	Yes	No	3	5to3	PRECISE	DNA SV: PRDM15--PFKL		21	41800000	21	44300000	3	5	0.95	0.31	0.27	3.1	2.4	0.9	0.7
EOF

cat > "${BASE}/precomputed/T2.data_sv.txt" << 'EOF'
Sample_Id	SV_Status	Site1_Hugo_Symbol	Site1_Ensembl_Transcript_Id	Site1_Region_Number	Site1_Region	Site2_Hugo_Symbol	Site2_Ensembl_Transcript_Id	Site2_Region_Number	Site2_Region	Site2_Effect_On_Frame	NCBI_Build	Class	DNA_Support	RNA_Support	Tumor_Variant_Count	Connection_Type	Breakpoint_Type	Event_Info	Annotation	Site1_Chromosome	Site1_Position	Site2_Chromosome	Site2_Position	Tumor_Split_Read_Count	Tumor_Paired_End_Read_Count	HMF.PURPLE_JCN	HMF.PURPLE_AF_SITE1	HMF.PURPLE_AF_SITE2	HMF.PURPLE_CN_SITE1	HMF.PURPLE_CN_SITE2	HMF.PURPLE_CN_CHANGE_SITE1	HMF.PURPLE_CN_CHANGE_SITE2
T2	SOMATIC	PRDM15		.	.	PFKL		.	.		GRCh38	INVERSION	Yes	No	4	5to3	PRECISE	DNA SV: PRDM15--PFKL		21	41800500	21	44300200	4	6	0.95	0.31	0.27	3.1	2.4	0.9	0.7
EOF

# Group-merged SV (used by FORMAT_PROCESS_ML_SV — named data_sv.txt)
cat > "${BASE}/precomputed/data_sv.txt" << 'EOF'
Sample_Id	SV_Status	Site1_Hugo_Symbol	Site1_Ensembl_Transcript_Id	Site1_Region_Number	Site1_Region	Site2_Hugo_Symbol	Site2_Ensembl_Transcript_Id	Site2_Region_Number	Site2_Region	Site2_Effect_On_Frame	NCBI_Build	Class	DNA_Support	RNA_Support	Tumor_Variant_Count	Connection_Type	Breakpoint_Type	Event_Info	Annotation	Site1_Chromosome	Site1_Position	Site2_Chromosome	Site2_Position	Tumor_Split_Read_Count	Tumor_Paired_End_Read_Count	HMF.PURPLE_JCN	HMF.PURPLE_AF_SITE1	HMF.PURPLE_AF_SITE2	HMF.PURPLE_CN_SITE1	HMF.PURPLE_CN_SITE2	HMF.PURPLE_CN_CHANGE_SITE1	HMF.PURPLE_CN_CHANGE_SITE2
T1-1DT	SOMATIC	PRDM15	NA	NA	NA	PFKL	NA	NA	NA	NA	GRCh38	DELETION	Yes	No	3	5to3	PRECISE	DNA WGS Fusion: PRDM15--PFKL	NA	21	41800000	21	44300000	3	5	0.95	0.31	0.27	3.1	2.4	0.9	0.7
T2-T	SOMATIC	PRDM15	NA	NA	NA	PFKL	NA	NA	NA	NA	GRCh38	DELETION	Yes	No	4	5to3	PRECISE	DNA WGS Fusion: PRDM15--PFKL	NA	21	41800500	21	44300200	4	6	0.95	0.31	0.27	3.1	2.4	0.9	0.7
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
Hugo_Symbol	Entrez_Gene_Id	T1-1DT	T2-T
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

# Group-merged mutations (used by FORMAT_ML_MUTATION — named data_mutations_dna_rna_germline.txt).
# Built from the hand-maintained pipeline-ID MAFs (T1-1DT, T2-T); skip:2 drops #version + header of the second.
{
    cat "${BASE}/precomputed/T1-1DT.somatic_rna_germline.maf"
    tail -n +3 "${BASE}/precomputed/T2-T.somatic_rna_germline.maf"
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
