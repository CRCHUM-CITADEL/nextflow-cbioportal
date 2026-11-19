"""Tests for bin/build_ensembl_annotations.py."""

import csv
import gzip
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).parent.parent / "bin" / "build_ensembl_annotations.py"

GFF3 = "\n".join([
    "##gff-version 3",
    "17\tensembl_havana\tgene\t43044292\t43170245\t.\t-\t.\tID=gene:ENSG00000012048;Name=BRCA1;"
    "biotype=protein_coding;description=BRCA1 DNA repair associated [Source:HGNC Symbol%3BAcc:HGNC:1100];"
    "gene_id=ENSG00000012048;version=27",
    "17\tensembl_havana\tmRNA\t43044292\t43125364\t.\t-\t.\tID=transcript:ENST00000357654;"
    "Parent=gene:ENSG00000012048;Name=BRCA1-203;transcript_id=ENST00000357654",
    "1\thavana\tncRNA_gene\t100\t200\t.\t+\t.\tID=gene:ENSG00000000001;biotype=lncRNA;gene_id=ENSG00000000001",
    "14\thavana\tV_gene_segment\t300\t400\t.\t+\t.\tID=transcript:ENST00000542354;"
    "Parent=gene:ENSG00000255569;Name=TRAV1-1-201;biotype=TR_V_gene;transcript_id=ENST00000542354",
    "",
])

ENTREZ = "\n".join([
    "gene_stable_id\ttranscript_stable_id\tprotein_stable_id\txref\tdb_name\tinfo_type",
    "ENSG00000012048\tENST00000357654\tENSP00000350283\t672\tEntrezGene\tDEPENDENT",
    "ENSG00000012048\tENST00000471181\t-\t672\tEntrezGene\tDEPENDENT",
    "ENSG00000012048\tENST00000471181\t-\t9999\tEntrezGene\tDEPENDENT",
    "ENSG00000012048\tENST00000357654\t-\tBRCA1\tEntrezGene_trans_name\tMISC",
    "",
])


def build(tmp_path):
    gff3 = tmp_path / "genes.gff3.gz"
    entrez = tmp_path / "entrez.tsv.gz"
    out = tmp_path / "out.tsv"
    with gzip.open(gff3, "wt") as fh:
        fh.write(GFF3)
    with gzip.open(entrez, "wt") as fh:
        fh.write(ENTREZ)
    subprocess.run([sys.executable, str(SCRIPT), str(gff3), str(entrez), str(out)], check=True)
    with open(out) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_header_matches_downstream_readers(tmp_path):
    rows = build(tmp_path)
    assert list(rows[0].keys()) == ["ensembl_id", "gene_symbol", "entrez_ncbi_id", "biotype",
                                    "chr", "start", "end", "strand", "description"]


def test_one_row_per_gene_entrez_pair(tmp_path):
    brca1 = [r for r in build(tmp_path) if r["ensembl_id"] == "ENSG00000012048"]
    # Duplicate xrefs collapse; EntrezGene_trans_name rows are not Entrez ids.
    assert [r["entrez_ncbi_id"] for r in brca1] == ["672", "9999"]
    assert brca1[0]["gene_symbol"] == "BRCA1"
    assert brca1[0]["strand"] == "-1"
    assert brca1[0]["description"] == "BRCA1 DNA repair associated [Source:HGNC Symbol;Acc:HGNC:1100]"


def test_gene_without_entrez_or_name_keeps_empty_fields(tmp_path):
    lnc = [r for r in build(tmp_path) if r["ensembl_id"] == "ENSG00000000001"]
    assert len(lnc) == 1
    assert lnc[0]["entrez_ncbi_id"] == ""
    assert lnc[0]["gene_symbol"] == ""
    assert lnc[0]["strand"] == "1"


def test_rows_sorted_by_ensembl_id(tmp_path):
    # The GFF3 lists BRCA1 (ENSG00000012048) before ENSG00000000001.
    ids = [r["ensembl_id"] for r in build(tmp_path)]
    assert ids == ["ENSG00000000001", "ENSG00000012048", "ENSG00000012048"]


def test_transcript_records_are_skipped(tmp_path):
    ids = {r["ensembl_id"] for r in build(tmp_path)}
    assert ids == {"ENSG00000012048", "ENSG00000000001"}
