#!/usr/bin/env python3
"""
Convert a fusion VCF (*-star-fusion.final.vcf, Illumina RNA fusion caller) to cBioPortal SV format
(full column set from cbio_sv.SV_COLUMNS).
Usage: fusion_vcf_to_sv.py <vcf_file> <sample_id>

Output: data_sv.txt (tab-separated, appended)

Only PASS records are used. Each fusion is a breakend pair with IDs <name>_1 (5' partner,
GENE_NAME) and <name>_2 (3' partner, FUSION_DRIVER_GENE); the pair becomes one SV row, with each
site's position and EXON_NUM taken from its own record. An unpaired record falls back to the
mate position in its BND ALT.
"""
import re
import sys

from cbio_sv import NCBI_BUILD, sv_length, strip_chr, write_rows

# BND ALT mate position: N[chr:pos[, N]chr:pos], [chr:pos[N, ]chr:pos]N (contig with or without chr)
ALT_MATE = re.compile(r"[\[\]]([^\[\]:]+):(\d+)[\[\]]")


def parse_info(info_str):
    """Parse VCF INFO field into a dict."""
    fields = {}
    for field in info_str.split(';'):
        if '=' in field:
            key, val = field.split('=', 1)
            fields[key] = val
        elif field:
            fields[field] = True
    return fields


def parse_alt_position(alt):
    """Mate (chrom, pos) from a breakend ALT, or (None, None)."""
    m = ALT_MATE.search(alt)
    return (m.group(1), int(m.group(2))) if m else (None, None)


def to_int(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <vcf_file> <sample_id>", file=sys.stderr)
        sys.exit(1)

    vcf_file, sample_id = sys.argv[1], sys.argv[2]

    # name -> {"1": record, "2": record}; records keep input order via `order`
    pairs, order = {}, []
    with open(vcf_file) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 8 or cols[6] != 'PASS':
                continue
            rec = {"chrom": cols[0], "pos": int(cols[1]), "id": cols[2], "alt": cols[4],
                   "info": parse_info(cols[7])}
            m = re.match(r"^(.*)_([12])$", rec["id"])
            name, end = (m.group(1), m.group(2)) if m else (f"{rec['id']}@{len(order)}", "1")
            if name not in pairs:
                pairs[name] = {}
                order.append(name)
            pairs[name].setdefault(end, rec)

    rows = []
    for name in order:
        site1, site2 = pairs[name].get("1"), pairs[name].get("2")
        lead = site1 or site2
        info = lead["info"]
        if site1 is None:            # only the 3' record passed: rebuild site1 from its ALT
            chrom1, pos1 = parse_alt_position(site2["alt"])
            chrom2, pos2 = site2["chrom"], site2["pos"]
        else:
            chrom1, pos1 = site1["chrom"], site1["pos"]
            chrom2, pos2 = (site2["chrom"], site2["pos"]) if site2 else parse_alt_position(site1["alt"])

        gene1 = info.get('GENE_NAME', '')
        gene2 = info.get('FUSION_DRIVER_GENE', '')
        split_reads = to_int(info.get('SPLIT_READS'))
        discordant = to_int(info.get('DISCORDANT_PAIRS'))
        support = sum(v for v in (split_reads, discordant) if v is not None) \
            if split_reads is not None or discordant is not None else None
        annotation = info.get('ANNOTATION', '')
        if annotation in ('.', True):
            annotation = ''

        rows.append({
            "Sample_Id":                   sample_id,
            "SV_Status":                   "SOMATIC",
            "Site1_Hugo_Symbol":           gene1,
            "Site1_Region_Number":         site1["info"].get("EXON_NUM", "") if site1 else "",
            "Site1_Chromosome":            strip_chr(chrom1) if chrom1 else "",
            "Site1_Position":              pos1,
            "Site2_Hugo_Symbol":           gene2,
            "Site2_Region_Number":         site2["info"].get("EXON_NUM", "") if site2 else "",
            "Site2_Chromosome":            strip_chr(chrom2) if chrom2 else "",
            "Site2_Position":              pos2,
            "NCBI_Build":                  NCBI_BUILD,
            "Class":                       "FUSION",
            "Length":                      sv_length(chrom1, pos1, chrom2, pos2),
            "DNA_Support":                 "No",
            "RNA_Support":                 "Yes",
            "Tumor_Variant_Count":         support,
            "Tumor_Split_Read_Count":      split_reads,
            "Tumor_Paired_End_Read_Count": discordant,
            "Connection_Type":             "5to3",
            "Breakpoint_Type":             "PRECISE",
            "Event_Info":                  f"RNA-seq Fusion: {gene1}--{gene2}",
            "Annotation":                  f"{gene1} - {gene2} fusion" if gene1 and gene2 else "",
            "External_Annotation":         annotation,
            "Comments":                    "Non-targeted fusion" if info.get("Non-Targeted") is True else "",
        })

    write_rows(rows)


if __name__ == '__main__':
    main()
