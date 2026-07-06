#!/usr/bin/env python3
"""
Convert a fusion VCF (*-star-fusion.final.vcf) to cBioPortal SV format.
Usage: fusion_vcf_to_sv.py <vcf_file> <sample_id>

Output: data_sv.txt (tab-separated, appended)

Only PASS records are included. Breakend pairs are deduplicated by
keeping only the first record of each pair (ID ending with _1).
"""
import re
import sys
import os


def parse_info(info_str):
    """Parse VCF INFO field into a dict."""
    fields = {}
    for field in info_str.split(';'):
        if '=' in field:
            key, val = field.split('=', 1)
            fields[key] = val
        else:
            fields[field] = True
    return fields


def parse_alt_position(alt):
    """Extract chrom and pos from breakend ALT field.

    Formats: .]chr8:128751238]  or  ]chr12:53586114].  or  [chr...[  etc.
    Returns (chrom, pos) or (None, None).
    """
    m = re.search(r'(chr[\dXYM]+):(\d+)', alt)
    if m:
        return m.group(1), int(m.group(2))
    return None, None


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <vcf_file> <sample_id>", file=sys.stderr)
        sys.exit(1)

    vcf_file = sys.argv[1]
    sample_id = sys.argv[2]
    out_file = "data_sv.txt"

    rows = []
    with open(vcf_file) as fh:
        for line in fh:
            if line.startswith('#'):
                continue

            cols = line.rstrip('\n').split('\t')
            if len(cols) < 8:
                continue

            chrom = cols[0]
            pos = int(cols[1])
            record_id = cols[2]
            alt = cols[4]
            filter_val = cols[6]
            info_str = cols[7]

            if filter_val != 'PASS':
                continue

            # Deduplicate breakend pairs: keep only _1 records
            if not record_id.endswith('_1'):
                continue

            info = parse_info(info_str)

            gene_name = info.get('GENE_NAME', '')
            driver_gene = info.get('FUSION_DRIVER_GENE', '')
            split_reads = info.get('SPLIT_READS', '0')

            # Strip chr prefix
            site1_chrom = chrom[3:] if chrom.startswith('chr') else chrom

            # Parse Site2 position from ALT for SV_Length
            _, site2_pos = parse_alt_position(alt)
            sv_length = abs(pos - site2_pos) if site2_pos is not None else 0

            rows.append((
                sample_id,
                "Somatic",
                gene_name,
                site1_chrom,
                pos,
                driver_gene,
                "FUSION",
                f"RNA-Seq FUSION : {gene_name}-{driver_gene}",
                split_reads,
                sv_length,
            ))

    write_header = not os.path.exists(out_file)
    with open(out_file, 'a') as out:
        if write_header:
            out.write('\t'.join([
                'Sample_Id', 'SV_Status', 'Site1_Hugo_Symbol',
                'Site1_Chromosome', 'Site1_Region', 'Site2_Hugo_Symbol',
                'Class', 'Event_Info', 'Tumor_Variant_Count', 'SV_Length',
            ]) + '\n')
        for row in rows:
            out.write('\t'.join(str(x) for x in row) + '\n')

    print(f"{'Written' if write_header else 'Appended'}: {out_file} ({len(rows)} fusion(s))")


if __name__ == '__main__':
    main()
