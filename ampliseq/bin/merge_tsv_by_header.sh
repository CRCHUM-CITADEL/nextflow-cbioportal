#!/usr/bin/env bash
# Merge TSV files that share a header line but may differ in columns (e.g. per-sample outputs
# written by different pipeline versions or writers). Pass 1 builds the union of columns in
# first-seen order; pass 2 writes every row mapped onto it, missing columns left empty.
# Header-only and empty (0-byte) files are allowed.
# Usage: merge_tsv_by_header.sh <output.tsv> <input.tsv>...
set -euo pipefail
out=$1; shift
[ $# -gt 0 ] || { echo "merge_tsv_by_header.sh: no input files" >&2; exit 1; }

awk -F'\t' -v OFS='\t' '
    pass == 1 {
        if (FNR == 1) for (i = 1; i <= NF; i++) if (!($i in seen)) { seen[$i] = 1; cols[++n] = $i }
        next
    }
    FNR == 1 {
        if (!printed) { line = cols[1]; for (j = 2; j <= n; j++) line = line OFS cols[j]; print line; printed = 1 }
        split("", idx); for (i = 1; i <= NF; i++) idx[$i] = i
        next
    }
    {
        line = ""
        for (j = 1; j <= n; j++) line = line (j > 1 ? OFS : "") ((cols[j] in idx) ? $(idx[cols[j]]) : "")
        print line
    }
    END { if (!printed && n) { line = cols[1]; for (j = 2; j <= n; j++) line = line OFS cols[j]; print line } }
' pass=1 "$@" pass=2 "$@" > "$out"
