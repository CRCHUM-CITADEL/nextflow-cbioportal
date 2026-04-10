#!/usr/bin/env bash
# Test incremental pipeline runs: verifies that already-processed samples are skipped
# when re-running with an updated samplesheet containing new samples.
#
# Usage: bash bin/test_incremental.sh
# Requires: -profile test (Apptainer + python-ampliseq.sif container)
set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUTDIR="$PIPELINE_DIR/test_incremental_results"

echo "Cleaning up previous test run..."
rm -rf "$OUTDIR" "$PIPELINE_DIR/.nextflow" "$PIPELINE_DIR/work"

echo ""
echo "=== Step 1: Run with SAMPLE_001 only ==="
nextflow run "$PIPELINE_DIR/main.nf" \
  -profile test \
  --input "$PIPELINE_DIR/assets/samplesheet_sample1.csv" \
  --outdir "$OUTDIR"

echo ""
echo "--- Checking SAMPLE_001 per-sample outputs ---"
for suffix in _sv.txt _cna.txt _seg.txt _mutations.txt; do
  FILE="$OUTDIR/samples/SAMPLE_001/SAMPLE_001$suffix"
  if [ -f "$FILE" ]; then
    echo "PASS: $FILE exists"
  else
    echo "FAIL: $FILE not found"
    exit 1
  fi
done

echo ""
echo "--- Checking merged outputs contain SAMPLE_001 ---"
for data_file in data_sv.txt data_cna.txt data_mutations.txt data_seg.txt; do
  FILE="$OUTDIR/$data_file"
  if [ -f "$FILE" ]; then
    echo "PASS: $FILE exists"
  else
    echo "FAIL: $FILE not found"
    exit 1
  fi
done

# Record mtime of SAMPLE_001 per-sample outputs before incremental run
MTIME_SV_BEFORE=$(stat -c %Y "$OUTDIR/samples/SAMPLE_001/SAMPLE_001_sv.txt")
MTIME_MUT_BEFORE=$(stat -c %Y "$OUTDIR/samples/SAMPLE_001/SAMPLE_001_mutations.txt")

echo ""
echo "=== Step 2: Run incrementally with SAMPLE_001 + SAMPLE_002 ==="
nextflow run "$PIPELINE_DIR/main.nf" \
  -profile test \
  --input "$PIPELINE_DIR/assets/samplesheet.csv" \
  --outdir "$OUTDIR"

echo ""
echo "--- Checking SAMPLE_002 per-sample outputs created ---"
for suffix in _sv.txt _cna.txt _seg.txt _mutations.txt; do
  FILE="$OUTDIR/samples/SAMPLE_002/SAMPLE_002$suffix"
  if [ -f "$FILE" ]; then
    echo "PASS: $FILE exists"
  else
    echo "FAIL: $FILE not found"
    exit 1
  fi
done

echo ""
echo "--- Checking SAMPLE_001 was NOT reprocessed (mtime unchanged) ---"
MTIME_SV_AFTER=$(stat -c %Y "$OUTDIR/samples/SAMPLE_001/SAMPLE_001_sv.txt")
MTIME_MUT_AFTER=$(stat -c %Y "$OUTDIR/samples/SAMPLE_001/SAMPLE_001_mutations.txt")

if [ "$MTIME_SV_BEFORE" = "$MTIME_SV_AFTER" ]; then
  echo "PASS: SAMPLE_001_sv.txt was not reprocessed (mtime unchanged)"
else
  echo "FAIL: SAMPLE_001_sv.txt was reprocessed (mtime changed)"
  exit 1
fi

if [ "$MTIME_MUT_BEFORE" = "$MTIME_MUT_AFTER" ]; then
  echo "PASS: SAMPLE_001_mutations.txt was not reprocessed (mtime unchanged)"
else
  echo "FAIL: SAMPLE_001_mutations.txt was reprocessed (mtime changed)"
  exit 1
fi

echo ""
echo "=== All incremental tests passed ==="
