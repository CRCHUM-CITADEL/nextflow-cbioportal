#!/usr/bin/env bash

# Build a Nextflow samplesheet (sample,quant_file) from oarfish *.quant files.

#

# Usage: ./make_samplesheet.sh <oarfish_output_dir> [output.csv]

# Example:

#   ./make_samplesheet.sh /lustre06/.../quantification/oarfish/output samplesheet.csv



set -euo pipefail



if [[ $# -lt 1 ]]; then

    echo "Usage: $0 <oarfish_output_dir> [output.csv]" >&2

    exit 1

fi



input_dir=$(realpath "$1")

output_csv="${2:-samplesheet.csv}"



if [[ ! -d "$input_dir" ]]; then

    echo "Error: '$input_dir' is not a directory" >&2

    exit 1

fi



# Collect .quant files (top level only), sorted by name

mapfile -t quant_files < <(find "$input_dir" -maxdepth 1 -type f -name '*.quant' | sort)



if [[ ${#quant_files[@]} -eq 0 ]]; then

    echo "Error: no .quant files found in '$input_dir'" >&2

    exit 1

fi



{

    echo "sample,quant_file"

    for f in "${quant_files[@]}"; do

        sample=$(basename "$f" .quant)

        echo "${sample},${f}"

    done

} > "$output_csv"



echo "Wrote ${#quant_files[@]} samples to ${output_csv}"
