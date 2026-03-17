#!/bin/bash

# Author : Suzie Tallon

# Date: February 2026

#script use to divide a gen file with many genes into files with one gene per file

# Input argument check
if [ $# -ne 2 ]; then
    echo "Usage: $0 <list_gen> <file_gen>"
    exit 1
fi

# Output directory (created if not already existing)
output_dir="gene_selection_11_03_25"
mkdir -p "$output_dir"

list_gen="$1"
file_gen="$2"

while IFS= read -r f; do
    echo "$f"
    nb_ligne=$(grep -n '>' "$file_gen" | grep -A 1 "$f" | awk -F ':' 'NR==2 {print $1 - prev  -1} {prev = $1}')
	#write output
	output_file="$output_dir/${f}.gen"
    grep -A ${nb_ligne} $f  $file_gen | tail -n +2 > "$output_file"
done < "$list_gen"


	

