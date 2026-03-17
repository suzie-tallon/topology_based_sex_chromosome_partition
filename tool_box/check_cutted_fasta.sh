#!/bin/bash

# Author : Suzie Tallon

# Date: February 2026

#script to check if a sequence is valid : length divisible by three and valid start and stop codon 

# check arguments
if [ $# -ne 1 ]; then
    echo "Usage: $0 <nom_dossier>"
    exit 1
fi

folder="$1"

# Define valid start and stop codon
start_patterns=("ATG" "NTG" "ANG" "ATN" "ANN" "NTN" "NNG" "NNN")
end_patterns=("TGA" "NGA" "TNA" "TGN" "TNN" "NGN" "NNA" "TAG" "NAG" "TNG" "TAN" "NNG" "TAA" "NAA" "NAN" "NNN")

# output file
output_file="cut_check"
> "$output_file"  # remove it if already existing

# on all fasta file in the folder
for file in "$folder"/*.fasta; do
    
    if [ -f "$file" ]; then
        echo "treating file $file :"
        
        # initialisation
        div_by_3_ok=true
        start_ok=true
        stop_ok=true
        
        while IFS= read -r line; do
            # ignore the line if it's the sequence name 
            if [[ "$line" == ">"* ]]; then
                continue
            fi
            
            # The length must be divisible by three
            seq_length=$(echo -n "$line" | wc -c)
            if (( seq_length % 3 != 0 )); then
                div_by_3_ok=false
            fi
            
            # valid start codon ?
            valid_start=false
            for pattern in "${start_patterns[@]}"; do
                if [[ "$line" == $pattern* ]]; then
                    valid_start=true
                    break
                fi
            done
            
            if ! $valid_start; then
                start_ok="${line:0:3} FALSE"  # print the three first letters and "FALSE"
            fi
            
            # valid stop codon ?
            valid_end=false
            for pattern in "${end_patterns[@]}"; do
                if [[ "$line" == *$pattern ]]; then
                    valid_end=true
                    break
                fi
            done
            
            if ! $valid_end; then
                stop_ok="${line: -3} FALSE"  # print the three last letters and "FALSE"
            fi
        done < "$file"
        
        # write output file
        echo "$(basename "$file") div 3 ${div_by_3_ok^^} start ${start_ok^^} stop ${stop_ok^^}" >> "$output_file"
        
        echo "checking done for $file."
        echo
    fi
done


