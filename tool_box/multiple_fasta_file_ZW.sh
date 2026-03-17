#!/bin/bash

# Author : Suzie Tallon

# Date: February 2026

#script use to divide a fasta file with many genes into files with one gene per file with the Z and W sequence

# Input argument check
if [ $# -ne 1 ]; then
    echo "Usage: $0 <fichier_donnees>"
    exit 1
fi

# Initialisation
block_count=0
current_filename=""

# Reading the input file line by line starting from the second line
tail -n +2 "$1" | while IFS= read -r line; do
    ((block_count++))

    # if new bloc extract file name 
    if ((block_count == 1)); then
        current_filename=$(echo "$line" | cut -d ' ' -f 1)
    fi

    # writing output 
    echo "$line" >> "${current_filename}.fasta"

    # re-initialisation
    if ((block_count == 4)); then
        block_count=0
    fi

done

# delete empty files
find . -type f -empty -delete

#move every thing into a new folder

mkdir fasta_to_cut
mv *_Z.fasta fasta_to_cut
cd fasta_to_cut

# rename file
for file in *_Z.fasta; do
    mv "$file" "${file/_Z/}"
done


for file in *.fasta; do
    mv "$file" "${file/>/}"
done

cd ..
