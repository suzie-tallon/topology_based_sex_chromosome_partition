#!/bin/bash

# Author : Suzie Tallon

# Date: February 2026

# create an input list in the right format for STAR

# Get the current folder path
current_folder=$(pwd)

# Output file name
OutputFile="input_file_list.txt"

# Remove the output file if it already exists
rm -f "$OutputFile"

# Iterate through files in the QC_trimming directory
for file in "$current_folder"/*_R1.fq.gz; do
    #Extract the file name without extension and path
    base_name=${file##*/}
    base_name=${base_name%_R[12]*}

    # Write lines to the output file
    echo "${base_name}_R1.fq.gz   ${base_name}_R2.fq.gz   ${base_name}" >> $OutputFile
done

echo "The file $OutputFile has been successfully created."
cat $OutputFile
