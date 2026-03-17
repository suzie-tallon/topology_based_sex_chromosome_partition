#!/bin/bash

# Author : Suzie Tallon

# Date: February 2026

#script to create the sample list 

# Get the current folder path
current_folder=$(pwd)

# Output file name
OutputFile="samples_list.txt"

# Remove the output file if it already exists
rm -f "$OutputFile"

# Iterate through files in the QC_trimming directory
for file in "$current_folder"/*_R1.fq.gz; do
    #Extract the file name without extension and path
    base_name=${file##*/}
    base_name=${base_name%_R[12]*}

    # Write output file
    echo "${base_name}" >> $OutputFile
done

echo "The file $OutputFile has been successfully created."
cat $OutputFile
