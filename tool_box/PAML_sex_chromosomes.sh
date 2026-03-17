#!/bin/bash

# Author : Aline Muyle

# Date: February 2024

##################################################
# 		Documentation
##################################################
# This code prepares X/Y sequences obtained with SEX-DETector or SDpop and
# runs PAML to compute pairwise dS values
# internal stop codons are replaced by NNN and the output gives the information
# on how many stop codons were found in the X and the Y and their position


# command : PAML_sex_chromosomes.sh

# input : the working directory must contain the input X and Y sequences, 1 file per gene
# in fasta format as such:

# >geneName_X
# ATGTGCGT...
# >geneName_Y
# ATGTTCGT...

# it is important that the X appears first and then the Y.
# The genes must each have their separate fasta files (not multiple genes per file).

# output : dS.txt file with the computed dS for each gene and number of stop codons 
# for the X and Y, 1 line per gene


##################################################
#	     parameters
##################################################
# create variable containing current folder path
WORKING_FOLDER=$(pwd)


##################################################
#		 Main code
##################################################
cd $WORKING_FOLDER

# clean up previous runs
rm -r $WORKING_FOLDER/PAML_format PAML_formatting.txt dS.txt

# create folder where aligments formated for PAML will be stored after replacing stop codons
# by 'NNN'
mkdir -p $WORKING_FOLDER/PAML_format

# create output folder for PAML analysis
mkdir -p $WORKING_FOLDER/PAML_format/PAML_output

# create header (column names) of the dS output file
echo -e "fasta_file_Name\tX_sequence_name\tX_number_internal_stop_codons\tX_positions_internal_stop_codons\tX_initial_sequence_length\tY_sequence_name\tY_number_stop_codons\tY_positions_stop_codons\tY_initial_sequence_length\tfinal_XY_sequence_length\tdS" > dS.txt ;

# run PAML on each fasta file present in the working directory 
# initiate pair number to 1, for each alignment the pair number will later be incremented by 1
pair=1
for f in *.fasta ;
	do

	filename=`echo $f | sed 's/.fasta//'` ;
	echo -e "\n$filename" ;
	
	# format for PAML, keep information on alignment in file PAML_formatting.txt
	PAML_formating_for_sex_chromosomes.pl $WORKING_FOLDER $f >> PAML_formatting.txt ;
	cp ./PAML_format/$f ./PAML_format/$pair.fasta
	sed -i 's/:/_/g' ./PAML_format/$pair.fasta
	sed -i 's/(/_/g' ./PAML_format/$pair.fasta
	sed -i 's/)/_/g' ./PAML_format/$pair.fasta
	sed -i 's/!/N/g' ./PAML_format/$pair.fasta
	
	# write PAML ctl file for contig
	echo -e "seqfile = $WORKING_FOLDER/PAML_format/$pair.fasta" > ./PAML_format/$pair.ctl ;
	echo -e "outfile = $WORKING_FOLDER/PAML_format/PAML_output/$pair.out" >> ./PAML_format/$pair.ctl; 
	echo -e "verbose = 0" >> ./PAML_format/$pair.ctl ;
	echo -e "icode = 0" >> ./PAML_format/$pair.ctl ;
	echo -e "weighting = 0" >> ./PAML_format/$pair.ctl ;
	echo -e "commonf3x4 = 0" >> ./PAML_format/$pair.ctl ;
	
	# run PAML on each contig
	yn00 ./PAML_format/$pair.ctl > ./PAML_format/PAML_output/$pair.console ;
	mv 2YN.dN ./PAML_format/PAML_output/$pair.dN ;
	mv 2YN.dS ./PAML_format/PAML_output/$pair.dS ;
	mv 2YN.t ./PAML_format/PAML_output/$pair.t ;
	rm rst rst1 rub ;
	
	# extract dS and add it to columns of ./PAML_formatting.txt
	dS=`sed -n '3p' ./PAML_format/PAML_output/$pair.dS | awk '{print $2}' | sed 's/\-0.0000/0/' | sed 's/-nan/NA/'`;
	reformat=`grep -w $filename PAML_formatting.txt`;
	echo -e "$reformat\t$dS" >> dS.txt ;
	pair=`expr $pair + 1`;

done


