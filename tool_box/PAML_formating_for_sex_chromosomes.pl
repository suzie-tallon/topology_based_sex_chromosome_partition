#!/usr/bin/perl -w

# Author : Aline Muyle

# Date: February 2024

##################################################
# 		Documentation
##################################################
# this code formats a pair of sequences for computation of the pairwise dS with PAML
# the sequences length are set to a multiple of 3
# stop codons are replaced by NNN and counted

# command : ./PAML_formating_for_sex_chromosomes.pl WORKING_FOLDER file.fasta >> PAML_formatting.txt
# example : ./PAML_formating_for_sex_chromosomes.pl /home/muyle/Documents/Etudiants_Stagiaires/Suzie_Tallon/dS_PAML 1.fasta >> PAML_formatting.txt

# input : a fasta file with X and Y sequences as such:

# >geneName_X
# ATGTGCGT...
# >geneName_Y
# ATGTTCGT...

# it is important that the X appears first and then the Y.

# output : a line that can be appended to a text file with the following columns:
# fasta_file_Name\tX_sequence_name\tX_number_internal_stop_codons\tX_positions_internal_stop_codons\tX_initial_sequence_length\tY_sequence_name\tY_number_stop_codons\tY_positions_stop_codons\tY_initial_sequence_length\tfinal_XY_sequence_length


##################################################
# 		  Modules
##################################################
# Perl modules
use strict;
use warnings;
use diagnostics;
use Tie::File;


##################################################
#		 Main code
##################################################
# opening of fasta file and sequences retrieval
my $folder_path = $ARGV[0] or die "Syntax error : ./PAML_formating_for_sex_chromosomes.pl WORKING_FOLDER file.fasta \n"; 
my $nom_fichier_entree = $ARGV[1] or die "Syntax error : ./PAML_formating_for_sex_chromosomes.pl WORKING_FOLDER file.fasta "; 
open(FICHIER_ENTREE, "<$folder_path/$nom_fichier_entree") or die("Error upon opening file $nom_fichier_entree ");
my @alignement = <FICHIER_ENTREE> ;
chomp(@alignement) ;
my $nb_ligne = scalar(@alignement);

# opening output file
my $gene = $nom_fichier_entree ;
$gene =~ s/\.fst// ;
$gene =~ s/\.fasta// ;
open(OUTPUT, ">$folder_path/PAML_format/$gene.fasta") ;

# empty lists of sequences, sequence names, number of stop codons and their positions,
# initial and final sequence length
my @sequences = () ;
my @names = () ;
my @number_stop_codons = () ;
my @positions_stop_codons = () ;
my @initial_sequence_length = () ;
my @final_sequence_length = () ;

# initialize output line with gene name
my $output_line = "$gene" ;


# loop to retrive sequences in case they are written on multiple lines
for (my $i=0 ; $i<$nb_ligne; $i+=1) {
	my $debut = substr($alignement[$i], 0, 1) ;
	if ($debut eq "\>") {
		my $nom_seq = $alignement[$i] ;
		$nom_seq =~s/\>// ;
		$nom_seq =~s/\-/\_/g ;
		chomp($nom_seq) ;
		push(@names, $nom_seq) ;
		my $seq = "" ;
		my $j = $i + 1 ;
		while ($j < $nb_ligne) {
			my $debut = substr($alignement[$j], 0, 1) ;
			last if $debut eq ">" ;
			my $data = $alignement[$j] ;
			chomp($data) ;
			$seq = $seq.$data ;
			$j++ ;
		}
		push(@sequences, $seq) ;
	}
}

# Stop codons are replaced by NNN and counted for each sequence
# if the sequence has length which is not a multiple of 3 the hanging bases are cut at the end
# initializing the minimum sequence length to the length of the first sequence:
my $minimun_length = length($sequences[0]) ;
for (my $i=0 ; $i<scalar(@names); $i+=1) {
	my $seq = $sequences[$i] ;
	
	# save initial sequence length
	$initial_sequence_length[$i] = length($seq);
	
	# cutting end of sequence if not a multiple of 3
	my $reste_euclide = $initial_sequence_length[$i] % 3 ;
	$seq = substr($seq, 0, $initial_sequence_length[$i]-$reste_euclide) ;
	$final_sequence_length[$i] = length($seq);
	
	# updating minimum length if needed, to keep track of the gene smallest sequence
	if ($minimun_length > $final_sequence_length[$i]) {
		$minimun_length = $final_sequence_length[$i] ;
	}

	# replacing internal stop codons by NNN, count them and register their position
	$number_stop_codons[$i] = 0 ;
	$positions_stop_codons[$i] = '' ;
	for (my $b=0 ; $b<($final_sequence_length[$i]-3); $b+=3) {
		my $current_codon = substr($seq, $b, 3) ;
		if (($current_codon eq "TAA")||($current_codon eq "TAG")||($current_codon eq "TGA")) {
			$seq = substr($seq, 0, $b) . 'NNN' . substr($seq, $b + 3, $final_sequence_length[$i]) ;
			$number_stop_codons[$i]++ ;
			my $current_position = $b+1 ;
			$positions_stop_codons[$i] .= "," . $current_position ;
		}
	}
	
	# replacing final stop codon by NNN if needed, not counting it because it's not internal
	my $final_exon_position = $final_sequence_length[$i] - 3 ;
	my $current_codon = substr($seq, $final_exon_position, 3) ;
	if (($current_codon eq "TAA")||($current_codon eq "TAG")||($current_codon eq "TGA")) {
		$seq = substr($seq, 0, $final_exon_position) . 'NNN' ;
	}

	# save sequence to scalar
	$sequences[$i] = $seq ;
}

# if one sequence is longer than the other, all sequences are cut to the same length
# preparation of output line
for (my $i=0 ; $i<scalar(@names); $i+=1) {
	my $seq = $sequences[$i] ;
	$seq = substr($seq, 0, $minimun_length) ;
	$final_sequence_length[$i] = length($seq);
	# save sequence to scalar
	$sequences[$i] = $seq ;
	
	# prepare output line
	$output_line .= "\t" . $names[$i] . "\t" . $number_stop_codons[$i] . "\t" . $positions_stop_codons[$i] . "\t" . $initial_sequence_length[$i] ;
}

$output_line .= "\t" . $minimun_length . "\n" ;
print($output_line) ;


# ecriture des sequences dans le fichier 
for (my $k=0 ; $k<scalar(@names); $k+=1) {
	my $seq = $sequences[$k] ;
	my $seq_name = $names[$k] ;
	print(OUTPUT "\>$seq_name\n$seq\n") ;
}

close(OUTPUT) ;


