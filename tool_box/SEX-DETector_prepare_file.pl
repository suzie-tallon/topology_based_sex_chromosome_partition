#!/usr/bin/perl

# Programmer : Aline Muyle aline.muyle@univ-lyon1.fr
# This code is provided under licence CeCILL B.
# This implies that it can be used by anyone, anyone who uses it has to cite the corresponding article and anyone who uses this code and modifies it to make a new code and publish it has to make the new code available to everyone.
# More precisely, this code cannot be used in any kind of commercial way. 

## Date and Version
my $version= '1.0';
my $version_date = '26th May 2016';

##################################################
# 		Documentation
##################################################
# see manual at https://lbbe.univ-lyon1.fr/Download-5251.html

# command : ./SEX-DETector/SEX-DETector_prepare_file.pl input.alr_gen output.alr_gen_summary -hom homogametic_progeny_name1,homogametic_progeny_name2 -het heterogametic_progeny_name1,heterogametic_progeny_name2 -hom_par homogametic_parent_name -het_par heterogametic_parent_name
# example : ./SEX-DETector/SEX-DETector_prepare_file.pl jeu_test.alr_gen jeu_test.alr_gen_summary -hom C1_26_female,C1_27_female,C1_29_female,C1_34_female -het C1_01_male,C1_3_male,C1_04_male,C1_05_male -hom_par U10_37_mother -het_par Leuk144-3_father 

# input : alr_gen file from reads2snp
# output header : >alr_gen_summary \n number_time_line_happens individual_names(what is important is that individual are placed in the right category (heterogametic or homogametic, the real genotype of a given individual doesn't matter here)
# other lines : number_time_line_happens	mother_genotype	father_genotype homogametic_progeny_genotypes_separated_by_'\t'	heterogametic_progeny_genotypes_separated_by_'\t'


##################################################
# 		  Modules
##################################################
# Perl modules
use strict;
use warnings;
use diagnostics;
use Tie::File;
use Getopt::Long ;
use File::Basename;

##################################################
#	     	   variables
##################################################
my %Parameters ;
my %Table ; # key=line, value=number of time the line appears in the file.

##################################################
#	     retrieving parameters
##################################################
$Parameters{'input_file_name'} = $ARGV[0] or die('Syntax error : ./SEX-DETector_prepare_file.pl input.alr_gen output.alr_gen_summary'."\n\n") ; 
$Parameters{'output_file_name'} = $ARGV[1] or die('Syntax error : ./SEX-DETector_prepare_file.pl input.alr_gen output.alr_gen_summary'."\n\n") ; 

my %alr_gen_columns ; # keys = individuals' names, values = column number in alr_gen file
$Parameters{'alr_gen_columns'} = \%alr_gen_columns ;
my $homogametic_input;
my $heterogametic_input;

GetOptions (
	'homogametic|hom=s'			=>	\$homogametic_input,
	'heterogametic|het=s'			=>	\$heterogametic_input,
	'homogametic_parent|hom_par=s'		=>	\$Parameters{'homogametic_parent_name'},
	'heterogametic_parent|het_par=s'	=>	\$Parameters{'heterogametic_parent_name'},
);

my @homogametic = split(/,/, $homogametic_input) ; # contains homogametic progeny individual names
$Parameters{'homogametic'} = \@homogametic ; 
my @heterogametic = split(/,/, $heterogametic_input) ; # contains heterogametic progeny individual names
$Parameters{'heterogametic'} = \@heterogametic ;

# Check parameters
checkParameters(\%Parameters) ;


##################################################
#	 opening input and output files
##################################################
open(INPUT_FILE, "<$Parameters{'input_file_name'}") or die('Error ! Cannot open input file '.$Parameters{'input_file_name'}."\n") ;
open(OUTPUT_FILE, '>'.$Parameters{'output_file_name'}) or die('Error ! Cannot open output file '.$Parameters{'output_file_name'}."\n") ;

my @INPUT_FILE ;
tie @INPUT_FILE, 'Tie::File', $Parameters{'input_file_name'} or die ('Error: Cannot Tie::File file: ' . $Parameters{'input_file_name'} . "\n");
my $INPUT_FILE_line_number = @INPUT_FILE ;

##################################################
#		 Main code
##################################################
# input file parsing :
for (my $i=0 ; $i<$INPUT_FILE_line_number; $i+=1) {
	# recover line
	my $line = $INPUT_FILE[$i] ;

	# Recovering parameters
	my $Parameters_ref = \%Parameters ;
	$Parameters_ref->{'current_alr_gen_line'} = $line ;
	$Parameters_ref->{'current_alr_gen_line_number'} = $i ;

	if (($line !~ m/\>/)&&($line !~ m/position/)) {

		# the line contains genotypes, retrieve them
		my %heterogametic_genotypes; # keys = heterogametic individuals' names, values = genotype
		$Parameters_ref->{'heterogametic_genotypes'} = \%heterogametic_genotypes ;
		my %homogametic_genotypes; # keys = homogametic individuals' names, values = genotype
		$Parameters_ref->{'homogametic_genotypes'} = \%homogametic_genotypes ;
		$Parameters_ref->{'homogametic_parent_genotype'} = "" ;
		$Parameters_ref->{'heterogametic_parent_genotype'} = "" ;
		retrieve_reads2snp_genotypes($Parameters_ref) ;

		# if line contains defined parental genotypes, study it
		if (($Parameters_ref->{'homogametic_parent_genotype'} ne "NN")&&($Parameters_ref->{'heterogametic_parent_genotype'} ne "NN")) {
			# segregation analysis
			add_line_to_table($Parameters_ref) ;
		}	
		
	} elsif ($line =~ m/position/) {
		# header contig line containing individuals names, retrieve their column number
		find_individuals_alr_gen_column_number($Parameters_ref) ;
	}
}

# printing output file
print(OUTPUT_FILE '>'.basename($Parameters{'input_file_name'})."_summary\n") ;
print(OUTPUT_FILE "number_time_line_happens\tSp|".$Parameters{'homogametic_parent_name'}."\tSp|".$Parameters{'heterogametic_parent_name'}) ;
for my $homogametic_individual (@homogametic) {
	print(OUTPUT_FILE "\tSp|".$homogametic_individual) ;
}
for my $heterogametic_individual (@heterogametic) {
	print(OUTPUT_FILE "\tSp|".$heterogametic_individual) ;
}
print(OUTPUT_FILE "\n") ;
for my $line (keys %Table) {
	print(OUTPUT_FILE $Table{$line}."\t".$line."\n") ;
}

##################################################
#	 closing input and output files
##################################################
close(INPUT_FILE) ;
close(OUTPUT_FILE) ;


##################################################
#	    	 Functions
##################################################

#---------------------------------------------------------------------------------
# Search column numbers for each individual in alr_gen file for the current contig
#---------------------------------------------------------------------------------
sub find_individuals_alr_gen_column_number {

	# Recovering parameters
	my $Parameters_ref = shift ;
	my $alr_gen_columns_ref = $Parameters_ref->{'alr_gen_columns'} ;
	my $heterogametic_ref = $Parameters_ref->{'heterogametic'} ;
	my $homogametic_ref = $Parameters_ref->{'homogametic'} ;

	# split line
	my @split_Alr_gen_line ;
	@split_Alr_gen_line = split(/\t/, $Parameters_ref->{'current_alr_gen_line'}) ;

	# print warning if more individuals than indicated are present for the given contig
#	if ($#split_Alr_gen_line > $#homogametic+$#heterogametic+6) {
#		print 'WARNING !!! There are individuals in contig ' . $Parameters_ref->{'contig_name'} . ' in alr_gen_ file that were not specified in command line, they will be ignored.' . "\n";
#	}

	# Search for each individual in the contig header and register their column number
	my @individuals = ($Parameters_ref->{'homogametic_parent_name'}, $Parameters_ref->{'heterogametic_parent_name'}, @{$homogametic_ref}, @{$heterogametic_ref}) ;
	foreach my $individual (@individuals) {
		if (!grep(/$individual$/, @split_Alr_gen_line)) {
			$alr_gen_columns_ref->{$individual} = 'NA' ;
		} else {
			foreach my $column (0 .. $#split_Alr_gen_line) {
				if ($split_Alr_gen_line[$column] =~ /$individual$/) {
					$alr_gen_columns_ref->{$individual} = $column ;
				}
			}
		}
	}

	return 0;
}


#----------------------------------------
# Retrieving reads2snp inferred genotypes 
#----------------------------------------
sub retrieve_reads2snp_genotypes {

	# Recovering parameters
	my $Parameters_ref = shift ;
	my $alr_gen_columns_ref = $Parameters_ref->{'alr_gen_columns'} ;
	my $heterogametic_genotypes_ref = $Parameters_ref->{'heterogametic_genotypes'} ;
	my $homogametic_genotypes_ref = $Parameters_ref->{'homogametic_genotypes'} ;
	my $heterogametic_ref = $Parameters_ref->{'heterogametic'} ;
	my $homogametic_ref = $Parameters_ref->{'homogametic'} ;

	# split line
	my @current_alr_gen_line_split = split(/\t/, $Parameters_ref->{'current_alr_gen_line'}) ;

	# Recovering genotypes for each individual
	# homogametic offspring
	foreach my $individual (@{$homogametic_ref}) {
		if ($alr_gen_columns_ref->{$individual} ne 'NA') {
			my $individual_genotype = $current_alr_gen_line_split[$alr_gen_columns_ref->{$individual}] ;
			$individual_genotype =~ s/\|\d{1}(\.\d+)?// ;
			$homogametic_genotypes_ref->{$individual} = $individual_genotype ;
		} else {
			$homogametic_genotypes_ref->{$individual} = 'NN' ;
		}
	}
	# heterogametic offspring
	foreach my $individual (@{$heterogametic_ref}) {
		if ($alr_gen_columns_ref->{$individual} ne 'NA') {
			my $individual_genotype = $current_alr_gen_line_split[$alr_gen_columns_ref->{$individual}] ;
			$individual_genotype =~ s/\|\d{1}(\.\d+)?// ;
			$heterogametic_genotypes_ref->{$individual} = $individual_genotype ;
		} else {
			$heterogametic_genotypes_ref->{$individual} = 'NN' ;
		}
	}
	# homogametic parent
	if ($alr_gen_columns_ref->{$Parameters_ref->{'homogametic_parent_name'}} ne 'NA') {
		$Parameters_ref->{'homogametic_parent_genotype'} = $current_alr_gen_line_split[$alr_gen_columns_ref->{$Parameters_ref->{'homogametic_parent_name'}}] ;
		$Parameters_ref->{'homogametic_parent_genotype'} =~ s/\|\d{1}(\.\d+)?// ;
	} else {
		$Parameters_ref->{'homogametic_parent_genotype'} = 'NN' ;
	}
	# heterogametic parent
	if ($alr_gen_columns_ref->{$Parameters_ref->{'heterogametic_parent_name'}} ne 'NA') {
		$Parameters_ref->{'heterogametic_parent_genotype'} = $current_alr_gen_line_split[$alr_gen_columns_ref->{$Parameters_ref->{'heterogametic_parent_name'}}] ;
		$Parameters_ref->{'heterogametic_parent_genotype'} =~ s/\|\d{1}(\.\d+)?// ;
	} else {
		$Parameters_ref->{'heterogametic_parent_genotype'} = 'NN' ;
	}

	return 0;
}

#-----------------------------------------------------
# Segregation analysis with parents and tagged progeny
#-----------------------------------------------------
sub add_line_to_table {

	# Recovering parameters
	my $Parameters_ref = shift ;
	my $heterogametic_genotypes_ref = $Parameters_ref->{'heterogametic_genotypes'} ;
	my $homogametic_genotypes_ref = $Parameters_ref->{'homogametic_genotypes'} ;

	# Do not consider heterogametic and homogametic with an undetermined genotype (NN)
	my $homogametic_length = scalar(keys %{$homogametic_genotypes_ref}) - grep(/NN/, values %{$homogametic_genotypes_ref}) ;
	my $heterogametic_length = scalar(keys %{$heterogametic_genotypes_ref}) - grep(/NN/, values %{$heterogametic_genotypes_ref}) ;

	# count if all genotypes are identical (monomorphic position)
	my $count = 0 ;
	if ($Parameters_ref->{'heterogametic_parent_genotype'} eq $Parameters_ref->{'homogametic_parent_genotype'}) {
		$count += 2 ;
		foreach my $genotype (values %{$heterogametic_genotypes_ref}) {
			if (($genotype ne 'NN')&&($genotype eq $Parameters_ref->{'heterogametic_parent_genotype'})) {
				$count ++ ;
			}
		}
		foreach my $genotype (values %{$homogametic_genotypes_ref}) {
			if (($genotype ne 'NN')&&($genotype eq $Parameters_ref->{'heterogametic_parent_genotype'})) {
				$count ++ ;
			}
		}
	}

	# Analysis of the line if at least 3 individuals of each sex in the progeny have a defined genotype and if position is not monomorphic
	if (($homogametic_length > 2)&&($heterogametic_length > 2)&&($count < ($homogametic_length + $heterogametic_length + 2))) {
		# writing line, ordering progeny genotypes...
		my $line_summary = $Parameters_ref->{'homogametic_parent_genotype'}."\t".$Parameters_ref->{'heterogametic_parent_genotype'} ;
		my @genoypes_order = ('AA', 'AC', 'AG', 'AT', 'CC', 'CG', 'CT', 'GG', 'GT', 'TT', 'NN') ;
		my %homogametic_genotypes_numbers = ('AA'=>0, 'AC'=>0, 'AG'=>0, 'AT'=>0, 'CC'=>0, 'CG'=>0, 'CT'=>0, 'GG'=>0, 'GT'=>0, 'TT'=>0, 'NN'=>0);
		my %heterogametic_genotypes_numbers = ('AA'=>0, 'AC'=>0, 'AG'=>0, 'AT'=>0, 'CC'=>0, 'CG'=>0, 'CT'=>0, 'GG'=>0, 'GT'=>0, 'TT'=>0, 'NN'=>0);
		for my  $genotype (values %{$homogametic_genotypes_ref}) {
			$homogametic_genotypes_numbers{$genotype} += 1 ;
		}
		for my  $genotype (values %{$heterogametic_genotypes_ref}) {
			$heterogametic_genotypes_numbers{$genotype} += 1 ;
		}
		for my  $genotype (@genoypes_order) {
			if ($homogametic_genotypes_numbers{$genotype} > 0) {
				for (my $i=0 ; $i<$homogametic_genotypes_numbers{$genotype}; $i+=1) {
					$line_summary .= "\t".$genotype ;
				}
			}
		}
		for my  $genotype (@genoypes_order) {
			if ($heterogametic_genotypes_numbers{$genotype} > 0) {
				for (my $i=0 ; $i<$heterogametic_genotypes_numbers{$genotype}; $i+=1) {
					$line_summary .= "\t".$genotype ;
				}
			}
		}

		# test if line already exists in %Table
		my $line_exists = "no" ;
		for my $line_in_Table (keys %Table) {
			if ($line_summary eq $line_in_Table) {
				$line_exists = "yes" ;
				$Table{$line_in_Table} += 1 ;
			}
		}
		if ($line_exists eq "no") {
			$Table{$line_summary} = 1 ;
		}
	}

	return 0;
}

#-----------------
# Check parameters
#-----------------
sub checkParameters {

	# Recovering parameters
	my $Parameters_ref = shift;
	my $heterogametic_ref = $Parameters_ref->{'heterogametic'} ;
	my $homogametic_ref = $Parameters_ref->{'homogametic'} ;
	my $homogametic_length = scalar(@{$homogametic_ref}) ;
	my $heterogametic_length = scalar(@{$heterogametic_ref}) ;
	

	# Check that individuals' names don't include special caracters
	if ((grep(/[\\\|\(\)\[\]\{\}\^\$\*\+\.\?]+/, @{$homogametic_ref}))||(grep(/[\\\|\(\)\[\]\{\}\^\$\*\+\.\?]+/, @{$heterogametic_ref}))||((defined($Parameters_ref->{'homogametic_parent_name'}))&&($Parameters_ref->{'homogametic_parent_name'} =~ /[\\\|\(\)\[\]\{\}\^\$\*\+\.\?]+/))||((defined($Parameters_ref->{'heterogametic_parent_name'}))&&($Parameters_ref->{'heterogametic_parent_name'} =~ /[\\\|\(\)\[\]\{\}\^\$\*\+\.\?]+/))) {
		print 'Error: please avoid using special caracters for individuals\' names.' . "\n" ;
		print 'Special caracters include \ | ( ) [ ] { } ^ $ * + ? . ' . "\n" ;
		print 'If columns in alr and alr_gen files are named after species_name|individual_name, only specify the individual name in command line' . "\n" ;
		exit(1) ;
	}
	# removing homogametic parent from @homogametic if necessary
	if ((@{$homogametic_ref})&&(defined($Parameters_ref->{'homogametic_parent_name'}))) {
		for (my $index=0 ; $index<$homogametic_length; $index+=1) {
			if ($Parameters_ref->{'homogametic_parent_name'} eq $homogametic_ref->[$index]) {
				delete($homogametic_ref->[$index]) ;
			}
		}
	}
	# removing heterogametic parent from @heterogametic if necessary
	if ((@{$heterogametic_ref})&&(defined($Parameters_ref->{'heterogametic_parent_name'}))) {
		for (my $index=0 ; $index<$heterogametic_length; $index+=1) {
			if ($Parameters_ref->{'heterogametic_parent_name'} eq $heterogametic_ref->[$index]) {
				delete($heterogametic_ref->[$index]) ;
			}
		}
	}
	# checking @heterogametic and @homogametic don't have individuals in common
	foreach my $individual (@{$heterogametic_ref}) {
		my $number_occurence = 0 ;
		foreach my $individual2 (@{$homogametic_ref}) {
			if ($individual eq $individual2) {
				$number_occurence ++;
			}
		}
		if ($number_occurence > 0) {
			print 'Error: individual ' . $individual . ' was assigned as homogametic and heterogametic' ."\n"; ;
			exit(1);
		}
	}
	# checking @heterogametic and @homogametic individuals are unique
	foreach my $individual (@{$heterogametic_ref}) {
		my $number_occurence = 0 ;
		foreach my $individual2 (@{$heterogametic_ref}) {
			if ($individual eq $individual2) {
				$number_occurence ++;
			}
		}
		if ($number_occurence > 1) {
			print 'Error: individual ' . $individual . ' was called multiple times' ."\n"; ;
			exit(1);
		}
	}
	foreach my $individual (@{$homogametic_ref}) {
		my $number_occurence = 0 ;
		foreach my $individual2 (@{$homogametic_ref}) {
			if ($individual eq $individual2) {
				$number_occurence ++;
			}
		}
		if ($number_occurence > 1) {
			print 'Error: individual ' . $individual . ' was called multiple times' ."\n"; ;
			exit(1);
		}
	}

 	# Check that if a homogametic parent is specified a heterogametic parent is specified as well and conversely
	if (((defined($Parameters_ref->{'homogametic_parent_name'}))&&(!defined($Parameters_ref->{'heterogametic_parent_name'})))||((!defined($Parameters_ref->{'homogametic_parent_name'}))&&(defined($Parameters_ref->{'heterogametic_parent_name'})))) {
		print 'Error: If a homogametic parent is specified you need to specify a heterogametic parent as well and conversely!' . "\n";
		exit(1) ;
	}
}



