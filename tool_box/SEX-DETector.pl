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

# example of command line : 
#./SEX-DETector/SEX-DETector.pl -alr jeu_test.alr -alr_gen jeu_test.alr_gen -alr_gen_sum jeu_test.alr_gen_summary -out test/jeu_test -hom C1_26_female,C1_27_female,C1_29_female,C1_34_female -het C1_01_male,C1_3_male,C1_04_male,C1_05_male -hom_par U10_37_mother -het_par Leuk144-3_father -seq -detail -detail-sex-linked -L -system xy -thr 0.8

##################################################
# 		  Modules
##################################################
# Perl modules
use strict ;
use warnings ;
use diagnostics ;
use Getopt::Long ;
use File::Basename;
use Tie::File;
use List::Util 'max'; 

##################################################
#	    variables initialization
##################################################
my $Start_time_hr = localtime ;
my $help = undef() ;
my %Parameters ;
my $homogametic_input;
my $heterogametic_input;
my %alr_columns ; # keys = individuals' names, values = column number in alr file
$Parameters{'alr_columns'} = \%alr_columns ;
my %alr_gen_columns ; # keys = individuals' names, values = column number in alr_gen file
$Parameters{'alr_gen_columns'} = \%alr_gen_columns ;
my %alpha; # keys = homogametic parent true genotype (and heterogametic with segregation type j=1), values = probability
$Parameters{'alpha'} = \%alpha ;
my %beta_2; # keys = heterogametic parent true genotype (segregation j=2), values = probability
$Parameters{'beta_2'} = \%beta_2 ;
my %beta_3; # keys = heterogametic parent true genotype (segregation j=3), values = probability
$Parameters{'beta_3'} = \%beta_3 ;
my %parents_genotypes_numbers_dataset = ('AA'=>0, 'AC'=>0, 'AG'=>0, 'AT'=>0, 'CC'=>0, 'CG'=>0, 'CT'=>0, 'GG'=>0, 'GT'=>0, 'TT'=>0) ;
$Parameters{'parents_genotypes_numbers_dataset'} = \%parents_genotypes_numbers_dataset ;
$Parameters{'non_optimized_parameters'} = 0 ;

##################################################
#	      options management 
##################################################
GetOptions (
	'help|h'				=>	\$help,
	'L'					=>	\$Parameters{'compute_Likelihood'},
	'debug'					=>	\$Parameters{'debug'},
	'no_sex_chr'				=>	\$Parameters{'no_sex_chr'},
	'SEM'					=>	\$Parameters{'SEM'},
	'skip_opt'				=>	\$Parameters{'skip_optimization'},
	'alrfile|alr=s'				=>	\$Parameters{'ALR_file'},
	'alr_gen_file|alr_gen=s'		=>	\$Parameters{'ALR_gen_file'},
	'alr_gen_summary_file|alr_gen_sum=s'	=>	\$Parameters{'ALR_gen_summary_file'},
	'output|out=s'				=>	\$Parameters{'output_file'},
	'system|syst=s'				=>	\$Parameters{'system'},
	'homogametic|hom=s'			=>	\$homogametic_input,
	'heterogametic|het=s'			=>	\$heterogametic_input,
	'homogametic_parent|hom_par=s'		=>	\$Parameters{'homogametic_parent_name'},
	'heterogametic_parent|het_par=s'	=>	\$Parameters{'heterogametic_parent_name'},
	'sequences|seq'				=>	\$Parameters{'sequences'},
	'detail-sex-linked|det-sex'		=>	\$Parameters{'detail-sex-linked'},
	'detail|det'				=>	\$Parameters{'detail'},
	'p=s'					=>	\$Parameters{'p'},
	'pi_1=s'				=>	\$Parameters{'pi_1'},
	'pi_2=s'				=>	\$Parameters{'pi_2'},
	'pi_3=s'				=>	\$Parameters{'pi_3'},
	'E=s'					=>	\$Parameters{'E'},
	'thr=s'					=>	\$Parameters{'threshold'},
	'param=s' 				=>	\$Parameters{'control_file'},
);

my @homogametic = split(/,/, $homogametic_input) ; # contains homogametic progeny individual names
$Parameters{'homogametic'} = \@homogametic ; 
my @heterogametic = split(/,/, $heterogametic_input) ; # contains heterogametic progeny individual names
$Parameters{'heterogametic'} = \@heterogametic ;

# Display help message if asked
if (defined($help)) {
	displayHelpMessage();
	exit(1) ;
}

# test if output directory exists, stop otherwise
my $output_dir = dirname($Parameters{'output_file'}) ;
print("\n\n$output_dir\n\n") ;
if (-d $output_dir) {
	#ok
} else {
	mkdir $output_dir;
	if (-d $output_dir) {
		# problem solved!
	} else {
	    die ('Error: output directory ' . $output_dir . " does not exist and cannot be created!\n");
	}
}

# retrieve parameters in control file
if (defined($Parameters{'control_file'})) {
	read_control_file(\%Parameters) ;
}

# Check parameters
checkParameters(\%Parameters) ;

# Generate output files names and headers
generateOutputfileNames(\%Parameters);

##################################################
#	           Main Code
##################################################

# Welcome message
print("\n");
print "###########################################################\n";
print "# 	  Welcome in SEX-DETector (Version $version)\n"; 
print "###########################################################\n";

# Verbose
if (defined($Parameters{'ALR_file'})) {
	print 'Selected alr input file : ' . $Parameters{'ALR_file'} . "\n";
}
print 'Selected alr_gen input file : ' . $Parameters{'ALR_gen_file'} . "\n";
print 'Chosen heterogametic system : ' . $Parameters{'system'} . "\n";
print 'homogametic parent name : ' . $Parameters{'homogametic_parent_name'} . "\n"; 
print 'heterogametic parent name : ' . $Parameters{'heterogametic_parent_name'} . "\n"; 
print "homogametic progeny names : @homogametic \n";
print "heterogametic progeny names : @heterogametic \n";
print 'initial Y genotyping error rate (p) : ' . $Parameters{'p'}."\n";
$Parameters{'p_fixed_to_zero'} = "no" ;
if ($Parameters{'p'} == 0) {
	$Parameters{'p_fixed_to_zero'} = "yes" ;
}
print 'initial genotyping error rate (epsilon) : ' . $Parameters{'E'}."\n";
my $pi_ref = $Parameters{'pi'} ;
print("initial autosomal probability : ".$pi_ref->{1}."\n"."initial XY or ZW probability : ".$pi_ref->{2}."\n"."initial hemizygous probability : ".$pi_ref->{3}."\n") ;

# Parameters initializations
# --------------------------
# Parents genotypes probabilities (computed) from the dataset
$Parameters{'alpha_beta_initializations_over'} = "no" ;
$Parameters{'parameters_estimations_over'} = "no" ;
if ((defined($Parameters{'alpha_AA'}))&&(defined($Parameters{'alpha_AC'}))&&(defined($Parameters{'alpha_AG'}))&&(defined($Parameters{'alpha_AT'}))&&(defined($Parameters{'alpha_CC'}))&&(defined($Parameters{'alpha_CG'}))&&(defined($Parameters{'alpha_CT'}))&&(defined($Parameters{'alpha_GG'}))&&(defined($Parameters{'alpha_GT'}))&&(defined($Parameters{'alpha_TT'}))&&(defined($Parameters{'beta_2_AC'}))&&(defined($Parameters{'beta_2_AG'}))&&(defined($Parameters{'beta_2_AT'}))&&(defined($Parameters{'beta_2_CG'}))&&(defined($Parameters{'beta_2_CT'}))&&(defined($Parameters{'beta_2_GT'}))&&(defined($Parameters{'beta_2_CA'}))&&(defined($Parameters{'beta_2_GA'}))&&(defined($Parameters{'beta_2_TA'}))&&(defined($Parameters{'beta_2_GC'}))&&(defined($Parameters{'beta_2_TC'}))&&(defined($Parameters{'beta_2_TG'}))&&(defined($Parameters{'beta_3_A'}))&&(defined($Parameters{'beta_3_C'}))&&(defined($Parameters{'beta_3_G'}))&&(defined($Parameters{'beta_3_T'}))) {
	%alpha = ('AA'=>$Parameters{'alpha_AA'}, 'AC'=>$Parameters{'alpha_AC'}, 'AG'=>$Parameters{'alpha_AG'}, 'AT'=>$Parameters{'alpha_AT'}, 'CC'=>$Parameters{'alpha_CC'}, 'CG'=>$Parameters{'alpha_CG'}, 'CT'=>$Parameters{'alpha_CT'}, 'GG'=>$Parameters{'alpha_GG'}, 'GT'=>$Parameters{'alpha_GT'}, 'TT'=>$Parameters{'alpha_TT'}) ;
	%beta_2 = ('AC'=>$Parameters{'beta_2_AC'}, 'CA'=>$Parameters{'beta_2_CA'}, 'AG'=>$Parameters{'beta_2_AG'}, 'GA'=>$Parameters{'beta_2_GA'}, 'AT'=>$Parameters{'beta_2_AT'}, 'TA'=>$Parameters{'beta_2_TA'}, 'CG'=>$Parameters{'beta_2_CG'}, 'GC'=>$Parameters{'beta_2_GC'}, 'CT'=>$Parameters{'beta_2_CT'}, 'TC'=>$Parameters{'beta_2_TC'}, 'GT'=>$Parameters{'beta_2_GT'}, 'TG'=>$Parameters{'beta_2_TG'}) ;
	%beta_3 = ('A'=>$Parameters{'beta_3_A'}, 'C'=>$Parameters{'beta_3_C'}, 'G'=>$Parameters{'beta_3_G'}, 'T'=>$Parameters{'beta_3_T'}) ;
} else {
	print "\n";
	print "Initializing parents genotypes frequencies from the dataset (parameters alpha and beta)\n" ;
	Expectation_Step(\%Parameters) ;
	initialize_alpha_and_beta(\%Parameters) ;
}
$Parameters{'alpha_beta_initializations_over'} = "yes" ;
print "alpha :\n";
for my $key (keys %alpha) {
	print($key."\t".$alpha{$key}."\n") ;
}
if (!defined($Parameters{'no_sex_chr'})) {
	print "beta_2 :\n";
	for my $key (keys %beta_2) {
		print($key."\t".$beta_2{$key}."\n") ;
	}
	print "beta_3 :\n";
	for my $key (keys %beta_3) {
		print($key."\t".$beta_3{$key}."\n") ;
	}
}
print "\n";

# printing current parameters values
open(PARAM, '>>'.$Parameters{'output_parameters'}) ;
print(PARAM "\nInitialization\t".$pi_ref->{1}."\t".$pi_ref->{2}."\t".$pi_ref->{3}."\t".$Parameters{'p'}."\t".$Parameters{'E'}."\t".$alpha{'AA'}."\t".$alpha{'AC'}."\t".$alpha{'AG'}."\t".$alpha{'AT'}."\t".$alpha{'CC'}."\t".$alpha{'CG'}."\t".$alpha{'CT'}."\t".$alpha{'GG'}."\t".$alpha{'GT'}."\t".$alpha{'TT'}."\t".$beta_2{'AC'}."\t".$beta_2{'CA'}."\t".$beta_2{'AG'}."\t".$beta_2{'GA'}."\t".$beta_2{'AT'}."\t".$beta_2{'TA'}."\t".$beta_2{'CG'}."\t".$beta_2{'GC'}."\t".$beta_2{'CT'}."\t".$beta_2{'TC'}."\t".$beta_2{'GT'}."\t".$beta_2{'TG'}."\t".$beta_3{'A'}."\t".$beta_3{'C'}."\t".$beta_3{'G'}."\t".$beta_3{'T'}) ;
close(PARAM) ;																																						    									   														        

# registering old values of parameters
$Parameters{'alpha_old'} = \%alpha ;
$Parameters{'beta_2_old'} = \%beta_2 ;
$Parameters{'beta_3_old'} = \%beta_3 ;
$Parameters{'pi_old'} = $pi_ref ;
$Parameters{'p_old'} = $Parameters{'p'} ;
$Parameters{'E_old'} = $Parameters{'E'} ;

# Generating segregation tables
create_segregation_table_tagged_progeny(\%Parameters) ;

# initializations for the next iteration
my %alpha_new = ('AA'=>0, 'AC'=>0, 'AG'=>0, 'AT'=>0, 'CC'=>0, 'CG'=>0, 'CT'=>0, 'GG'=>0, 'GT'=>0, 'TT'=>0);
$Parameters{'alpha_new'} = \%alpha_new ;
my %beta_2_new = ('AC'=>0, 'CA'=>0, 'AG'=>0, 'GA'=>0, 'AT'=>0, 'TA'=>0, 'CG'=>0, 'GC'=>0, 'CT'=>0, 'TC'=>0, 'GT'=>0, 'TG'=>0);
$Parameters{'beta_2_new'} = \%beta_2_new ;
my %beta_3_new = ('A'=>0, 'C'=>0, 'G'=>0, 'T'=>0);
$Parameters{'beta_3_new'} = \%beta_3_new ;
$Parameters{'p_new_numerator'} = 0 ;
$Parameters{'E_new_numerator'} = 0 ;
$Parameters{'sum_k_t'} = 0 ;
my %sum_k_t_S_k_t = (1=>0,2=>0,3=>0) ;
$Parameters{'sum_k_t_S_k_t'} = \%sum_k_t_S_k_t ;
$Parameters{'p_new_denominator'} = 0 ;
$Parameters{'E_new_denominator'} = 0 ;
$Parameters{'Q'} = 0 ;
$Parameters{'H'} = 0 ;
my %TMG_1_counts = ('AA'=>0, 'AC'=>0, 'AG'=>0, 'AT'=>0, 'CC'=>0, 'CG'=>0, 'CT'=>0, 'GG'=>0, 'GT'=>0, 'TT'=>0);
$Parameters{'TMG_1_counts'} = \%TMG_1_counts ;
my %TMG_2_counts = ('AA'=>0, 'AC'=>0, 'AG'=>0, 'AT'=>0, 'CC'=>0, 'CG'=>0, 'CT'=>0, 'GG'=>0, 'GT'=>0, 'TT'=>0);
$Parameters{'TMG_2_counts'} = \%TMG_2_counts ;
my %TMG_3_counts = ('AA'=>0, 'AC'=>0, 'AG'=>0, 'AT'=>0, 'CC'=>0, 'CG'=>0, 'CT'=>0, 'GG'=>0, 'GT'=>0, 'TT'=>0);
$Parameters{'TMG_3_counts'} = \%TMG_3_counts ;
my %TFG_1_counts = ('AA'=>0, 'AC'=>0, 'AG'=>0, 'AT'=>0, 'CC'=>0, 'CG'=>0, 'CT'=>0, 'GG'=>0, 'GT'=>0, 'TT'=>0);
$Parameters{'TFG_1_counts'} = \%TFG_1_counts ;
my %TFG_2_counts = ('AC'=>0, 'CA'=>0, 'AG'=>0, 'GA'=>0, 'AT'=>0, 'TA'=>0, 'CG'=>0, 'GC'=>0, 'CT'=>0, 'TC'=>0, 'GT'=>0, 'TG'=>0);
$Parameters{'TFG_2_counts'} = \%TFG_2_counts ;
my %TFG_3_counts = ('A'=>0, 'C'=>0, 'G'=>0, 'T'=>0);
$Parameters{'TFG_3_counts'} = \%TFG_3_counts ;
my %S_counts = (1=>0,2=>0,3=>0) ;
$Parameters{'S_counts'} = \%S_counts ;
$Parameters{'sample_size'} = 0 ;

if (!defined($Parameters{'skip_optimization'})) {
	# alr and alr_gen files parsing in order to estimate new values of parameters
	$Parameters{'run'} = 1 ;
	print "EM algorithm iteration ".$Parameters{'run'}."\n";
	Expectation_Step(\%Parameters) ;

	# estimations of new parameters values and comparison with old ones to test convergence
	my $convergence = Maximization_Step(\%Parameters) ;

	# further runs of parameters estimations if necessary until values converge
	while ($convergence eq "no") {

		# initializations for the next iteration
		%alpha_new = ('AA'=>0, 'AC'=>0, 'AG'=>0, 'AT'=>0, 'CC'=>0, 'CG'=>0, 'CT'=>0, 'GG'=>0, 'GT'=>0, 'TT'=>0);
		%beta_2_new = ('AC'=>0, 'CA'=>0, 'AG'=>0, 'GA'=>0, 'AT'=>0, 'TA'=>0, 'CG'=>0, 'GC'=>0, 'CT'=>0, 'TC'=>0, 'GT'=>0, 'TG'=>0);
		%beta_3_new = ('A'=>0, 'C'=>0, 'G'=>0, 'T'=>0);
		$Parameters{'p_new_numerator'} = 0 ;
		$Parameters{'E_new_numerator'} = 0 ;
		$Parameters{'sum_k_t'} = 0 ;
		%sum_k_t_S_k_t = (1=>0,2=>0,3=>0) ;
		$Parameters{'p_new_denominator'} = 0 ;
		$Parameters{'E_new_denominator'} = 0 ;
		$Parameters{'Q'} = 0 ;
		$Parameters{'H'} = 0 ;
		%TMG_1_counts = ('AA'=>0, 'AC'=>0, 'AG'=>0, 'AT'=>0, 'CC'=>0, 'CG'=>0, 'CT'=>0, 'GG'=>0, 'GT'=>0, 'TT'=>0);
		%TMG_2_counts = ('AA'=>0, 'AC'=>0, 'AG'=>0, 'AT'=>0, 'CC'=>0, 'CG'=>0, 'CT'=>0, 'GG'=>0, 'GT'=>0, 'TT'=>0);
		%TMG_3_counts = ('AA'=>0, 'AC'=>0, 'AG'=>0, 'AT'=>0, 'CC'=>0, 'CG'=>0, 'CT'=>0, 'GG'=>0, 'GT'=>0, 'TT'=>0);
		%TFG_1_counts = ('AA'=>0, 'AC'=>0, 'AG'=>0, 'AT'=>0, 'CC'=>0, 'CG'=>0, 'CT'=>0, 'GG'=>0, 'GT'=>0, 'TT'=>0);
		%TFG_2_counts = ('AC'=>0, 'CA'=>0, 'AG'=>0, 'GA'=>0, 'AT'=>0, 'TA'=>0, 'CG'=>0, 'GC'=>0, 'CT'=>0, 'TC'=>0, 'GT'=>0, 'TG'=>0);
		%TFG_3_counts = ('A'=>0, 'C'=>0, 'G'=>0, 'T'=>0);
		%S_counts = (1=>0,2=>0,3=>0) ;
		$Parameters{'sample_size'} = 0 ;
		$Parameters{'information'} = 0 ;

		# alr and alr_gen files parsing in order to estimate expectations of hidden variables
		print "EM algorithm iteration ".$Parameters{'run'}."\n";
		Expectation_Step(\%Parameters) ;

		# estimations of new parameters values and comparison with old ones to test convergence
		$convergence = Maximization_Step(\%Parameters) ;
	}

	$Parameters{'parameters_estimations_over'} = "yes" ;
	print "\nEstimation of parameters done.\n\n";
} else {
	$Parameters{'parameters_estimations_over'} = "yes" ;
	$Parameters{'run'} = 1 ;
	print "\nNo estimation of parameters.\n\n";
}

# initializations for the next iteration
%alpha_new = ('AA'=>0, 'AC'=>0, 'AG'=>0, 'AT'=>0, 'CC'=>0, 'CG'=>0, 'CT'=>0, 'GG'=>0, 'GT'=>0, 'TT'=>0);
%beta_2_new = ('AC'=>0, 'CA'=>0, 'AG'=>0, 'GA'=>0, 'AT'=>0, 'TA'=>0, 'CG'=>0, 'GC'=>0, 'CT'=>0, 'TC'=>0, 'GT'=>0, 'TG'=>0);
%beta_3_new = ('A'=>0, 'C'=>0, 'G'=>0, 'T'=>0);
$Parameters{'p_new_numerator'} = 0 ;
$Parameters{'E_new_numerator'} = 0 ;
$Parameters{'sum_k_t'} = 0 ;
%sum_k_t_S_k_t = (1=>0,2=>0,3=>0) ;
$Parameters{'p_new_denominator'} = 0 ;
$Parameters{'E_new_denominator'} = 0 ;
$Parameters{'Q'} = 0 ;
$Parameters{'H'} = 0 ;
%TMG_1_counts = ('AA'=>0, 'AC'=>0, 'AG'=>0, 'AT'=>0, 'CC'=>0, 'CG'=>0, 'CT'=>0, 'GG'=>0, 'GT'=>0, 'TT'=>0);
%TMG_2_counts = ('AA'=>0, 'AC'=>0, 'AG'=>0, 'AT'=>0, 'CC'=>0, 'CG'=>0, 'CT'=>0, 'GG'=>0, 'GT'=>0, 'TT'=>0);
%TMG_3_counts = ('AA'=>0, 'AC'=>0, 'AG'=>0, 'AT'=>0, 'CC'=>0, 'CG'=>0, 'CT'=>0, 'GG'=>0, 'GT'=>0, 'TT'=>0);
%TFG_1_counts = ('AA'=>0, 'AC'=>0, 'AG'=>0, 'AT'=>0, 'CC'=>0, 'CG'=>0, 'CT'=>0, 'GG'=>0, 'GT'=>0, 'TT'=>0);
%TFG_2_counts = ('AC'=>0, 'CA'=>0, 'AG'=>0, 'GA'=>0, 'AT'=>0, 'TA'=>0, 'CG'=>0, 'GC'=>0, 'CT'=>0, 'TC'=>0, 'GT'=>0, 'TG'=>0);
%TFG_3_counts = ('A'=>0, 'C'=>0, 'G'=>0, 'T'=>0);
%S_counts = (1=>0,2=>0,3=>0) ;
$Parameters{'sample_size'} = 0 ;

# attribution of contigs to categories and printing results in output file
print "Assigning contigs to segregation categories and writing output(s)...\n\n";
Expectation_Step(\%Parameters) ;

# estimations of new parameters values and final Likelihood of the model
my $convergence = Maximization_Step(\%Parameters) ;


print "\n";
# output start and end hour :
my $End_time_hr = localtime ;
print("Computation started : \t\t$Start_time_hr\n") ;
print("Computation ended : \t\t$End_time_hr\n") ;
print "###########################################################\n";
print "#                     End of execution                     \n"; 
print "###########################################################\n";
print("\n");


##################################################
#	    	 Functions
##################################################

#------------------------------
# alr and alr_gen files parsing
#------------------------------
sub Expectation_Step {

	# Recovering parameters
	my $Parameters_ref = shift ;

	# Initialization
	$Parameters_ref->{'contig_name'} = "" ;

	# Tie::File input files
	my @alr_line ;
	my @alr_gen_line ;
	if ($Parameters_ref->{'parameters_estimations_over'} eq "yes") {
		# use the whole alr_gen file to infer contig status
		tie @alr_gen_line, 'Tie::File', $Parameters_ref->{'ALR_gen_file'} or die ('Error: Cannot Tie::File file: ' . $Parameters_ref->{'ALR_gen_file'} . "\n");
		if (defined($Parameters_ref->{'ALR_file'})) {
			tie @alr_line, 'Tie::File', $Parameters_ref->{'ALR_file'} or die ('Error: Cannot Tie::File file: ' . $Parameters_ref->{'ALR_file'} . "\n");
		}
	} else {
		# estimation of parameters, use the alr_gen_summary file to go faster!!!
		tie @alr_gen_line, 'Tie::File', $Parameters_ref->{'ALR_gen_summary_file'} or die ('Error: Cannot Tie::File file: ' . $Parameters_ref->{'ALR_gen_summary_file'} . "\n");
	}
	# Check that alr and alr_gen files are of same length
	my $alr_line_number = @alr_gen_line ;
	if (defined($Parameters_ref->{'ALR_file'})) {
		my $alr_gen_line_number = @alr_gen_line ;
		if ($alr_line_number!= $alr_gen_line_number) {
			print "\n\n" . 'Error: alr and alr_gen files do not have the same number of line!' . "\n\n" ;
			exit(1) ;		
		}
	}

	# parsing each line of both input files to study SNPs
	for (my $line_number=0 ; $line_number<$alr_line_number; $line_number+=1) {
		$Parameters_ref->{'current_alr_line'} = $alr_line[$line_number] ;
		$Parameters_ref->{'current_alr_gen_line'} = $alr_gen_line[$line_number] ;

		if ($Parameters_ref->{'current_alr_gen_line'} =~ m/^>/) {
			# the line corresponds to a new contig
			# if necessary analyze all the results of the previous contig
			if ($Parameters_ref->{'contig_name'} ne "") {
				if ($Parameters_ref->{'parameters_estimations_over'} eq "yes") {
					analyze_contig_results(\%Parameters) ;
				}
			}

			# Register new contig name
			$Parameters_ref->{'contig_name'} = $Parameters_ref->{'current_alr_gen_line'} ;
			$Parameters_ref->{'contig_name'} =~ s/>// ;

			# initializations for the new contig
			if ($Parameters_ref->{'alpha_beta_initializations_over'} eq "yes") {
				$Parameters_ref->{'SNP_calling_sex_linked_number'} = 0 ;
				$Parameters_ref->{'S_k_1'} = 0 ;
				$Parameters_ref->{'S_k_2'} = 0 ;
				$Parameters_ref->{'S_k_3'} = 0 ;
				$Parameters_ref->{'autosomal_without_error'} = 0 ;
				$Parameters_ref->{'autosomal_with_error'} = 0 ;
				$Parameters_ref->{'XY_without_error'} = 0 ;
				$Parameters_ref->{'XY_with_error'} = 0 ;
				$Parameters_ref->{'hemizygous_without_error'} = 0 ;
				$Parameters_ref->{'hemizygous_with_error'} = 0 ;
				if ($Parameters{'parameters_estimations_over'} eq "yes") {
					if (defined($Parameters_ref->{'sequences'})) {
						$Parameters_ref->{'current_contig_sequence_X1'} = "" ;
						$Parameters_ref->{'current_contig_sequence_X2'} = "" ;
						$Parameters_ref->{'current_contig_sequence_X3'} = "" ;
						$Parameters_ref->{'current_contig_sequence_Y'} = "" ;
					}
					if (defined($Parameters_ref->{'ALR_file'})) {
						$Parameters_ref->{'number_clean_XY_SNPs_without_error'} = 0 ;
						$Parameters_ref->{'number_clean_hemizygous_SNPs_without_error'} = 0 ;
					}
					if (defined($Parameters_ref->{'detail-sex-linked'})) {
						$Parameters_ref->{'detail_sex-linked_SNPs'} = "" ;
					}
				}
			}

		} elsif (($Parameters_ref->{'current_alr_gen_line'} =~ m/^position/)||($Parameters_ref->{'current_alr_gen_line'} =~ m/^number_time_line_happens/)) {
			# the line corresponds to headers, search the colomn number of each individual
			find_individuals_alr_gen_column_number($Parameters_ref) ;
			if (($Parameters_ref->{'parameters_estimations_over'} eq "yes")&&(defined($Parameters_ref->{'ALR_file'}))) {
				find_individuals_alr_column_number($Parameters_ref) ;
			}
		} else {
			# the line corresponds to a position in a contig

			# retrieve current contig position
			my @split_Alr_gen_line = split(/\t/, $Parameters_ref->{'current_alr_gen_line'}) ;
			$Parameters_ref->{'current_contig_position'} = $split_Alr_gen_line[0] ;

			# initialise genotypes
			my %heterogametic_genotypes; # keys = heterogametic individuals' names, values = genotype
			$Parameters_ref->{'heterogametic_genotypes'} = \%heterogametic_genotypes ;
			my %homogametic_genotypes; # keys = homogametic individuals' names, values = genotype
			$Parameters_ref->{'homogametic_genotypes'} = \%homogametic_genotypes ;
			$Parameters_ref->{'homogametic_parent_genotype'} = "" ;
			$Parameters_ref->{'heterogametic_parent_genotype'} = "" ;
			# retrieving reads2snp inferred genotypes
			retrieve_reads2snp_genotypes($Parameters_ref) ;

			if ($Parameters_ref->{'alpha_beta_initializations_over'} eq "yes") {
				# Check that both parents have a defined genotype (not NN)
				if (($Parameters_ref->{'homogametic_parent_genotype'} ne "NN")&&($Parameters_ref->{'heterogametic_parent_genotype'} ne "NN")) {
					# segregation analysis
					segregation_analysis_parents_tagged_progeny($Parameters_ref) ;
				} else {
					if ( (defined($Parameters_ref->{'ALR_file'}))&&(defined($Parameters_ref->{'sequences'}))&&($Parameters{'parameters_estimations_over'} eq "yes")) {
						# add base to sequences
						my @current_alr_line_split = split(/\t/, $Parameters_ref->{'current_alr_line'}) ;
						$Parameters_ref->{'current_contig_sequence_X1'} .= $current_alr_line_split[0] ;
						$Parameters_ref->{'current_contig_sequence_X2'} .= $current_alr_line_split[0] ;
						$Parameters_ref->{'current_contig_sequence_X3'} .= $current_alr_line_split[0] ;
						$Parameters_ref->{'current_contig_sequence_Y'} .= $current_alr_line_split[0] ;
					}
				}
			}
		}
	}

	# analyze the results of the last contig if necesary
	if ($Parameters_ref->{'parameters_estimations_over'} eq "yes") {
		analyze_contig_results(\%Parameters) ;
	}

	return 0;
}


#---------------------------------------------------------------------------------
# Search column numbers for each individual in alr_gen file for the current contig
#---------------------------------------------------------------------------------
sub find_individuals_alr_gen_column_number {

	# Recovering parameters
	my $Parameters_ref = shift ;
	my $alr_gen_columns_ref = $Parameters_ref->{'alr_gen_columns'} ;
	my $heterogametic_ref = $Parameters_ref->{'heterogametic'} ;
	my $homogametic_ref = $Parameters_ref->{'homogametic'} ;
	my @heterogametic = @{$heterogametic_ref} ;
	my @homogametic = @{$homogametic_ref} ;

	# split line
	my @split_Alr_gen_line ;
	@split_Alr_gen_line = split(/\t/, $Parameters_ref->{'current_alr_gen_line'}) ;

	# print warning if more individuals than indicated are present for the given contig
#	if (scalar(@split_Alr_gen_line) > (scalar(@homogametic)+scalar(@heterogametic)+3)) {
#		print 'WARNING!!! There are individuals in contig ' . $Parameters_ref->{'contig_name'} . ' in alr_gen_ file that were not specified in command line, they will be ignored.' . "\n";
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


#-----------------------------------------------------------------------------------------
# Search column numbers for each individual in alr and alr_gen file for the current contig
#-----------------------------------------------------------------------------------------
sub find_individuals_alr_column_number {

	# Recovering parameters
	my $Parameters_ref = shift ;
	my $alr_columns_ref = $Parameters_ref->{'alr_columns'} ;
	my $heterogametic_ref = $Parameters_ref->{'heterogametic'} ;
	my $homogametic_ref = $Parameters_ref->{'homogametic'} ;
	my @heterogametic = @{$heterogametic_ref} ;
	my @homogametic = @{$homogametic_ref} ;

	# split line
	my @split_Alr_line = split(/\t/, $Parameters_ref->{'current_alr_line'}) ;

	# print warning if more individuals than indicated are present for the given contig
#	if (scalar(@split_Alr_line) > (scalar(@homogametic)+scalar(@heterogametic)+4)) {
#		print 'WARNING!!! There are individuals in contig ' . $Parameters_ref->{'contig_name'} . ' in alr file that were not specified in command line, they will be ignored.' . "\n";
#	}

	# Search for each individual in the contig header and register their column number
	my @individuals = ($Parameters_ref->{'homogametic_parent_name'}, $Parameters_ref->{'heterogametic_parent_name'}, @{$homogametic_ref}, @{$heterogametic_ref}) ;
	foreach my $individual (@individuals) {
		if (!grep(/$individual$/, @split_Alr_line)) {
			#print 'Individual ' . $individual . ' is absent from contig ' . $Parameters_ref->{'contig_name'} . " in alr file\n";
			$alr_columns_ref->{$individual} = 'NA' ;
		} else {
			foreach my $column (0 .. $#split_Alr_line) {
				if ($split_Alr_line[$column] =~ /$individual$/) {
					$alr_columns_ref->{$individual} = $column ;
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
	my $parents_genotypes_numbers_dataset_ref = $Parameters{'parents_genotypes_numbers_dataset'} ;

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
	# adding this genotype to counts in order to initiate alpha_m and beta_n_j
	if ($Parameters_ref->{'alpha_beta_initializations_over'} eq "no") {
		if ($Parameters_ref->{'homogametic_parent_genotype'} ne 'NN') {
			$parents_genotypes_numbers_dataset_ref->{$Parameters_ref->{'homogametic_parent_genotype'}} += $current_alr_gen_line_split[0] ;
		}
		if ($Parameters_ref->{'heterogametic_parent_genotype'} ne 'NN') {
			$parents_genotypes_numbers_dataset_ref->{$Parameters_ref->{'heterogametic_parent_genotype'}} += $current_alr_gen_line_split[0] ;
		}
	}

	return 0;
}


#----------------------------------
# Retrieving reads2snp reads counts
#----------------------------------
sub retrieve_reads2snp_expression_levels {

	# Recovering parameters
	my $Parameters_ref = shift ;
	my $alr_columns_ref = $Parameters_ref->{'alr_columns'} ;
	my $heterogametic_expression_ref = $Parameters_ref->{'heterogametic_expression'} ;
	my $homogametic_expression_ref = $Parameters_ref->{'homogametic_expression'} ;
	my $heterogametic_ref = $Parameters_ref->{'heterogametic'} ;
	my $homogametic_ref = $Parameters_ref->{'homogametic'} ;
	
	# split line
	my @current_alr_line_split = split(/\t/, $Parameters_ref->{'current_alr_line'}) ;
	
	# Recovering expression levels for each individual
	# homogametic parent
	if ($alr_columns_ref->{$Parameters_ref->{'homogametic_parent_name'}} ne 'NA') {
		$Parameters_ref->{'homogametic_parent_expression'} = $current_alr_line_split[$alr_columns_ref->{$Parameters_ref->{'homogametic_parent_name'}}] ;
		$Parameters_ref->{'homogametic_parent_expression'} =~ s/\d+\[// ;
		$Parameters_ref->{'homogametic_parent_expression'} =~ s/]// ;
	} else {
		$Parameters_ref->{'homogametic_parent_expression'} = '0/0/0/0' ;
	}
	# heterogametic parent
	if ($alr_columns_ref->{$Parameters_ref->{'heterogametic_parent_name'}} ne 'NA') {
		$Parameters_ref->{'heterogametic_parent_expression'} = $current_alr_line_split[$alr_columns_ref->{$Parameters_ref->{'heterogametic_parent_name'}}] ;
		$Parameters_ref->{'heterogametic_parent_expression'} =~ s/\d+\[// ;
		$Parameters_ref->{'heterogametic_parent_expression'} =~ s/]// ;
	} else {
		$Parameters_ref->{'heterogametic_parent_expression'} = '0/0/0/0' ;
	}
	# homogametic sex
	foreach my $individual (@{$homogametic_ref}) {
		if ($alr_columns_ref->{$individual} ne 'NA') {
			my $individual_expression = $current_alr_line_split[$alr_columns_ref->{$individual}] ;
			$individual_expression =~ s/\d+\[// ;
			$individual_expression =~ s/]// ;
			$homogametic_expression_ref->{$individual} = $individual_expression ;
		} else {
			$homogametic_expression_ref->{$individual} = 'NN' ;
		}
	}
	# heterogametic sex
	foreach my $individual (@{$heterogametic_ref}) {
		if ($alr_columns_ref->{$individual} ne 'NA') {
			my $individual_expression = $current_alr_line_split[$alr_columns_ref->{$individual}] ;
			$individual_expression =~ s/\d+\[// ;
			$individual_expression =~ s/]// ;
			$heterogametic_expression_ref->{$individual} = $individual_expression ;
		} else {
			$heterogametic_expression_ref->{$individual} = 'NN' ;
		}
	}

	return 0;
}


#-----------------------------------------------------
# Segregation analysis with parents and tagged progeny
#-----------------------------------------------------
sub segregation_analysis_parents_tagged_progeny {

	# Recovering parameters
	my $Parameters_ref = shift ;
	my $heterogametic_genotypes_ref = $Parameters_ref->{'heterogametic_genotypes'} ;
	my $homogametic_genotypes_ref = $Parameters_ref->{'homogametic_genotypes'} ;
	my $alpha_ref = $Parameters_ref->{'alpha'} ;
	my $beta_2_ref = $Parameters_ref->{'beta_2'} ;
	my $beta_3_ref = $Parameters_ref->{'beta_3'} ;
	my $pi_ref = $Parameters_ref->{'pi'} ;
	my $sum_k_t_S_k_t_ref = $Parameters_ref->{'sum_k_t_S_k_t'} ;
	my $alpha_new_ref = $Parameters_ref->{'alpha_new'} ;
	my $beta_2_new_ref = $Parameters_ref->{'beta_2_new'} ;
	my $beta_3_new_ref = $Parameters_ref->{'beta_3_new'} ;
	my $lambda_1_ref = $Parameters_ref->{'lambda_1'} ;
	my $lambda_2_ref = $Parameters_ref->{'lambda_2'} ;
	my $lambda_3_ref = $Parameters_ref->{'lambda_3'} ;
	my $TMG_1_counts_ref = $Parameters_ref->{'TMG_1_counts'} ;
	my $TMG_2_counts_ref = $Parameters_ref->{'TMG_2_counts'} ;
	my $TMG_3_counts_ref = $Parameters_ref->{'TMG_3_counts'} ;
	my $TFG_1_counts_ref = $Parameters_ref->{'TFG_1_counts'} ;
	my $TFG_2_counts_ref = $Parameters_ref->{'TFG_2_counts'} ;
	my $TFG_3_counts_ref = $Parameters_ref->{'TFG_3_counts'} ;
	my $S_counts_ref = $Parameters_ref->{'S_counts'} ;

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

	# Analysis of segregation using genotypes inferred by Reads2snp if at least 3 individuals of each sex in the progeny have a defined genotype and if position is not monomorphic
	if (($homogametic_length > 2)&&($heterogametic_length > 2)&&($count < ($homogametic_length + $heterogametic_length + 2))) {

		# Initializations
		my @current_alr_gen_line_split = split(/\t/, $Parameters_ref->{'current_alr_gen_line'}) ;
		my $line_occurence_number ;
		if ($Parameters{'parameters_estimations_over'} eq "no") {
			$line_occurence_number = $current_alr_gen_line_split[0] ;
		} else {
			$line_occurence_number = 1 ;
		}
		my %individuals ;
		my $individuals_ref = \%individuals ;
		$individuals{$Parameters_ref->{'homogametic_parent_name'}} = 0 ;
		$individuals{$Parameters_ref->{'heterogametic_parent_name'}} = 0 ;
		my %het_individuals ;
		my $het_individuals_ref = \%het_individuals ;
		$het_individuals{$Parameters_ref->{'heterogametic_parent_name'}} = 0 ;
		foreach my $ind (keys %{$heterogametic_genotypes_ref}) {
			if ($heterogametic_genotypes_ref->{$ind} ne 'NN') {
				$individuals{$ind} = 0 ;
				$het_individuals{$ind} = 0 ;
			}
		}
		foreach my $ind (keys %{$homogametic_genotypes_ref}) {
			if ($homogametic_genotypes_ref->{$ind} ne 'NN') {
				$individuals{$ind} = 0 ;
			}
		}
		$Parameters_ref->{'individuals'} = $individuals_ref ;
		$Parameters_ref->{'het_individuals'} = $het_individuals_ref ;

		# increment sample size
		$Parameters_ref->{'sample_size'} += $line_occurence_number ;

		# First iteration
		my ($GE_k_t_1_ref, $GE_k_t_2_ref, $GE_k_t_3_ref, $YGE_k_t_2_ref, $factor_ref, $S_k_t_ref, $TMG_k_t_1_ref, $TMG_k_t_2_ref, $TMG_k_t_3_ref, $TFG_k_t_1_ref, $TFG_k_t_2_ref, $TFG_k_t_3_ref) ;
		if ((defined($Parameters_ref->{'SEM'}))&&($Parameters_ref->{'run'} == 1)&&(!defined($Parameters{'skip_optimization'}))) {
			# S step
			($S_k_t_ref->{1}, $S_k_t_ref->{2}, $S_k_t_ref->{3}) = random_multinomial_3($pi_ref->{1}, $pi_ref->{2}, $pi_ref->{3}) ;
			# TMG_k_t_1
			if ($S_k_t_ref->{1}!= 0) {
				($TMG_k_t_1_ref->{'AA'}, $TMG_k_t_1_ref->{'AC'}, $TMG_k_t_1_ref->{'AG'}, $TMG_k_t_1_ref->{'AT'}, $TMG_k_t_1_ref->{'CC'}, $TMG_k_t_1_ref->{'CG'}, $TMG_k_t_1_ref->{'CT'}, $TMG_k_t_1_ref->{'GG'}, $TMG_k_t_1_ref->{'GT'}, $TMG_k_t_1_ref->{'TT'}) = random_multinomial_10($alpha_ref->{'AA'}, $alpha_ref->{'AC'}, $alpha_ref->{'AG'}, $alpha_ref->{'AT'}, $alpha_ref->{'CC'}, $alpha_ref->{'CG'}, $alpha_ref->{'CT'}, $alpha_ref->{'GG'}, $alpha_ref->{'GT'}, $alpha_ref->{'TT'}) ;
			} else {
				($TMG_k_t_1_ref->{'AA'}, $TMG_k_t_1_ref->{'AC'}, $TMG_k_t_1_ref->{'AG'}, $TMG_k_t_1_ref->{'AT'}, $TMG_k_t_1_ref->{'CC'}, $TMG_k_t_1_ref->{'CG'}, $TMG_k_t_1_ref->{'CT'}, $TMG_k_t_1_ref->{'GG'}, $TMG_k_t_1_ref->{'GT'}, $TMG_k_t_1_ref->{'TT'}) = (0,0,0,0,0,0,0,0,0,0) ;
			}
			# TMG_k_t_2
			if ($S_k_t_ref->{2}!= 0) {
				($TMG_k_t_2_ref->{'AA'}, $TMG_k_t_2_ref->{'AC'}, $TMG_k_t_2_ref->{'AG'}, $TMG_k_t_2_ref->{'AT'}, $TMG_k_t_2_ref->{'CC'}, $TMG_k_t_2_ref->{'CG'}, $TMG_k_t_2_ref->{'CT'}, $TMG_k_t_2_ref->{'GG'}, $TMG_k_t_2_ref->{'GT'}, $TMG_k_t_2_ref->{'TT'}) = random_multinomial_10($alpha_ref->{'AA'}, $alpha_ref->{'AC'}, $alpha_ref->{'AG'}, $alpha_ref->{'AT'}, $alpha_ref->{'CC'}, $alpha_ref->{'CG'}, $alpha_ref->{'CT'}, $alpha_ref->{'GG'}, $alpha_ref->{'GT'}, $alpha_ref->{'TT'}) ;
			} else {
				($TMG_k_t_2_ref->{'AA'}, $TMG_k_t_2_ref->{'AC'}, $TMG_k_t_2_ref->{'AG'}, $TMG_k_t_2_ref->{'AT'}, $TMG_k_t_2_ref->{'CC'}, $TMG_k_t_2_ref->{'CG'}, $TMG_k_t_2_ref->{'CT'}, $TMG_k_t_2_ref->{'GG'}, $TMG_k_t_2_ref->{'GT'}, $TMG_k_t_2_ref->{'TT'}) = (0,0,0,0,0,0,0,0,0,0) ;
			}
			# TMG_k_t_3
			if ($S_k_t_ref->{3}!= 0) {
				($TMG_k_t_3_ref->{'AA'}, $TMG_k_t_3_ref->{'AC'}, $TMG_k_t_3_ref->{'AG'}, $TMG_k_t_3_ref->{'AT'}, $TMG_k_t_3_ref->{'CC'}, $TMG_k_t_3_ref->{'CG'}, $TMG_k_t_3_ref->{'CT'}, $TMG_k_t_3_ref->{'GG'}, $TMG_k_t_3_ref->{'GT'}, $TMG_k_t_3_ref->{'TT'}) = random_multinomial_10($alpha_ref->{'AA'}, $alpha_ref->{'AC'}, $alpha_ref->{'AG'}, $alpha_ref->{'AT'}, $alpha_ref->{'CC'}, $alpha_ref->{'CG'}, $alpha_ref->{'CT'}, $alpha_ref->{'GG'}, $alpha_ref->{'GT'}, $alpha_ref->{'TT'}) ;
			} else {
				($TMG_k_t_3_ref->{'AA'}, $TMG_k_t_3_ref->{'AC'}, $TMG_k_t_3_ref->{'AG'}, $TMG_k_t_3_ref->{'AT'}, $TMG_k_t_3_ref->{'CC'}, $TMG_k_t_3_ref->{'CG'}, $TMG_k_t_3_ref->{'CT'}, $TMG_k_t_3_ref->{'GG'}, $TMG_k_t_3_ref->{'GT'}, $TMG_k_t_3_ref->{'TT'}) = (0,0,0,0,0,0,0,0,0,0) ;
			}
			# TFG_k_t_1
			if ($S_k_t_ref->{1}!= 0) {
				($TFG_k_t_1_ref->{'AA'}, $TFG_k_t_1_ref->{'AC'}, $TFG_k_t_1_ref->{'AG'}, $TFG_k_t_1_ref->{'AT'}, $TFG_k_t_1_ref->{'CC'}, $TFG_k_t_1_ref->{'CG'}, $TFG_k_t_1_ref->{'CT'}, $TFG_k_t_1_ref->{'GG'}, $TFG_k_t_1_ref->{'GT'}, $TFG_k_t_1_ref->{'TT'}) = random_multinomial_10($alpha_ref->{'AA'}, $alpha_ref->{'AC'}, $alpha_ref->{'AG'}, $alpha_ref->{'AT'}, $alpha_ref->{'CC'}, $alpha_ref->{'CG'}, $alpha_ref->{'CT'}, $alpha_ref->{'GG'}, $alpha_ref->{'GT'}, $alpha_ref->{'TT'}) ;
			} else {
				($TFG_k_t_1_ref->{'AA'}, $TFG_k_t_1_ref->{'AC'}, $TFG_k_t_1_ref->{'AG'}, $TFG_k_t_1_ref->{'AT'}, $TFG_k_t_1_ref->{'CC'}, $TFG_k_t_1_ref->{'CG'}, $TFG_k_t_1_ref->{'CT'}, $TFG_k_t_1_ref->{'GG'}, $TFG_k_t_1_ref->{'GT'}, $TFG_k_t_1_ref->{'TT'}) = (0,0,0,0,0,0,0,0,0,0) ;
			}
			# TFG_k_t_2
			if ($S_k_t_ref->{2}!= 0) {
				($TFG_k_t_2_ref->{'AC'}, $TFG_k_t_2_ref->{'CA'}, $TFG_k_t_2_ref->{'AG'}, $TFG_k_t_2_ref->{'GA'}, $TFG_k_t_2_ref->{'AT'}, $TFG_k_t_2_ref->{'TA'}, $TFG_k_t_2_ref->{'CG'}, $TFG_k_t_2_ref->{'GC'}, $TFG_k_t_2_ref->{'CT'}, $TFG_k_t_2_ref->{'TC'}, $TFG_k_t_2_ref->{'GT'}, $TFG_k_t_2_ref->{'TG'}) = random_multinomial_12($beta_2_ref->{'AC'}, $beta_2_ref->{'CA'}, $beta_2_ref->{'AG'}, $beta_2_ref->{'GA'}, $beta_2_ref->{'AT'}, $beta_2_ref->{'TA'}, $beta_2_ref->{'CG'}, $beta_2_ref->{'GC'}, $beta_2_ref->{'CT'}, $beta_2_ref->{'TC'}, $beta_2_ref->{'GT'}, $beta_2_ref->{'TG'}) ;
			} else {
				($TFG_k_t_2_ref->{'AC'}, $TFG_k_t_2_ref->{'CA'}, $TFG_k_t_2_ref->{'AG'}, $TFG_k_t_2_ref->{'GA'}, $TFG_k_t_2_ref->{'AT'}, $TFG_k_t_2_ref->{'TA'}, $TFG_k_t_2_ref->{'CG'}, $TFG_k_t_2_ref->{'GC'}, $TFG_k_t_2_ref->{'CT'}, $TFG_k_t_2_ref->{'TC'}, $TFG_k_t_2_ref->{'GT'}, $TFG_k_t_2_ref->{'TG'}) = (0,0,0,0,0,0,0,0,0,0,0,0) ;
			}
			# TFG_k_t_3
			if ($S_k_t_ref->{3}!= 0) {
				($TFG_k_t_3_ref->{'A'}, $TFG_k_t_3_ref->{'C'}, $TFG_k_t_3_ref->{'G'}, $TFG_k_t_3_ref->{'T'}) = random_multinomial_4($beta_3_ref->{'A'}, $beta_3_ref->{'C'}, $beta_3_ref->{'G'}, $beta_3_ref->{'T'}) ;
			} else {
				($TFG_k_t_3_ref->{'A'}, $TFG_k_t_3_ref->{'C'}, $TFG_k_t_3_ref->{'G'}, $TFG_k_t_3_ref->{'T'}) = (0,0,0,0) ;
			}
			# GE_k_t, YGE_k_t
			for my $m (keys %{$TMG_k_t_1_ref}) {
				for my $n (keys %{$TFG_k_t_1_ref}) {
					for my $i (keys %{$individuals_ref}) {
						if (($TMG_k_t_1_ref->{$m}!= 0)&&($TFG_k_t_1_ref->{$n}!= 0)) {
							$GE_k_t_1_ref->{$n.','.$m}{$i} = random_Bernoulli($Parameters_ref->{'E'}) ;
						}
					}
				}
				for my $n (keys %{$TFG_k_t_2_ref}) {
					for my $i (keys %{$individuals_ref}) {
						if (($TMG_k_t_2_ref->{$m}!= 0)&&($TFG_k_t_2_ref->{$n}!= 0)) {
							$GE_k_t_2_ref->{$n.','.$m}{$i} = random_Bernoulli($Parameters_ref->{'E'}) ;
						}
					}
					for my $i_het (keys %{$het_individuals_ref}) {
						if (($TMG_k_t_2_ref->{$m}!= 0)&&($TFG_k_t_2_ref->{$n}!= 0)) {
							$YGE_k_t_2_ref->{$n.','.$m}{$i_het} = random_Bernoulli($Parameters_ref->{'p'}) ;
						}
					}
				}
				for my $n (keys %{$TFG_k_t_3_ref}) {
					for my $i (keys %{$individuals_ref}) {
						if (($TMG_k_t_3_ref->{$m}!= 0)&&($TFG_k_t_3_ref->{$n}!= 0)) {
							$GE_k_t_3_ref->{$n.','.$m}{$i} = random_Bernoulli($Parameters_ref->{'E'}) ;
						}
					}
				}
			}
		} else {
			# E step
			# Computation of %GE_k_t{i} and %YGE_k_t_2{i_het} and factorise to compute \prod_{i}\mathbb{P}( G_{kt}^{i} | S_{ktj}=1, \TMG_{ktj}^{m}=1, \TFG_{ktj}^{n}=1 )
			($GE_k_t_1_ref, $GE_k_t_2_ref, $GE_k_t_3_ref, $YGE_k_t_2_ref, $factor_ref) = Compute_GE_YGE_and_factorise($Parameters_ref,'new') ;
			# Computation of %S_k_t{j}, %TMG_k_t{m}, %TFG_k_t_j{n}
			($S_k_t_ref, $TMG_k_t_1_ref, $TMG_k_t_2_ref, $TMG_k_t_3_ref, $TFG_k_t_1_ref, $TFG_k_t_2_ref, $TFG_k_t_3_ref) = Compute_S_TMG_TFG_k_t_j($Parameters_ref,'new', $factor_ref) ;
		}

		# Compute Likelihood if parameters estimation over, or if debug mode on
		if ($pi_ref->{2} == 0) {
			$S_k_t_ref->{2} = 0 ;
		}
		if ($pi_ref->{3} == 0) {
			$S_k_t_ref->{3} = 0 ;
		}
		if ((defined($Parameters{'debug'})&&($Parameters_ref->{'run'} > 1))||(($Parameters_ref->{'parameters_estimations_over'} eq "yes")&&(defined($Parameters_ref->{'compute_Likelihood'})))) {
			# Computing Q_k_t and incrementing Q
			$Parameters_ref->{'Q'} += $line_occurence_number * compute_Q_k_t($Parameters_ref, $S_k_t_ref, $TMG_k_t_1_ref, $TMG_k_t_2_ref, $TMG_k_t_3_ref, $TFG_k_t_1_ref, $TFG_k_t_2_ref, $TFG_k_t_3_ref, $GE_k_t_1_ref, $GE_k_t_2_ref, $GE_k_t_3_ref, $YGE_k_t_2_ref) ;
			# Computing H_k_t and incrementing H
			$Parameters_ref->{'H'} += $line_occurence_number * compute_H_k_t($Parameters_ref, $S_k_t_ref, $TMG_k_t_1_ref, $TMG_k_t_2_ref, $TMG_k_t_3_ref, $TFG_k_t_1_ref, $TFG_k_t_2_ref, $TFG_k_t_3_ref, $GE_k_t_1_ref, $GE_k_t_2_ref, $GE_k_t_3_ref, $YGE_k_t_2_ref) ;
		}

		# initializations
		my $sum_YGE_k_t_2_m_n = 0 ;
		my $sum_GE_k_t_1_m_n = 0 ;
		my $sum_GE_k_t_2_m_n = 0 ;
		my $sum_GE_k_t_3_m_n = 0 ;
		my $sum_TMG_TFG_1 = 0 ;
		my $sum_TMG_TFG_2 = 0 ;
		my $sum_TMG_TFG_3 = 0 ;
		my $number_i = ($heterogametic_length + $homogametic_length + 2) ;
		my $number_i_het = ($heterogametic_length + 1) ;
		my $S_k_1_factor = 0 ;
		my $S_k_2_factor = 0 ;
		my $S_k_3_factor = 0 ;
		$Parameters_ref->{'number_individuals_with_aberrant_reads'} = 0 ;

		# S step
		if ((defined($Parameters_ref->{'SEM'}))&&($Parameters_ref->{'parameters_estimations_over'} eq "no")&&($Parameters_ref->{'run'} > 1)&&($Parameters_ref->{'run'} < 11)&&(!defined($Parameters{'skip_optimization'}))) {
			# record old values of variables
			my ($S_k_t_old_ref, $TMG_k_t_1_old_ref, $TMG_k_t_2_old_ref, $TMG_k_t_3_old_ref, $TFG_k_t_1_old_ref, $TFG_k_t_2_old_ref, $TFG_k_t_3_old_ref, $GE_k_t_1_old_ref, $GE_k_t_2_old_ref, $GE_k_t_3_old_ref, $YGE_k_t_2_old_ref) = ($S_k_t_ref, $TMG_k_t_1_ref, $TMG_k_t_2_ref, $TMG_k_t_3_ref, $TFG_k_t_1_ref, $TFG_k_t_2_ref, $TFG_k_t_3_ref, $GE_k_t_1_ref, $GE_k_t_2_ref, $GE_k_t_3_ref, $YGE_k_t_2_ref) ;
			for (my $i = 0; $i < $line_occurence_number; $i++) {
				# S
				($S_k_t_ref->{1}, $S_k_t_ref->{2}, $S_k_t_ref->{3}) = random_multinomial_3($S_k_t_old_ref->{1}, $S_k_t_old_ref->{2}, $S_k_t_old_ref->{3}) ;
				# TMG_k_t_1
				if ($S_k_t_ref->{1}!= 0) {
					($TMG_k_t_1_ref->{'AA'}, $TMG_k_t_1_ref->{'AC'}, $TMG_k_t_1_ref->{'AG'}, $TMG_k_t_1_ref->{'AT'}, $TMG_k_t_1_ref->{'CC'}, $TMG_k_t_1_ref->{'CG'}, $TMG_k_t_1_ref->{'CT'}, $TMG_k_t_1_ref->{'GG'}, $TMG_k_t_1_ref->{'GT'}, $TMG_k_t_1_ref->{'TT'}) = random_multinomial_10($TMG_k_t_1_old_ref->{'AA'}, $TMG_k_t_1_old_ref->{'AC'}, $TMG_k_t_1_old_ref->{'AG'}, $TMG_k_t_1_old_ref->{'AT'}, $TMG_k_t_1_old_ref->{'CC'}, $TMG_k_t_1_old_ref->{'CG'}, $TMG_k_t_1_old_ref->{'CT'}, $TMG_k_t_1_old_ref->{'GG'}, $TMG_k_t_1_old_ref->{'GT'}, $TMG_k_t_1_old_ref->{'TT'}) ;
				} else {
					($TMG_k_t_1_ref->{'AA'}, $TMG_k_t_1_ref->{'AC'}, $TMG_k_t_1_ref->{'AG'}, $TMG_k_t_1_ref->{'AT'}, $TMG_k_t_1_ref->{'CC'}, $TMG_k_t_1_ref->{'CG'}, $TMG_k_t_1_ref->{'CT'}, $TMG_k_t_1_ref->{'GG'}, $TMG_k_t_1_ref->{'GT'}, $TMG_k_t_1_ref->{'TT'}) = (0,0,0,0,0,0,0,0,0,0) ;
				}
				# TMG_k_t_2
				if ($S_k_t_ref->{2}!= 0) {
					($TMG_k_t_2_ref->{'AA'}, $TMG_k_t_2_ref->{'AC'}, $TMG_k_t_2_ref->{'AG'}, $TMG_k_t_2_ref->{'AT'}, $TMG_k_t_2_ref->{'CC'}, $TMG_k_t_2_ref->{'CG'}, $TMG_k_t_2_ref->{'CT'}, $TMG_k_t_2_ref->{'GG'}, $TMG_k_t_2_ref->{'GT'}, $TMG_k_t_2_ref->{'TT'}) = random_multinomial_10($TMG_k_t_2_old_ref->{'AA'}, $TMG_k_t_2_old_ref->{'AC'}, $TMG_k_t_2_old_ref->{'AG'}, $TMG_k_t_2_old_ref->{'AT'}, $TMG_k_t_2_old_ref->{'CC'}, $TMG_k_t_2_old_ref->{'CG'}, $TMG_k_t_2_old_ref->{'CT'}, $TMG_k_t_2_old_ref->{'GG'}, $TMG_k_t_2_old_ref->{'GT'}, $TMG_k_t_2_old_ref->{'TT'}) ;
				} else {
					($TMG_k_t_2_ref->{'AA'}, $TMG_k_t_2_ref->{'AC'}, $TMG_k_t_2_ref->{'AG'}, $TMG_k_t_2_ref->{'AT'}, $TMG_k_t_2_ref->{'CC'}, $TMG_k_t_2_ref->{'CG'}, $TMG_k_t_2_ref->{'CT'}, $TMG_k_t_2_ref->{'GG'}, $TMG_k_t_2_ref->{'GT'}, $TMG_k_t_2_ref->{'TT'}) = (0,0,0,0,0,0,0,0,0,0) ;
				}
				# TMG_k_t_3
				if ($S_k_t_ref->{3}!= 0) {
					($TMG_k_t_3_ref->{'AA'}, $TMG_k_t_3_ref->{'AC'}, $TMG_k_t_3_ref->{'AG'}, $TMG_k_t_3_ref->{'AT'}, $TMG_k_t_3_ref->{'CC'}, $TMG_k_t_3_ref->{'CG'}, $TMG_k_t_3_ref->{'CT'}, $TMG_k_t_3_ref->{'GG'}, $TMG_k_t_3_ref->{'GT'}, $TMG_k_t_3_ref->{'TT'}) = random_multinomial_10($TMG_k_t_3_old_ref->{'AA'}, $TMG_k_t_3_old_ref->{'AC'}, $TMG_k_t_3_old_ref->{'AG'}, $TMG_k_t_3_old_ref->{'AT'}, $TMG_k_t_3_old_ref->{'CC'}, $TMG_k_t_3_old_ref->{'CG'}, $TMG_k_t_3_old_ref->{'CT'}, $TMG_k_t_3_old_ref->{'GG'}, $TMG_k_t_3_old_ref->{'GT'}, $TMG_k_t_3_old_ref->{'TT'}) ;
				} else {
					($TMG_k_t_3_ref->{'AA'}, $TMG_k_t_3_ref->{'AC'}, $TMG_k_t_3_ref->{'AG'}, $TMG_k_t_3_ref->{'AT'}, $TMG_k_t_3_ref->{'CC'}, $TMG_k_t_3_ref->{'CG'}, $TMG_k_t_3_ref->{'CT'}, $TMG_k_t_3_ref->{'GG'}, $TMG_k_t_3_ref->{'GT'}, $TMG_k_t_3_ref->{'TT'}) = (0,0,0,0,0,0,0,0,0,0) ;
				}
				# TFG_k_t_1
				if ($S_k_t_ref->{1}!= 0) {
					($TFG_k_t_1_ref->{'AA'}, $TFG_k_t_1_ref->{'AC'}, $TFG_k_t_1_ref->{'AG'}, $TFG_k_t_1_ref->{'AT'}, $TFG_k_t_1_ref->{'CC'}, $TFG_k_t_1_ref->{'CG'}, $TFG_k_t_1_ref->{'CT'}, $TFG_k_t_1_ref->{'GG'}, $TFG_k_t_1_ref->{'GT'}, $TFG_k_t_1_ref->{'TT'}) = random_multinomial_10($TFG_k_t_1_old_ref->{'AA'}, $TFG_k_t_1_old_ref->{'AC'}, $TFG_k_t_1_old_ref->{'AG'}, $TFG_k_t_1_old_ref->{'AT'}, $TFG_k_t_1_old_ref->{'CC'}, $TFG_k_t_1_old_ref->{'CG'}, $TFG_k_t_1_old_ref->{'CT'}, $TFG_k_t_1_old_ref->{'GG'}, $TFG_k_t_1_old_ref->{'GT'}, $TFG_k_t_1_old_ref->{'TT'}) ;
				} else {
					($TFG_k_t_1_ref->{'AA'}, $TFG_k_t_1_ref->{'AC'}, $TFG_k_t_1_ref->{'AG'}, $TFG_k_t_1_ref->{'AT'}, $TFG_k_t_1_ref->{'CC'}, $TFG_k_t_1_ref->{'CG'}, $TFG_k_t_1_ref->{'CT'}, $TFG_k_t_1_ref->{'GG'}, $TFG_k_t_1_ref->{'GT'}, $TFG_k_t_1_ref->{'TT'}) = (0,0,0,0,0,0,0,0,0,0) ;
				}
				# TFG_k_t_2
				if ($S_k_t_ref->{2}!= 0) {
					($TFG_k_t_2_ref->{'AC'}, $TFG_k_t_2_ref->{'CA'}, $TFG_k_t_2_ref->{'AG'}, $TFG_k_t_2_ref->{'GA'}, $TFG_k_t_2_ref->{'AT'}, $TFG_k_t_2_ref->{'TA'}, $TFG_k_t_2_ref->{'CG'}, $TFG_k_t_2_ref->{'GC'}, $TFG_k_t_2_ref->{'CT'}, $TFG_k_t_2_ref->{'TC'}, $TFG_k_t_2_ref->{'GT'}, $TFG_k_t_2_ref->{'TG'}) = random_multinomial_12($TFG_k_t_2_old_ref->{'AC'}, $TFG_k_t_2_old_ref->{'CA'}, $TFG_k_t_2_old_ref->{'AG'}, $TFG_k_t_2_old_ref->{'GA'}, $TFG_k_t_2_old_ref->{'AT'}, $TFG_k_t_2_old_ref->{'TA'}, $TFG_k_t_2_old_ref->{'CG'}, $TFG_k_t_2_old_ref->{'GC'}, $TFG_k_t_2_old_ref->{'CT'}, $TFG_k_t_2_old_ref->{'TC'}, $TFG_k_t_2_old_ref->{'GT'}, $TFG_k_t_2_old_ref->{'TG'}) ;
				} else {
					($TFG_k_t_2_ref->{'AC'}, $TFG_k_t_2_ref->{'CA'}, $TFG_k_t_2_ref->{'AG'}, $TFG_k_t_2_ref->{'GA'}, $TFG_k_t_2_ref->{'AT'}, $TFG_k_t_2_ref->{'TA'}, $TFG_k_t_2_ref->{'CG'}, $TFG_k_t_2_ref->{'GC'}, $TFG_k_t_2_ref->{'CT'}, $TFG_k_t_2_ref->{'TC'}, $TFG_k_t_2_ref->{'GT'}, $TFG_k_t_2_ref->{'TG'}) = (0,0,0,0,0,0,0,0,0,0,0,0) ;
				}
				# TFG_k_t_3
				if ($S_k_t_ref->{3}!= 0) {
					($TFG_k_t_3_ref->{'A'}, $TFG_k_t_3_ref->{'C'}, $TFG_k_t_3_ref->{'G'}, $TFG_k_t_3_ref->{'T'}) = random_multinomial_4($TFG_k_t_3_old_ref->{'A'}, $TFG_k_t_3_old_ref->{'C'}, $TFG_k_t_3_old_ref->{'G'}, $TFG_k_t_3_old_ref->{'T'}) ;
				} else {
					($TFG_k_t_3_ref->{'A'}, $TFG_k_t_3_ref->{'C'}, $TFG_k_t_3_ref->{'G'}, $TFG_k_t_3_ref->{'T'}) = (0,0,0,0) ;
				}
				# GE_k_t, YGE_k_t
				for my $m (keys %{$TMG_k_t_1_ref}) {
					for my $n (keys %{$TFG_k_t_1_ref}) {
						for my $i (keys %{$individuals_ref}) {
							if (($TMG_k_t_1_ref->{$m}!= 0)&&($TFG_k_t_1_ref->{$n}!= 0)) {
								$GE_k_t_1_ref->{$n.','.$m}{$i} = random_Bernoulli($GE_k_t_1_old_ref->{$n.','.$m}{$i}) ;
							} else {
								$GE_k_t_1_ref->{$n.','.$m}{$i} = 0 ;
							}
						}
					}
					for my $n (keys %{$TFG_k_t_2_ref}) {
						for my $i (keys %{$individuals_ref}) {
							if (($TMG_k_t_2_ref->{$m}!= 0)&&($TFG_k_t_2_ref->{$n}!= 0)) {
								$GE_k_t_2_ref->{$n.','.$m}{$i} = random_Bernoulli($GE_k_t_2_old_ref->{$n.','.$m}{$i}) ;
							} else {
								$GE_k_t_2_ref->{$n.','.$m}{$i} = 0 ;
							}
						}
						for my $i_het (keys %{$het_individuals_ref}) {
							if (($TMG_k_t_2_ref->{$m}!= 0)&&($TFG_k_t_2_ref->{$n}!= 0)) {
								$YGE_k_t_2_ref->{$n.','.$m}{$i_het} = random_Bernoulli($YGE_k_t_2_old_ref->{$n.','.$m}{$i_het}) ;
							} else {
								$YGE_k_t_2_ref->{$n.','.$m}{$i_het} = 0 ;
							}
						}
					}
					for my $n (keys %{$TFG_k_t_3_ref}) {
						for my $i (keys %{$individuals_ref}) {
							if (($TMG_k_t_3_ref->{$m}!= 0)&&($TFG_k_t_3_ref->{$n}!= 0)) {
								$GE_k_t_3_ref->{$n.','.$m}{$i} = random_Bernoulli($GE_k_t_3_old_ref->{$n.','.$m}{$i}) ;
							} else {
								$GE_k_t_3_ref->{$n.','.$m}{$i} = 0 ;
							}
						}
					}
				}

				# incrementing sums for new parameters value computation
				$Parameters_ref->{'sum_k_t'} += 1 ;
				$sum_k_t_S_k_t_ref->{1} += $S_k_t_ref->{1} ;
				$sum_k_t_S_k_t_ref->{2} += $S_k_t_ref->{2} ;
				$sum_k_t_S_k_t_ref->{3} += $S_k_t_ref->{3} ;
				for my $m (keys %{$TMG_k_t_1_ref}) {
					if ($S_k_t_ref->{1}!= 0) {
						$alpha_new_ref->{$m} += $S_k_t_ref->{1} * $TMG_k_t_1_ref->{$m} ;
						$alpha_new_ref->{$m} += $S_k_t_ref->{1} * $TFG_k_t_1_ref->{$m} ;
					}
					if ($S_k_t_ref->{2}!= 0) {
						$alpha_new_ref->{$m} += $S_k_t_ref->{2} * $TMG_k_t_2_ref->{$m} ;
					}
					if ($S_k_t_ref->{3}!= 0) {
						$alpha_new_ref->{$m} += $S_k_t_ref->{3} * $TMG_k_t_3_ref->{$m} ;
					}
				}
				if ($S_k_t_ref->{2}!= 0) {
					for my $n (keys %{$TFG_k_t_2_ref}) {
						$beta_2_new_ref->{$n} += $S_k_t_ref->{2} * $TFG_k_t_2_ref->{$n} ;
					}
				}
				if ($S_k_t_ref->{3}!= 0) {
					for my $n (keys %{$TFG_k_t_3_ref}) {
						$beta_3_new_ref->{$n} += $S_k_t_ref->{3} * $TFG_k_t_3_ref->{$n} ;
					}
				}
				for my $m (keys %{$TMG_k_t_1_ref}) {
					if ($S_k_t_ref->{1}!= 0) {
						for my $n (keys %{$TFG_k_t_1_ref}) {
							if (($TMG_k_t_1_ref->{$m}!= 0)&&($TFG_k_t_1_ref->{$n}!= 0)) {
								my $sum_GE_k_t_1 = 0 ;
								my $one_minus_GE_k_t_1_sum = 0 ;
								for my $i (keys %{$individuals_ref}) {
									$sum_GE_k_t_1 += $GE_k_t_1_ref->{$n.','.$m}{$i} ;
									$one_minus_GE_k_t_1_sum += 1 - $GE_k_t_1_ref->{$n.','.$m}{$i} ;
								}
								$sum_GE_k_t_1_m_n += $sum_GE_k_t_1 * $TMG_k_t_1_ref->{$m} * $TFG_k_t_1_ref->{$n} ;
								$sum_TMG_TFG_1 += $TMG_k_t_1_ref->{$m} * $TFG_k_t_1_ref->{$n} ;
								$S_k_1_factor += $one_minus_GE_k_t_1_sum * $TMG_k_t_1_ref->{$m} * $TFG_k_t_1_ref->{$n} ;
							}
						}
					}
					if ($S_k_t_ref->{2}!= 0) {
						for my $n (keys %{$TFG_k_t_2_ref}) {
							if (($TMG_k_t_2_ref->{$m}!= 0)&&($TFG_k_t_2_ref->{$n}!= 0)) {
								my $sum_GE_k_t_2 = 0 ;
								my $one_minus_GE_k_t_2_sum = 0 ;
								for my $i (keys %{$individuals_ref}) {
									$sum_GE_k_t_2 += $GE_k_t_2_ref->{$n.','.$m}{$i} ;
									$one_minus_GE_k_t_2_sum += 1 - $GE_k_t_2_ref->{$n.','.$m}{$i} ;
								}
								$sum_GE_k_t_2_m_n += $sum_GE_k_t_2 * $TMG_k_t_2_ref->{$m} * $TFG_k_t_2_ref->{$n} ;
								my $sum_YGE_k_t_2 = 0 ;
								my $one_minus_YGE_k_t_2_sum = 0 ;
								for my $i_het (keys %{$het_individuals_ref}) {
									$sum_YGE_k_t_2 += $YGE_k_t_2_ref->{$n.','.$m}{$i_het} ;
									$one_minus_YGE_k_t_2_sum += 1 - $YGE_k_t_2_ref->{$n.','.$m}{$i_het} ;
								}
								$sum_YGE_k_t_2_m_n += $sum_YGE_k_t_2 * $TMG_k_t_2_ref->{$m} * $TFG_k_t_2_ref->{$n} ;
								$sum_TMG_TFG_2 += $TMG_k_t_2_ref->{$m} * $TFG_k_t_2_ref->{$n} ;
								$S_k_2_factor += $one_minus_GE_k_t_2_sum * $one_minus_YGE_k_t_2_sum * $TMG_k_t_2_ref->{$m} * $TFG_k_t_2_ref->{$n} ;
							}
						}
					}
					if ($S_k_t_ref->{3}!= 0) {
						for my $n (keys %{$TFG_k_t_3_ref}) {
							if (($TMG_k_t_3_ref->{$m}!= 0)&&($TFG_k_t_3_ref->{$n}!= 0)) {
								my $sum_GE_k_t_3 = 0 ;
								my $one_minus_GE_k_t_3_sum = 0 ;
								for my $i (keys %{$individuals_ref}) {
									$sum_GE_k_t_3 += $GE_k_t_3_ref->{$n.','.$m}{$i} ;
									$one_minus_GE_k_t_3_sum += 1 - $GE_k_t_3_ref->{$n.','.$m}{$i} ;
								}
								$sum_GE_k_t_3_m_n += $sum_GE_k_t_3 * $TMG_k_t_3_ref->{$m} * $TFG_k_t_3_ref->{$n} ;
								$sum_TMG_TFG_3 += $TMG_k_t_3_ref->{$m} * $TFG_k_t_3_ref->{$n} ;
								$S_k_3_factor += $one_minus_GE_k_t_3_sum * $TMG_k_t_3_ref->{$m} * $TFG_k_t_3_ref->{$n} ;
							}
						}
					}
				}
				if ($S_k_t_ref->{1}!= 0) {
					$Parameters_ref->{'E_new_numerator'} += $S_k_t_ref->{1} * $sum_GE_k_t_1_m_n ;
					$Parameters_ref->{'E_new_denominator'} += $S_k_t_ref->{1} * $sum_TMG_TFG_1 * $number_i ;
				}
				if ($S_k_t_ref->{2}!= 0) {
					$Parameters_ref->{'E_new_numerator'} += $S_k_t_ref->{2} * $sum_GE_k_t_2_m_n ;
					$Parameters_ref->{'E_new_denominator'} += $S_k_t_ref->{2} * $sum_TMG_TFG_2 * $number_i ;
					$Parameters_ref->{'p_new_numerator'} += $S_k_t_ref->{2} * $sum_YGE_k_t_2_m_n ;
					$Parameters_ref->{'p_new_denominator'} += $S_k_t_ref->{2} * $sum_TMG_TFG_2 * $number_i_het ;
				}
				if ($S_k_t_ref->{3}!= 0) {
					$Parameters_ref->{'E_new_numerator'} += $S_k_t_ref->{3} * $sum_GE_k_t_3_m_n ;
					$Parameters_ref->{'E_new_denominator'} += $S_k_t_ref->{3} * $sum_TMG_TFG_3 * $number_i ;
				}
			}
		} else {
			# no S step
			# incrementing sums for new parameters value computation
			$Parameters_ref->{'sum_k_t'} += $line_occurence_number ;
			$sum_k_t_S_k_t_ref->{1} += $line_occurence_number * $S_k_t_ref->{1} ;
			$sum_k_t_S_k_t_ref->{2} += $line_occurence_number * $S_k_t_ref->{2} ;
			$sum_k_t_S_k_t_ref->{3} += $line_occurence_number * $S_k_t_ref->{3} ;
			for my $m (keys %{$TMG_k_t_1_ref}) {
				if ($S_k_t_ref->{1}!= 0) {
					$alpha_new_ref->{$m} += $line_occurence_number * $S_k_t_ref->{1} * $TMG_k_t_1_ref->{$m} ;
					$alpha_new_ref->{$m} += $line_occurence_number * $S_k_t_ref->{1} * $TFG_k_t_1_ref->{$m} ;
				}
				if ($S_k_t_ref->{2}!= 0) {
					$alpha_new_ref->{$m} += $line_occurence_number * $S_k_t_ref->{2} * $TMG_k_t_2_ref->{$m} ;
				}
				if ($S_k_t_ref->{3}!= 0) {
					$alpha_new_ref->{$m} += $line_occurence_number * $S_k_t_ref->{3} * $TMG_k_t_3_ref->{$m} ;
				}
			}
			if ($S_k_t_ref->{2}!= 0) {
				for my $n (keys %{$TFG_k_t_2_ref}) {
					$beta_2_new_ref->{$n} += $line_occurence_number * $S_k_t_ref->{2} * $TFG_k_t_2_ref->{$n} ;
				}
			}
			if ($S_k_t_ref->{3}!= 0) {
				for my $n (keys %{$TFG_k_t_3_ref}) {
					$beta_3_new_ref->{$n} += $line_occurence_number * $S_k_t_ref->{3} * $TFG_k_t_3_ref->{$n} ;
				}
			}
			for my $m (keys %{$TMG_k_t_1_ref}) {
				if ($S_k_t_ref->{1}!= 0) {
					for my $n (keys %{$TFG_k_t_1_ref}) {
						if (($TMG_k_t_1_ref->{$m}!= 0)&&($TFG_k_t_1_ref->{$n}!= 0)) {
							my $sum_GE_k_t_1 = 0 ;
							my $one_minus_GE_k_t_1_sum = 0 ;
							for my $i (keys %{$individuals_ref}) {
								$sum_GE_k_t_1 += $GE_k_t_1_ref->{$n.','.$m}{$i} ;
								$one_minus_GE_k_t_1_sum += 1 - $GE_k_t_1_ref->{$n.','.$m}{$i} ;
							}
							$sum_GE_k_t_1_m_n += $sum_GE_k_t_1 * $TMG_k_t_1_ref->{$m} * $TFG_k_t_1_ref->{$n} ;
							$sum_TMG_TFG_1 += $TMG_k_t_1_ref->{$m} * $TFG_k_t_1_ref->{$n} ;
							$S_k_1_factor += $one_minus_GE_k_t_1_sum * $TMG_k_t_1_ref->{$m} * $TFG_k_t_1_ref->{$n} ;
						}
					}
				}
				if ($S_k_t_ref->{2}!= 0) {
					for my $n (keys %{$TFG_k_t_2_ref}) {
						if (($TMG_k_t_2_ref->{$m}!= 0)&&($TFG_k_t_2_ref->{$n}!= 0)) {
							my $sum_GE_k_t_2 = 0 ;
							my $one_minus_GE_k_t_2_sum = 0 ;
							for my $i (keys %{$individuals_ref}) {
								$sum_GE_k_t_2 += $GE_k_t_2_ref->{$n.','.$m}{$i} ;
								$one_minus_GE_k_t_2_sum += 1 - $GE_k_t_2_ref->{$n.','.$m}{$i} ;
							}
							$sum_GE_k_t_2_m_n += $sum_GE_k_t_2 * $TMG_k_t_2_ref->{$m} * $TFG_k_t_2_ref->{$n} ;
							my $sum_YGE_k_t_2 = 0 ;
							my $one_minus_YGE_k_t_2_sum = 0 ;
							for my $i_het (keys %{$het_individuals_ref}) {
								$sum_YGE_k_t_2 += $YGE_k_t_2_ref->{$n.','.$m}{$i_het} ;
								$one_minus_YGE_k_t_2_sum += 1 - $YGE_k_t_2_ref->{$n.','.$m}{$i_het} ;
							}
							$sum_YGE_k_t_2_m_n += $sum_YGE_k_t_2 * $TMG_k_t_2_ref->{$m} * $TFG_k_t_2_ref->{$n} ;
							$sum_TMG_TFG_2 += $TMG_k_t_2_ref->{$m} * $TFG_k_t_2_ref->{$n} ;
							$S_k_2_factor += $one_minus_GE_k_t_2_sum * $one_minus_YGE_k_t_2_sum * $TMG_k_t_2_ref->{$m} * $TFG_k_t_2_ref->{$n} ;
						}
					}
				}
				if ($S_k_t_ref->{3}!= 0) {
					for my $n (keys %{$TFG_k_t_3_ref}) {
						if (($TMG_k_t_3_ref->{$m}!= 0)&&($TFG_k_t_3_ref->{$n}!= 0)) {
							my $sum_GE_k_t_3 = 0 ;
							my $one_minus_GE_k_t_3_sum = 0 ;
							for my $i (keys %{$individuals_ref}) {
								$sum_GE_k_t_3 += $GE_k_t_3_ref->{$n.','.$m}{$i} ;
								$one_minus_GE_k_t_3_sum += 1 - $GE_k_t_3_ref->{$n.','.$m}{$i} ;
							}
							$sum_GE_k_t_3_m_n += $sum_GE_k_t_3 * $TMG_k_t_3_ref->{$m} * $TFG_k_t_3_ref->{$n} ;
							$sum_TMG_TFG_3 += $TMG_k_t_3_ref->{$m} * $TFG_k_t_3_ref->{$n} ;
							$S_k_3_factor += $one_minus_GE_k_t_3_sum * $TMG_k_t_3_ref->{$m} * $TFG_k_t_3_ref->{$n} ;
						}
					}
				}
			}
			if ($S_k_t_ref->{1}!= 0) {
				$Parameters_ref->{'E_new_numerator'} += $line_occurence_number * $S_k_t_ref->{1} * $sum_GE_k_t_1_m_n ;
				$Parameters_ref->{'E_new_denominator'} += $line_occurence_number * $S_k_t_ref->{1} * $sum_TMG_TFG_1 * $number_i ;
			}
			if ($S_k_t_ref->{2}!= 0) {
				$Parameters_ref->{'E_new_numerator'} += $line_occurence_number * $S_k_t_ref->{2} * $sum_GE_k_t_2_m_n ;
				$Parameters_ref->{'E_new_denominator'} += $line_occurence_number * $S_k_t_ref->{2} * $sum_TMG_TFG_2 * $number_i ;
				$Parameters_ref->{'p_new_numerator'} += $line_occurence_number * $S_k_t_ref->{2} * $sum_YGE_k_t_2_m_n ;
				$Parameters_ref->{'p_new_denominator'} += $line_occurence_number * $S_k_t_ref->{2} * $sum_TMG_TFG_2 * $number_i_het ;
			}
			if ($S_k_t_ref->{3}!= 0) {
				$Parameters_ref->{'E_new_numerator'} += $line_occurence_number * $S_k_t_ref->{3} * $sum_GE_k_t_3_m_n ;
				$Parameters_ref->{'E_new_denominator'} += $line_occurence_number * $S_k_t_ref->{3} * $sum_TMG_TFG_3 * $number_i ;
			}
		}

		# retrieve the most probable parents genotypes and segregation type
		my $higher_TMG_k_t_1 = "" ;
		my $higher_TMG_k_t_1_proba = 0 ;
		for my $m (keys %{$TMG_k_t_1_ref}) {
			if ($TMG_k_t_1_ref->{$m} > $higher_TMG_k_t_1_proba) {
				$higher_TMG_k_t_1 = $m ;
				$higher_TMG_k_t_1_proba = $TMG_k_t_1_ref->{$m} ;
			}
		}
		my $higher_TMG_k_t_2 = "" ;
		my $higher_TMG_k_t_2_proba = 0 ;
		for my $m (keys %{$TMG_k_t_2_ref}) {
			if ($TMG_k_t_2_ref->{$m} > $higher_TMG_k_t_2_proba) {
				$higher_TMG_k_t_2 = $m ;
				$higher_TMG_k_t_2_proba = $TMG_k_t_2_ref->{$m} ;
			}
		}
		my $higher_TMG_k_t_3 = "" ;
		my $higher_TMG_k_t_3_proba = 0 ;
		for my $m (keys %{$TMG_k_t_3_ref}) {
			if ($TMG_k_t_3_ref->{$m} > $higher_TMG_k_t_3_proba) {
				$higher_TMG_k_t_3 = $m ;
				$higher_TMG_k_t_3_proba = $TMG_k_t_3_ref->{$m} ;
			}
		}
		my $higher_TFG_k_t_1 = "" ;
		my $higher_TFG_k_t_1_proba = 0 ;
		for my $n (keys %{$TFG_k_t_1_ref}) {
			if ($TFG_k_t_1_ref->{$n} > $higher_TFG_k_t_1_proba) {
				$higher_TFG_k_t_1 = $n ;
				$higher_TFG_k_t_1_proba = $TFG_k_t_1_ref->{$n} ;
			}
		}
		my $higher_TFG_k_t_2 = "" ;
		my $higher_TFG_k_t_2_proba = 0 ;
		for my $n (keys %{$TFG_k_t_2_ref}) {
			if ($TFG_k_t_2_ref->{$n} > $higher_TFG_k_t_2_proba) {
				$higher_TFG_k_t_2 = $n ;
				$higher_TFG_k_t_2_proba = $TFG_k_t_2_ref->{$n} ;
			}
		}
		my $higher_TFG_k_t_3 = "" ;
		my $higher_TFG_k_t_3_proba = 0 ;
		for my $n (keys %{$TFG_k_t_3_ref}) {
			if ($TFG_k_t_3_ref->{$n} > $higher_TFG_k_t_3_proba) {
				$higher_TFG_k_t_3 = $n ;
				$higher_TFG_k_t_3_proba = $TFG_k_t_3_ref->{$n} ;
			}
		}
		my $higher_S_k_t = "" ;
		my $higher_S_k_t_proba = 0 ;
		for my $j (keys %{$S_k_t_ref}) {
			if ($S_k_t_ref->{$j} > $higher_S_k_t_proba) {
				$higher_S_k_t = $j ;
				$higher_S_k_t_proba = $S_k_t_ref->{$j} ;
			}
		}

		# add counts of most probable parent genotypes and segregation type
		if ($higher_S_k_t  == 1) {
			$TMG_1_counts_ref->{$higher_TMG_k_t_1} += $line_occurence_number ;
			$TFG_1_counts_ref->{$higher_TFG_k_t_1} += $line_occurence_number ;
		} elsif ($higher_S_k_t  == 2) {
			$TMG_2_counts_ref->{$higher_TMG_k_t_2} += $line_occurence_number ;
			$TFG_2_counts_ref->{$higher_TFG_k_t_2} += $line_occurence_number ;
		} else {
			$TMG_3_counts_ref->{$higher_TMG_k_t_3} += $line_occurence_number ;
			$TFG_3_counts_ref->{$higher_TFG_k_t_3} += $line_occurence_number ;
		}
		$S_counts_ref->{$higher_S_k_t} += $line_occurence_number ;

		# retrieve the most probable SNP types for each j
		my $higher_1_type = "" ;
		if ($S_k_t_ref->{1} > 0) {
			$higher_1_type = $lambda_1_ref->{$higher_TFG_k_t_1.','.$higher_TMG_k_t_1}{'SNP_type'} ;
		}
		my $higher_2_type = "" ;
		if ($S_k_t_ref->{2}!= 0) {
			$higher_2_type = $lambda_2_ref->{$higher_TFG_k_t_2.','.$higher_TMG_k_t_2}{'SNP_type'} ;
		}
		my $higher_3_type = "" ;
		if ($S_k_t_ref->{3}!= 0) {
			$higher_3_type = $lambda_3_ref->{$higher_TFG_k_t_3.','.$higher_TMG_k_t_3}{'SNP_type'} ;
		}

		# retrieve the number of genotyping errors that were done for the most probable parent genotypes
		my $number_error_1 = 0 ;
		my $number_error_2 = 0 ;
		my $number_error_3 = 0 ;
		my $number_Y_error = 0 ;
		for my $i (keys %{$individuals_ref}) {
			if ($S_k_t_ref->{1}!= 0) {
				if ($GE_k_t_1_ref->{$higher_TFG_k_t_1.','.$higher_TMG_k_t_1}{$i} > 0.5) {
					$number_error_1++ ;
				}
			}
			if ($S_k_t_ref->{2}!= 0) {
				if ($GE_k_t_2_ref->{$higher_TFG_k_t_2.','.$higher_TMG_k_t_2}{$i} > 0.5) {
					$number_error_2++ ;
				}
			}
			if ($S_k_t_ref->{3}!= 0) {
				if ($GE_k_t_3_ref->{$higher_TFG_k_t_3.','.$higher_TMG_k_t_3}{$i} > 0.5) {
					$number_error_3++ ;
				}
			}
		}
		for my $i_het (keys %{$het_individuals_ref}) {
			if ($S_k_t_ref->{2}!= 0) {
				if ($YGE_k_t_2_ref->{$higher_TFG_k_t_2.','.$higher_TMG_k_t_2}{$i_het} > 0.5) {
					$number_Y_error++ ;
				}
			}
		}

		# print SNP details if requested and look at Y expression in homogametic individuals if possible
		if ($Parameters_ref->{'parameters_estimations_over'} eq "yes") {
			my $homogametic_ref = $Parameters{'homogametic'} ;
			my @homogametic = @{$homogametic_ref} ;
			my $heterogametic_ref = $Parameters{'heterogametic'} ;
			my @heterogametic = @{$heterogametic_ref} ;
			if (defined($Parameters_ref->{'ALR_file'})) {
				# initialise expression levels
				my %heterogametic_expression; # keys = heterogametic individuals' names, values = reads counts
				$Parameters{'heterogametic_expression'} = \%heterogametic_expression ;
				my %homogametic_expression; # keys = homogametic individuals' names, values = reads counts
				$Parameters{'homogametic_expression'} = \%homogametic_expression ;
				$Parameters_ref->{'homogametic_parent_expression'} = "" ;
				$Parameters_ref->{'heterogametic_parent_expression'} = "" ;
				# retrieve expression levels
				retrieve_reads2snp_expression_levels($Parameters_ref) ;
				# count number of homogametic individuals with Y expression
				$Parameters_ref->{'number_homogametic_with_Y_expr'} = 0 ;
				if ((($S_k_t_ref->{2} > $S_k_t_ref->{1})&&($S_k_t_ref->{2} > $S_k_t_ref->{3}))&&(($higher_2_type eq 'XY')||($higher_2_type eq 'XXY')||($higher_2_type eq 'XXXY'))) {
					my $Y = substr($higher_TFG_k_t_2, 1, 1) ;
					if ($Parameters_ref->{'homogametic_parent_expression'} ne "NN") {
						my $expr = retrieve_base_expression($Y, $Parameters_ref->{'homogametic_parent_expression'}) ;
						my $total_expr = retrieve_total_expression($Parameters_ref->{'homogametic_parent_expression'}) ;
						if (($total_expr != 0)&&(($expr/$total_expr) > 0.02)) {
							$Parameters_ref->{'number_homogametic_with_Y_expr'} ++ ;
						}
					}
					foreach my $individual (@homogametic) {
						if ($homogametic_expression{$individual} ne "NN") {
							my $expr = retrieve_base_expression($Y, $homogametic_expression{$individual}) ;
							my $total_expr = retrieve_total_expression($homogametic_expression{$individual}) ;
							if (($total_expr != 0)&&(($expr/$total_expr) > 0.02)) {
								$Parameters_ref->{'number_homogametic_with_Y_expr'} ++ ;
							}
						}
					}
				}
				# count number of heterogametic individuals with paternal X or maternal Z expression
				$Parameters_ref->{'number_heterogametic_with_paternalX_expr'} = 0 ;
				if ((($S_k_t_ref->{2} > $S_k_t_ref->{1})&&($S_k_t_ref->{2} > $S_k_t_ref->{3}))&&((($higher_2_type eq 'XX')||($higher_2_type eq 'XXX')||($higher_2_type eq 'XXY')||($higher_2_type eq 'XXXY')))) {
					my $paternalX = substr($higher_TFG_k_t_2, 0, 1) ;
					my $maternalX1 = substr($higher_TMG_k_t_2, 0, 1) ;
					my $maternalX2 = substr($higher_TMG_k_t_2, 1, 1) ;
					if (($paternalX ne $maternalX1)&&($paternalX ne $maternalX2)) {
						# Heterogametic progeny should not have any reads bearing the paternal X allele
						foreach my $individual (@heterogametic) {
							if ($heterogametic_expression{$individual} ne "NN") {
								my $expr = retrieve_base_expression($paternalX, $heterogametic_expression{$individual}) ;
								my $total_expr = retrieve_total_expression($heterogametic_expression{$individual}) ;
								if (($total_expr != 0)&&(($expr/$total_expr) > 0.02)) {
									$Parameters_ref->{'number_heterogametic_with_paternalX_expr'} ++ ;
								}
							}
						}
					}
				}
				if (($S_k_t_ref->{3} > $S_k_t_ref->{1})&&($S_k_t_ref->{3} > $S_k_t_ref->{2})) {
					my $paternalX = substr($higher_TFG_k_t_3, 0, 1) ;
					my $maternalX1 = substr($higher_TMG_k_t_3, 0, 1) ;
					my $maternalX2 = substr($higher_TMG_k_t_3, 1, 1) ;
					if (($paternalX ne $maternalX1)&&($paternalX ne $maternalX2)) {
						# Heterogametic progeny should not have any reads bearing the paternal X allele
						foreach my $individual (@heterogametic) {
							if ($heterogametic_expression{$individual} ne "NN") {
								my $expr = retrieve_base_expression($paternalX, $heterogametic_expression{$individual}) ;
								my $total_expr = retrieve_total_expression($heterogametic_expression{$individual}) ;
								if (($total_expr != 0)&&(($expr/$total_expr) > 0.02)) {
									$Parameters_ref->{'number_heterogametic_with_paternalX_expr'} ++ ;
								}
							}
						}
					}
				}
				$Parameters_ref->{'number_individuals_with_aberrant_reads'} = $Parameters_ref->{'number_homogametic_with_Y_expr'} + $Parameters_ref->{'number_heterogametic_with_paternalX_expr'} ;
				# print raw SNP output
				if (defined($Parameters_ref->{'detail'})) {
					open(DETAIL, ">>".$Parameters_ref->{'detail_file_name'}) ;
					print(DETAIL $Parameters_ref->{'contig_name'}."\t".$Parameters_ref->{'current_contig_position'}."\t".$S_k_t_ref->{1}."\t".$S_k_t_ref->{2}."\t".$S_k_t_ref->{3}."\t".$higher_TFG_k_t_1.','.$higher_TMG_k_t_1."\t".$higher_1_type."\t".$higher_TFG_k_t_2.','.$higher_TMG_k_t_2."\t".$higher_2_type."\t".$higher_TFG_k_t_3.','.$higher_TMG_k_t_3."\t".$higher_3_type."\t".$number_error_1."\t".$number_error_2."\t".$number_error_3."\t".$number_Y_error."\t".$Parameters_ref->{'number_individuals_with_aberrant_reads'}) ;
					# print genotypes and expressions levels
					print(DETAIL "\t".$Parameters_ref->{'heterogametic_parent_genotype'}."\t".$Parameters_ref->{'heterogametic_parent_expression'}."\t".$Parameters_ref->{'homogametic_parent_genotype'}."\t".$Parameters_ref->{'homogametic_parent_expression'}) ;
					foreach my $hom_ind (@homogametic) {
						print(DETAIL "\t".$homogametic_genotypes_ref->{$hom_ind}."\t".$homogametic_expression{$hom_ind}) ;
					}
					foreach my $het_ind (@heterogametic) {
						print(DETAIL "\t".$heterogametic_genotypes_ref->{$het_ind}."\t".$heterogametic_expression{$het_ind}) ;
					}
					print(DETAIL "\n") ;
					close(DETAIL) ;
				}
				# prepare sex-linked details
				if (defined($Parameters_ref->{'detail-sex-linked'})) {
					prepare_sex_linked_details($Parameters_ref, $S_k_t_ref, $higher_2_type, $higher_3_type, $higher_TFG_k_t_2, $higher_TMG_k_t_2, $higher_TFG_k_t_3, $higher_TMG_k_t_3, \%homogametic_expression, \%heterogametic_expression, $homogametic_genotypes_ref, $heterogametic_genotypes_ref, $number_error_1, $number_error_2, $number_error_3, $number_Y_error) ;
				}
			} else {
				if (defined($Parameters_ref->{'detail'})) {
					open(DETAIL, ">>".$Parameters_ref->{'detail_file_name'}) ;
					print(DETAIL $Parameters_ref->{'contig_name'}."\t".$Parameters_ref->{'current_contig_position'}."\t".$S_k_t_ref->{1}."\t".$S_k_t_ref->{2}."\t".$S_k_t_ref->{3}."\t".$higher_TFG_k_t_1.','.$higher_TMG_k_t_1."\t".$higher_1_type."\t".$higher_TFG_k_t_2.','.$higher_TMG_k_t_2."\t".$higher_2_type."\t".$higher_TFG_k_t_3.','.$higher_TMG_k_t_3."\t".$higher_3_type."\t".$number_error_1."\t".$number_error_2."\t".$number_error_3."\t".$number_Y_error) ;
					# print genotype
					print(DETAIL "\t".$Parameters_ref->{'heterogametic_parent_genotype'}."\t".$Parameters_ref->{'homogametic_parent_genotype'}) ;
					foreach my $hom_ind (@homogametic) {
						print(DETAIL "\t".$homogametic_genotypes_ref->{$hom_ind}) ;
					}
					foreach my $het_ind (@heterogametic) {
						print(DETAIL "\t".$heterogametic_genotypes_ref->{$het_ind}) ;
					}
					print(DETAIL "\n") ;
					close(DETAIL) ;
				}
			}
		}

		# if needed the sequences are retrieved
		if ((!defined($Parameters_ref->{'no_sex_chr'}))&&(defined($Parameters_ref->{'sequences'}))&&($Parameters{'parameters_estimations_over'} eq "yes")) {
			if (($S_k_t_ref->{2} > $S_k_t_ref->{1})&&($S_k_t_ref->{2} > $S_k_t_ref->{3})&&($higher_2_type ne 'monomorphic')&&($higher_2_type ne 'not_informative')) {
				# XY informative SNP
				$Parameters_ref->{'current_contig_sequence_X1'} .= $lambda_2_ref->{$higher_TFG_k_t_2.','.$higher_TMG_k_t_2}->{'X_alleles'}->[0] ;
				$Parameters_ref->{'current_contig_sequence_X2'} .= $lambda_2_ref->{$higher_TFG_k_t_2.','.$higher_TMG_k_t_2}->{'X_alleles'}->[1] ;
				$Parameters_ref->{'current_contig_sequence_X3'} .= $lambda_2_ref->{$higher_TFG_k_t_2.','.$higher_TMG_k_t_2}->{'X_alleles'}->[2] ;
				$Parameters_ref->{'current_contig_sequence_Y'} .= $lambda_2_ref->{$higher_TFG_k_t_2.','.$higher_TMG_k_t_2}->{'Y_allele'} ;
			} elsif (($S_k_t_ref->{3} > $S_k_t_ref->{1})&&($S_k_t_ref->{3} > $S_k_t_ref->{2})&&($higher_3_type ne 'monomorphic')&&($higher_3_type ne 'not_informative')) {
				# X hemizygous informative SNP
				$Parameters_ref->{'current_contig_sequence_X1'} .= $lambda_3_ref->{$higher_TFG_k_t_3.','.$higher_TMG_k_t_3}->{'X_alleles'}->[0] ;
				$Parameters_ref->{'current_contig_sequence_X2'} .= $lambda_3_ref->{$higher_TFG_k_t_3.','.$higher_TMG_k_t_3}->{'X_alleles'}->[1] ;
				$Parameters_ref->{'current_contig_sequence_X3'} .= $lambda_3_ref->{$higher_TFG_k_t_3.','.$higher_TMG_k_t_3}->{'X_alleles'}->[2] ;
				$Parameters_ref->{'current_contig_sequence_Y'} .= 'N' ;
			} else {
				# non informative SNP or autosomal SNP, do not write it in sex-linked sequences
				if (defined($Parameters_ref->{'ALR_file'})) {
					my @current_alr_line_split = split(/\t/, $Parameters_ref->{'current_alr_line'}) ;
					$Parameters_ref->{'current_contig_sequence_X1'} .= $current_alr_line_split[0] ;
					$Parameters_ref->{'current_contig_sequence_X2'} .= $current_alr_line_split[0] ;
					$Parameters_ref->{'current_contig_sequence_X3'} .= $current_alr_line_split[0] ;
					$Parameters_ref->{'current_contig_sequence_Y'} .= $current_alr_line_split[0] ;
				}
			}
		}

		# Check if position is polymorphic and informative
		if ((($higher_1_type eq 'monomorphic')&&($S_k_t_ref->{1} > $S_k_t_ref->{2})&&($S_k_t_ref->{1} > $S_k_t_ref->{3}))||(($higher_2_type eq 'monomorphic')&&($S_k_t_ref->{2} > $S_k_t_ref->{1})&&($S_k_t_ref->{2} > $S_k_t_ref->{3}))||(($higher_3_type eq 'monomorphic')&&($S_k_t_ref->{3} > $S_k_t_ref->{1})&&($S_k_t_ref->{3} > $S_k_t_ref->{2}))) {
			# monomorphic position
		} else {
			if ((($higher_1_type eq 'not_informative')&&($S_k_t_ref->{1} > $S_k_t_ref->{2})&&($S_k_t_ref->{1} > $S_k_t_ref->{3}))||(($higher_2_type eq 'not_informative')&&($S_k_t_ref->{2} > $S_k_t_ref->{1})&&($S_k_t_ref->{2} > $S_k_t_ref->{3}))) {
				# non informative
			} else {
				# informative position for contig segregation pattern
				$Parameters_ref->{'S_k_1'} += $S_k_1_factor*$S_k_t_ref->{1} ;
				$Parameters_ref->{'S_k_2'} += $S_k_2_factor*$S_k_t_ref->{2} ;
				$Parameters_ref->{'S_k_3'} += $S_k_3_factor*$S_k_t_ref->{3} ;
				if (($S_k_t_ref->{1} > $S_k_t_ref->{2})&&($S_k_t_ref->{1} > $S_k_t_ref->{3})) {
					# position inferred as autosomal
					if ($number_error_1 == 0) {
						$Parameters_ref->{'autosomal_without_error'} ++ ;
					} else {
						$Parameters_ref->{'autosomal_with_error'} ++ ;
					}
				} elsif (($S_k_t_ref->{2} > $S_k_t_ref->{1})&&($S_k_t_ref->{2} > $S_k_t_ref->{3})) {
					# position inferred as X/Y
					if (($number_Y_error == 0)&&($number_error_2 == 0)) {
						$Parameters_ref->{'XY_without_error'} ++ ;
						if (defined($Parameters_ref->{'ALR_file'})) {
							if ($Parameters_ref->{'number_individuals_with_aberrant_reads'} == 0) {
								$Parameters_ref->{'number_clean_XY_SNPs_without_error'} ++ ;
							}
						} else {
							$Parameters_ref->{'number_clean_XY_SNPs_without_error'} ++ ;
						}
					} else {
						$Parameters_ref->{'XY_with_error'} ++ ;
					}
				} elsif (($S_k_t_ref->{3} > $S_k_t_ref->{1})&&($S_k_t_ref->{3} > $S_k_t_ref->{2})) {
					# position inferred as hemizygous
					if ($number_error_3 == 0) {
						$Parameters_ref->{'hemizygous_without_error'} ++ ;
						if (defined($Parameters_ref->{'ALR_file'})) {
							if ($Parameters_ref->{'number_individuals_with_aberrant_reads'} == 0) {
								$Parameters_ref->{'number_clean_hemizygous_SNPs_without_error'} ++ ;
							}
						} else {
							$Parameters_ref->{'number_clean_hemizygous_SNPs_without_error'} ++ ;
						}
					} else {
						$Parameters_ref->{'hemizygous_with_error'} ++ ;
					}
				}
			}
		}
	} else {
		# The position could not be studied due to a lack of individuals or monomorphic position, if needed the sequence is retrieved
		if ((defined($Parameters_ref->{'sequences'}))&&($Parameters{'parameters_estimations_over'} eq "yes")) {
			if (defined($Parameters_ref->{'ALR_file'})) {
				my @current_alr_line_split = split(/\t/, $Parameters_ref->{'current_alr_line'}) ;
				$Parameters_ref->{'current_contig_sequence_X1'} .= $current_alr_line_split[0] ;
				$Parameters_ref->{'current_contig_sequence_X2'} .= $current_alr_line_split[0] ;
				$Parameters_ref->{'current_contig_sequence_X3'} .= $current_alr_line_split[0] ;
				$Parameters_ref->{'current_contig_sequence_Y'} .= $current_alr_line_split[0] ;
			}
		}
	}

	return 0;
}


#--------------------------
# retrieve base expression
#--------------------------
sub retrieve_base_expression {

	# retrieve parameters
	my $base = shift ;
	my $alr_expression = shift ;

	# initialize
	my $expr ;
	my @alr_expression_split = split(/\//, $alr_expression) ;

	if ($base eq 'A') {
		$expr = $alr_expression_split[0] ;
	} elsif ($base eq 'C') {
		$expr = $alr_expression_split[1] ;
	} elsif ($base eq 'G') {
		$expr = $alr_expression_split[2] ;
	} else {
		# $base eq T
		$expr = $alr_expression_split[3] ;
	}

	return $expr ;
}


#---------------------------
# retrieve total expression
#---------------------------
sub retrieve_total_expression {

	# retrieve parameters
	my $alr_expression = shift ;

	# split
	my @alr_expression_split = split(/\//, $alr_expression) ;

	my $total_expr = $alr_expression_split[0] + $alr_expression_split[1] + $alr_expression_split[2] + $alr_expression_split[3] ;

	return $total_expr ;
}


#-----------------
# random Bernoulli
#-----------------
sub random_Bernoulli {

	# retrieve parameters
	my $arg = shift ;

	# create ouptut random Bernoulli
	my $output ;

	my $random = rand(1) ;
	if ($random < $arg) {
		$output = 1 ;
	} else {
		$output = 0 ;
	}

	return $output ;
}


#---------------------
# random multinomial 3
#---------------------
sub random_multinomial_3 {

	# retrieve parameters
	my @args = @_ ;

	# create ouptut random multinomial
	my @output ;

	my $random = rand(1) ;
	if ($random < $args[0]) {
		@output = (1, 0, 0) ;
	} elsif ($random < ($args[0] + $args[1])) {
		@output = (0, 1, 0) ;
	} else {
		@output = (0, 0, 1) ;
	}

	return @output ;
}


#---------------------
# random multinomial 4
#---------------------
sub random_multinomial_4 {

	# retrieve parameters
	my @args = @_ ;

	# create ouptut random multinomial
	my @output ;

	my $random = rand(1) ;
	if ($random < $args[0]) {
		@output = (1, 0, 0, 0) ;
	} elsif ($random < ($args[0] + $args[1])) {
		@output = (0, 1, 0, 0) ;
	} elsif ($random < ($args[0] + $args[1] + $args[2])) {
		@output = (0, 0, 1, 0) ;
	} else {
		@output = (0, 0, 0, 1) ;
	}

	return @output ;
}


#----------------------
# random multinomial 10
#----------------------
sub random_multinomial_10 {

	# retrieve parameters
	my @args = @_ ;

	# create ouptut random multinomial
	my @output ;

	my $random = rand(1) ;
	if ($random < $args[0]) {
		@output = (1, 0, 0, 0, 0, 0, 0, 0, 0, 0) ;
	} elsif ($random < ($args[0] + $args[1])) {
		@output = (0, 1, 0, 0, 0, 0, 0, 0, 0, 0) ;
	} elsif ($random < ($args[0] + $args[1] + $args[2])) {
		@output = (0, 0, 1, 0, 0, 0, 0, 0, 0, 0) ;
	} elsif ($random < ($args[0] + $args[1] + $args[2] + $args[3])) {
		@output = (0, 0, 0, 1, 0, 0, 0, 0, 0, 0) ;
	} elsif ($random < ($args[0] + $args[1] + $args[2] + $args[3] + $args[4])) {
		@output = (0, 0, 0, 0, 1, 0, 0, 0, 0, 0) ;
	} elsif ($random < ($args[0] + $args[1] + $args[2] + $args[3] + $args[4] + $args[5])) {
		@output = (0, 0, 0, 0, 0, 1, 0, 0, 0, 0) ;
	} elsif ($random < ($args[0] + $args[1] + $args[2] + $args[3] + $args[4] + $args[5] + $args[6])) {
		@output = (0, 0, 0, 0, 0, 0, 1, 0, 0, 0) ;
	} elsif ($random < ($args[0] + $args[1] + $args[2] + $args[3] + $args[4] + $args[5] + $args[6] + $args[7])) {
		@output = (0, 0, 0, 0, 0, 0, 0, 1, 0, 0) ;
	} elsif ($random < ($args[0] + $args[1] + $args[2] + $args[3] + $args[4] + $args[5] + $args[6] + $args[7] + $args[8])) {
		@output = (0, 0, 0, 0, 0, 0, 0, 0, 1, 0) ;
	} else {
		@output = (0, 0, 0, 0, 0, 0, 0, 0, 0, 1) ;
	}

	return @output ;
}


#----------------------
# random multinomial 12
#----------------------
sub random_multinomial_12 {

	# retrieve parameters
	my @args = @_ ;

	# create ouptut random multinomial
	my @output ;

	my $random = rand(1) ;
	if ($random < $args[0]) {
		@output = (1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0) ;
	} elsif ($random < ($args[0] + $args[1])) {
		@output = (0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0) ;
	} elsif ($random < ($args[0] + $args[1] + $args[2])) {
		@output = (0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0) ;
	} elsif ($random < ($args[0] + $args[1] + $args[2] + $args[3])) {
		@output = (0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0) ;
	} elsif ($random < ($args[0] + $args[1] + $args[2] + $args[3] + $args[4])) {
		@output = (0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0) ;
	} elsif ($random < ($args[0] + $args[1] + $args[2] + $args[3] + $args[4] + $args[5])) {
		@output = (0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0) ;
	} elsif ($random < ($args[0] + $args[1] + $args[2] + $args[3] + $args[4] + $args[5] + $args[6])) {
		@output = (0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0) ;
	} elsif ($random < ($args[0] + $args[1] + $args[2] + $args[3] + $args[4] + $args[5] + $args[6] + $args[7])) {
		@output = (0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0) ;
	} elsif ($random < ($args[0] + $args[1] + $args[2] + $args[3] + $args[4] + $args[5] + $args[6] + $args[7] + $args[8])) {
		@output = (0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0) ;
	} elsif ($random < ($args[0] + $args[1] + $args[2] + $args[3] + $args[4] + $args[5] + $args[6] + $args[7] + $args[8] + $args[9])) {
		@output = (0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0) ;
	} elsif ($random < ($args[0] + $args[1] + $args[2] + $args[3] + $args[4] + $args[5] + $args[6] + $args[7] + $args[8] + $args[9] + $args[10])) {
		@output = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0) ;
	} else {
		@output = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1) ;
	}

	return @output ;
}


#----------------
# Compute %GE_k_t
#----------------
sub Compute_GE_YGE_and_factorise {

	# Recovering parameters
	my $Parameters_ref = shift ;
	my $type = shift ;
	my ($alpha_ref, $beta_2_ref, $beta_3_ref, $pi_ref, $E, $p) ;
	if ($type eq 'old') {
		$alpha_ref = $Parameters_ref->{'alpha_old'} ;
		$beta_2_ref = $Parameters_ref->{'beta_2_old'} ;
		$beta_3_ref = $Parameters_ref->{'beta_3_old'} ;
		$pi_ref = $Parameters_ref->{'pi_old'} ;
		$p = $Parameters_ref->{'p_old'} ;
		$E = $Parameters_ref->{'E_old'} ;
	} else {
		$alpha_ref = $Parameters_ref->{'alpha'} ;
		$beta_2_ref = $Parameters_ref->{'beta_2'} ;
		$beta_3_ref = $Parameters_ref->{'beta_3'} ;
		$pi_ref = $Parameters_ref->{'pi'} ;
		$p = $Parameters_ref->{'p'} ;
		$E = $Parameters_ref->{'E'} ;
	}
	my $lambda_1_ref = $Parameters_ref->{'lambda_1'} ;
	my $lambda_1_E_ref = $Parameters_ref->{'lambda_1_E'} ;
	my $lambda_2_ref = $Parameters_ref->{'lambda_2'} ;
	my $lambda_2_E_ref = $Parameters_ref->{'lambda_2_E'} ;
	my $lambda_2_p_ref = $Parameters_ref->{'lambda_2_p'} ;
	my $lambda_2_E_p_ref = $Parameters_ref->{'lambda_2_E_p'} ;
	my $lambda_3_ref = $Parameters_ref->{'lambda_3'} ;
	my $lambda_3_E_ref = $Parameters_ref->{'lambda_3_E'} ;
	my $mu_2_het_ref = $Parameters_ref->{'mu_2_het'} ;
	my $mu_2_p_het_ref = $Parameters_ref->{'mu_2_p_het'} ;
	my $mu_2_E_het_ref = $Parameters_ref->{'mu_2_E_het'} ;
	my $mu_2_E_p_het_ref = $Parameters_ref->{'mu_2_E_p_het'} ;
	my $mu_3_het_ref = $Parameters_ref->{'mu_3_het'} ;
	my $mu_3_E_het_ref = $Parameters_ref->{'mu_3_E_het'} ;
	my $mu_hom_ref = $Parameters_ref->{'mu_hom'} ;
	my $mu_E_hom_ref = $Parameters_ref->{'mu_E_hom'} ;
	my $heterogametic_genotypes_ref = $Parameters_ref->{'heterogametic_genotypes'} ;
	my $homogametic_genotypes_ref = $Parameters_ref->{'homogametic_genotypes'} ;

	my (%GE_k_t_1, %GE_k_t_2, %GE_k_t_3, %YGE_k_t_2, %factor) ;

	# heterogametic progeny
	foreach my $i_het (keys %{$heterogametic_genotypes_ref}) {
		my $l = $heterogametic_genotypes_ref->{$i_het} ;
		if ($l ne 'NN') {
			# j=1
			foreach my $n (keys %{$alpha_ref}) {
				foreach my $m (keys %{$alpha_ref}) {
					my @terms = ((1-$E)*$lambda_1_ref->{$n.','.$m}->{$l}, $E*$lambda_1_E_ref->{$n.','.$m}->{$l}) ;
					my $denominator = $terms[0] + $terms[1] ;
					# GE
					$GE_k_t_1{$n.','.$m}{$i_het} = $terms[1] / $denominator ;
					# factorisation
					my $max = max(@terms) ;
					$factor{'1'}{$n.','.$m} += log($max) + log($denominator/$max) ;
				}
			}
			# j=2
			if ($pi_ref->{2}!= 0) {
				foreach my $n (keys %{$beta_2_ref}) {
					foreach my $m (keys %{$alpha_ref}) {
						my @terms = ( (1-$E)*(1-$p)*$lambda_2_ref->{$n.','.$m}->{'het,'.$l}, $E*(1-$p)*$lambda_2_E_ref->{$n.','.$m}->{'het,'.$l}, (1-$E)*$p*$lambda_2_p_ref->{$n.','.$m}->{'het,'.$l}, $E*$p*$lambda_2_E_p_ref->{$n.','.$m}->{'het,'.$l}) ;
						my $denominator = $terms[0] + $terms[1] + $terms[2] + $terms[3] ;
						# GE, YGE
						if ($p!= 0) {
							if (($lambda_2_ref->{$n.','.$m}->{'het,'.$l} == 0)&&($lambda_2_E_ref->{$n.','.$m}->{'het,'.$l}!= 0)&&($lambda_2_p_ref->{$n.','.$m}->{'het,'.$l}!= 0)) {
								# a Y genotyping error is the most likely
								$GE_k_t_2{$n.','.$m}{$i_het} = 0 ;
								$YGE_k_t_2{$n.','.$m}{$i_het} = 1 ;
							} elsif (($lambda_2_ref->{$n.','.$m}->{'het,'.$l}!= 0)&&($lambda_2_E_ref->{$n.','.$m}->{'het,'.$l}!= 0)&&($lambda_2_p_ref->{$n.','.$m}->{'het,'.$l} == 0)&&($lambda_2_E_p_ref->{$n.','.$m}->{'het,'.$l}!= 0)) {
								# No genotyping error is the most likely
								$GE_k_t_2{$n.','.$m}{$i_het} = 0 ;
								$YGE_k_t_2{$n.','.$m}{$i_het} = 0 ;
							} elsif (($lambda_2_ref->{$n.','.$m}->{'het,'.$l}!= 0)&&($lambda_2_E_ref->{$n.','.$m}->{'het,'.$l} == 0)&&($lambda_2_p_ref->{$n.','.$m}->{'het,'.$l} == 0)&&($lambda_2_E_p_ref->{$n.','.$m}->{'het,'.$l}!= 0)) {
								# No genotyping error is the most likely
								$GE_k_t_2{$n.','.$m}{$i_het} = 0 ;
								$YGE_k_t_2{$n.','.$m}{$i_het} = 0 ;
							} else {
								$GE_k_t_2{$n.','.$m}{$i_het} = ($terms[1] + $terms[3]) / $denominator ;
								$YGE_k_t_2{$n.','.$m}{$i_het} = ($terms[2] + $terms[3]) / $denominator ;
							}
						} else {
							$GE_k_t_2{$n.','.$m}{$i_het} = ($terms[1] + $terms[3]) / $denominator ;
							$YGE_k_t_2{$n.','.$m}{$i_het} = ($terms[2] + $terms[3]) / $denominator ;
						}
						# factorisation
						my $max = max(@terms) ;
						$factor{'2'}{$n.','.$m} += log($max) + log($denominator/$max) ;
					}
				}
			}
			# j=3
			if ($pi_ref->{3}!= 0) {
				foreach my $n (keys %{$beta_3_ref}) {
					foreach my $m (keys %{$alpha_ref}) {
						my @terms = ( (1-$E)*$lambda_3_ref->{$n.','.$m}->{'het,'.$l}, $E*$lambda_3_E_ref->{$n.','.$m}->{'het,'.$l} ) ;
						my $denominator = $terms[0] + $terms[1] ;
						# GE
						$GE_k_t_3{$n.','.$m}{$i_het} = $terms[1] / $denominator ;
						# factorisation
						my $max = max(@terms) ;
						$factor{'3'}{$n.','.$m} += log($max) + log($denominator/$max) ;
					}
				}
			}
		}
	}

	# homogametic progeny
	foreach my $i_hom (keys %{$homogametic_genotypes_ref}) {
		my $l = $homogametic_genotypes_ref->{$i_hom} ;
		if ($l ne 'NN') {
			# j=1
			foreach my $n (keys %{$alpha_ref}) {
				foreach my $m (keys %{$alpha_ref}) {
					my @terms = ( (1-$E)*$lambda_1_ref->{$n.','.$m}->{$l}, $E*$lambda_1_E_ref->{$n.','.$m}->{$l} ) ;
					my $denominator = $terms[0] + $terms[1] ;
					# GE
					$GE_k_t_1{$n.','.$m}{$i_hom} = $terms[1] / $denominator ;
					# factorisation
					my $max = max(@terms) ;
					$factor{'1'}{$n.','.$m} += log($max) + log($denominator/$max) ;
				}
			}
			# j=2
			if ($pi_ref->{2}!= 0) {
				foreach my $n (keys %{$beta_2_ref}) {
					foreach my $m (keys %{$alpha_ref}) {
						my @terms = ( (1-$E)*$lambda_2_ref->{$n.','.$m}->{'hom,'.$l}, $E*$lambda_2_E_ref->{$n.','.$m}->{'hom,'.$l} ) ;
						my $denominator = $terms[0] + $terms[1] ;
						# GE
						$GE_k_t_2{$n.','.$m}{$i_hom} = $terms[1] / $denominator ;
						# factorisation
						my $max = max(@terms) ;
						$factor{'2'}{$n.','.$m} += log($max) + log($denominator/$max) ;
					}
				}
			}
			# j=3
			if ($pi_ref->{3}!= 0) {
				foreach my $n (keys %{$beta_3_ref}) {
					foreach my $m (keys %{$alpha_ref}) {
						my @terms = ( (1-$E)*$lambda_3_ref->{$n.','.$m}->{'hom,'.$l}, $E*$lambda_3_E_ref->{$n.','.$m}->{'hom,'.$l} ) ;
						my $denominator = $terms[0] + $terms[1] ;
						# GE
						$GE_k_t_3{$n.','.$m}{$i_hom} = $terms[1] / $denominator ;
						# factorisation
						my $max = max(@terms) ;
						$factor{'3'}{$n.','.$m} += log($max) + log($denominator/$max) ;
					}
				}
			}
		}
	}
	# heterogametic parent
	my $l = $Parameters_ref->{'heterogametic_parent_genotype'} ;
	# j = 1
	foreach my $n (keys %{$alpha_ref}) {
		foreach my $m (keys %{$alpha_ref}) {
			my @terms = ( (1-$E)*$mu_hom_ref->{$n}->{$l}, $E*$mu_E_hom_ref->{$n}->{$l} ) ;
			my $denominator = $terms[0] + $terms[1] ;
			# GE
			$GE_k_t_1{$n.','.$m}{$Parameters_ref->{'heterogametic_parent_name'}} = $terms[1] / $denominator ;
			# factorisation
			my $max = max(@terms) ;
			$factor{'1'}{$n.','.$m} += log($max) + log($denominator/$max) ;
		}
	}
	# j = 2
	if ($pi_ref->{2}!= 0) {
		foreach my $n (keys %{$beta_2_ref}) {
			foreach my $m (keys %{$alpha_ref}) {
				my @terms = ( (1-$E)*(1-$p)*$mu_2_het_ref->{$n}->{$l}, $E*(1-$p)*$mu_2_E_het_ref->{$n}->{$l}, (1-$E)*$p*$mu_2_p_het_ref->{$n}->{$l}, $E*$p*$mu_2_E_p_het_ref->{$n}->{$l}) ;
				my $denominator = $terms[0] + $terms[1] + $terms[2] + $terms[3] ;
				# GE, YGE
				if (($mu_2_het_ref->{$n}->{$l} == 0)&&($mu_2_E_het_ref->{$n}->{$l}!= 0)&&($mu_2_p_het_ref->{$n}->{$l}!= 0)) {
					# a Y genotyping error is the most likely
					$GE_k_t_2{$n.','.$m}{$Parameters_ref->{'heterogametic_parent_name'}} = 0 ;
					$YGE_k_t_2{$n.','.$m}{$Parameters_ref->{'heterogametic_parent_name'}} = 1 ;
				} elsif (($mu_2_het_ref->{$n}->{$l}!= 0)&&($mu_2_E_p_het_ref->{$n}->{$l}!= 0)) {
					# No genotyping error is the most likely
					$GE_k_t_2{$n.','.$m}{$Parameters_ref->{'heterogametic_parent_name'}} = 0 ;
					$YGE_k_t_2{$n.','.$m}{$Parameters_ref->{'heterogametic_parent_name'}} = 0 ;
				} else {
					$GE_k_t_2{$n.','.$m}{$Parameters_ref->{'heterogametic_parent_name'}} = ($terms[1] + $terms[3]) / $denominator ;
					$YGE_k_t_2{$n.','.$m}{$Parameters_ref->{'heterogametic_parent_name'}} = ($terms[2] + $terms[3]) / $denominator ;
				}
				# factorisation
				my $max = max(@terms) ;
				$factor{'2'}{$n.','.$m} += log($max) + log($denominator/$max) ;
			}
		}
	}
	# j = 3
	if ($pi_ref->{3}!= 0) {
		foreach my $n (keys %{$beta_3_ref}) {
			foreach my $m (keys %{$alpha_ref}) {
				my @terms = ( (1-$E)*$mu_3_het_ref->{$n}->{$l}, $E*$mu_3_E_het_ref->{$n}->{$l} ) ;
				my $denominator = $terms[0] + $terms[1] ;
				# GE
				$GE_k_t_3{$n.','.$m}{$Parameters_ref->{'heterogametic_parent_name'}} = $terms[1] / $denominator ;
				# factorisation
				my $max = max(@terms) ;
				$factor{'3'}{$n.','.$m} += log($max) + log($denominator/$max) ;
			}
		}
	}
	# homogametic parent
	$l = $Parameters_ref->{'homogametic_parent_genotype'} ;
	# j = 1
	foreach my $n (keys %{$alpha_ref}) {
		foreach my $m (keys %{$alpha_ref}) {
			my @terms = ( (1-$E)*$mu_hom_ref->{$m}->{$l}, $E*$mu_E_hom_ref->{$m}->{$l} ) ;
			my $denominator = $terms[0] + $terms[1] ;
			# GE
			$GE_k_t_1{$n.','.$m}{$Parameters_ref->{'homogametic_parent_name'}} = $terms[1] / $denominator ;
			# factorisation
			my $max = max(@terms) ;
			$factor{'1'}{$n.','.$m} += log($max) + log($denominator/$max) ;
		}
	}
	# j = 2
	if ($pi_ref->{2}!= 0) {
		foreach my $n (keys %{$beta_2_ref}) {
			foreach my $m (keys %{$alpha_ref}) {
				my @terms = ( (1-$E)*$mu_hom_ref->{$m}->{$l}, $E*$mu_E_hom_ref->{$m}->{$l} ) ;
				my $denominator = $terms[0] + $terms[1] ;
				# GE
				$GE_k_t_2{$n.','.$m}{$Parameters_ref->{'homogametic_parent_name'}} = $terms[1] / $denominator ;
				# factorisation
				my $max = max(@terms) ;
				$factor{'2'}{$n.','.$m} += log($max) + log($denominator/$max) ;
			}
		}
	}
	# j = 3
	if ($pi_ref->{3}!= 0) {
		foreach my $n (keys %{$beta_3_ref}) {
			foreach my $m (keys %{$alpha_ref}) {
				my @terms = ( (1-$E)*$mu_hom_ref->{$m}->{$l}, $E*$mu_E_hom_ref->{$m}->{$l} ) ;
				my $denominator = $terms[0] + $terms[1] ;
				# GE
				$GE_k_t_3{$n.','.$m}{$Parameters_ref->{'homogametic_parent_name'}} = $terms[1] / $denominator ;
				# factorisation
				my $max = max(@terms) ;
				$factor{'3'}{$n.','.$m} += log($max) + log($denominator/$max) ;
			}
		}
	}

	return (\%GE_k_t_1, \%GE_k_t_2, \%GE_k_t_3, \%YGE_k_t_2, \%factor) ;
}


#----------------------
# Compute %S_k_t{j}
#----------------------
sub Compute_S_TMG_TFG_k_t_j {

	# Recovering parameters
	my $Parameters_ref = shift ;
	my $type = shift ;
	my $factor_ref = shift ;
	my ($alpha_ref, $beta_2_ref, $beta_3_ref, $pi_ref, $E, $p) ;
	if ($type eq 'old') {
		$alpha_ref = $Parameters_ref->{'alpha_old'} ;
		$beta_2_ref = $Parameters_ref->{'beta_2_old'} ;
		$beta_3_ref = $Parameters_ref->{'beta_3_old'} ;
		$pi_ref = $Parameters_ref->{'pi_old'} ;
		$p = $Parameters_ref->{'p_old'} ;
		$E = $Parameters_ref->{'E_old'} ;
	} else {
		$alpha_ref = $Parameters_ref->{'alpha'} ;
		$beta_2_ref = $Parameters_ref->{'beta_2'} ;
		$beta_3_ref = $Parameters_ref->{'beta_3'} ;
		$pi_ref = $Parameters_ref->{'pi'} ;
		$p = $Parameters_ref->{'p'} ;
		$E = $Parameters_ref->{'E'} ;
	}

	# initialisations
	my (%S_k_t, %TMG_k_t_1, %TMG_k_t_2, %TMG_k_t_3, %TFG_k_t_1, %TFG_k_t_2, %TFG_k_t_3, %Bjm, %Ajn) ;
	my %Dj = (1=>0, 2=>0, 3=>0) ;
	my $D = 0 ;

	# Compute B_{jm} = \sum_{n} \beta_{jn} exp(P_{jnm})  where P_{jnm} = factor{j}{n,m}
	# j = 1
	foreach my $m (keys %{$alpha_ref}) {
		my $n_max = 'AC' ;
		foreach my $n (keys %{$alpha_ref}) {
			if (($alpha_ref->{$n}*exp($factor_ref->{'1'}{$n.','.$m})) > ($alpha_ref->{$n_max}*exp($factor_ref->{'1'}{$n_max.','.$m}))) {
				$n_max = $n ;
			}
		}
		my $sum = 0 ;
		foreach my $n (keys %{$alpha_ref}) {
			$sum += ($alpha_ref->{$n} / $alpha_ref->{$n_max}) * exp($factor_ref->{'1'}{$n.','.$m} - $factor_ref->{'1'}{$n_max.','.$m}) ;
		}
		$Bjm{'1'}{$m} = $alpha_ref->{$n_max}*exp($factor_ref->{'1'}{$n_max.','.$m}) * $sum ;
	}
	# j = 2
	if ($pi_ref->{2}!= 0) {
		foreach my $m (keys %{$alpha_ref}) {
			my $n_max = 'AC' ;
			foreach my $n (keys %{$beta_2_ref}) {
				if (($beta_2_ref->{$n}*exp($factor_ref->{'2'}{$n.','.$m})) > ($beta_2_ref->{$n_max}*exp($factor_ref->{'2'}{$n_max.','.$m}))) {
					$n_max = $n ;
				}
			}
			my $sum = 0 ;
			foreach my $n (keys %{$beta_2_ref}) {
				$sum += ($beta_2_ref->{$n} / $beta_2_ref->{$n_max}) * exp($factor_ref->{'2'}{$n.','.$m} - $factor_ref->{'2'}{$n_max.','.$m}) ;
			}
			$Bjm{'2'}{$m} = $beta_2_ref->{$n_max}*exp($factor_ref->{'2'}{$n_max.','.$m}) * $sum ;
		}
	}
	# j = 3
	if ($pi_ref->{3}!= 0) {
		foreach my $m (keys %{$alpha_ref}) {
			my $n_max = 'A' ;
			foreach my $n (keys %{$beta_3_ref}) {
				if (($beta_3_ref->{$n}*exp($factor_ref->{'3'}{$n.','.$m})) > ($beta_3_ref->{$n_max}*exp($factor_ref->{'3'}{$n_max.','.$m}))) {
					$n_max = $n ;
				}
			}
			my $sum = 0 ;
			foreach my $n (keys %{$beta_3_ref}) {
				$sum += ($beta_3_ref->{$n} / $beta_3_ref->{$n_max}) * exp($factor_ref->{'3'}{$n.','.$m} - $factor_ref->{'3'}{$n_max.','.$m}) ;
			}
			$Bjm{'3'}{$m} = $beta_3_ref->{$n_max}*exp($factor_ref->{'3'}{$n_max.','.$m}) * $sum ;
		}
	}

	# Compute A_{jn} = \sum_{m} \alpha_{m} exp(P_{jnm})  where P_{jnm} = factor{j}{n,m}
	# j=1
	foreach my $n (keys %{$alpha_ref}) {
		my $m_max = 'AA' ;
		foreach my $m (keys %{$alpha_ref}) {
			if (($alpha_ref->{$m} * exp($factor_ref->{'1'}{$n.','.$m})) > ($alpha_ref->{$m_max} * exp($factor_ref->{'1'}{$n.','.$m_max}))) {
				$m_max = $m ;
			}
		}
		my $sum = 0 ;
		foreach my $m (keys %{$alpha_ref}) {
			$sum += ($alpha_ref->{$m} / $alpha_ref->{$m_max}) * exp($factor_ref->{'1'}{$n.','.$m} - $factor_ref->{'1'}{$n.','.$m_max}) ;
		}
		$Ajn{'1'}{$n} = $alpha_ref->{$m_max} * exp($factor_ref->{'1'}{$n.','.$m_max}) * $sum ;
	}
	# j=2
	if ($pi_ref->{2}!= 0) {
		foreach my $n (keys %{$beta_2_ref}) {
			my $m_max = 'AA' ;
			foreach my $m (keys %{$alpha_ref}) {
				if (($alpha_ref->{$m} * exp($factor_ref->{'2'}{$n.','.$m})) > ($alpha_ref->{$m_max} * exp($factor_ref->{'2'}{$n.','.$m_max}))) {
					$m_max = $m ;
				}
			}
			my $sum = 0 ;
			foreach my $m (keys %{$alpha_ref}) {
				$sum += ($alpha_ref->{$m} / $alpha_ref->{$m_max}) * exp($factor_ref->{'2'}{$n.','.$m} - $factor_ref->{'2'}{$n.','.$m_max}) ;
			}
			$Ajn{'2'}{$n} = $alpha_ref->{$m_max} * exp($factor_ref->{'2'}{$n.','.$m_max}) * $sum ;
		}
	}
	# j=3
	if ($pi_ref->{3}!= 0) {
		foreach my $n (keys %{$beta_3_ref}) {
			my $m_max = 'AA' ;
			foreach my $m (keys %{$alpha_ref}) {
				if (($alpha_ref->{$m} * exp($factor_ref->{'3'}{$n.','.$m})) > ($alpha_ref->{$m_max} * exp($factor_ref->{'3'}{$n.','.$m_max}))) {
					$m_max = $m ;
				}
			}
			my $sum = 0 ;
			foreach my $m (keys %{$alpha_ref}) {
				$sum += ($alpha_ref->{$m} / $alpha_ref->{$m_max}) * exp($factor_ref->{'3'}{$n.','.$m} - $factor_ref->{'3'}{$n.','.$m_max}) ;
			}
			$Ajn{'3'}{$n} = $alpha_ref->{$m_max} * exp($factor_ref->{'3'}{$n.','.$m_max}) * $sum ;
		}
	}

	# Compute D_j = \sum_{m} \alpha_m Bjm
	# j = 1
	my $m_max = 'AA' ;
	foreach my $m (keys %{$alpha_ref}) {
		if (($alpha_ref->{$m} * $Bjm{'1'}{$m}) > ($alpha_ref->{$m_max} * $Bjm{'1'}{$m_max})) {
			$m_max = $m ;
		}
	}
	foreach my $m (keys %{$alpha_ref}) {
		$Dj{'1'} += ($alpha_ref->{$m} * $Bjm{'1'}{$m}) / ($alpha_ref->{$m_max} * $Bjm{'1'}{$m_max}) ;
	}
	$Dj{'1'} *= $alpha_ref->{$m_max} * $Bjm{'1'}{$m_max} ;
	# j = 2
	if ($pi_ref->{2}!= 0) {
		foreach my $m (keys %{$alpha_ref}) {
			if (($alpha_ref->{$m} * $Bjm{'2'}{$m}) > ($alpha_ref->{$m_max} * $Bjm{'2'}{$m_max})) {
				$m_max = $m ;
			}
		}
		foreach my $m (keys %{$alpha_ref}) {
			$Dj{'2'} += ($alpha_ref->{$m} * $Bjm{'2'}{$m}) / ($alpha_ref->{$m_max} * $Bjm{'2'}{$m_max}) ;
		}
		$Dj{'2'} *= $alpha_ref->{$m_max} * $Bjm{'2'}{$m_max} ;
	}
	# j = 3
	if ($pi_ref->{3}!= 0) {
		foreach my $m (keys %{$alpha_ref}) {
			if (($alpha_ref->{$m} * $Bjm{'3'}{$m}) > ($alpha_ref->{$m_max} * $Bjm{'3'}{$m_max})) {
				$m_max = $m ;
			}
		}
		foreach my $m (keys %{$alpha_ref}) {
			$Dj{'3'} += ($alpha_ref->{$m} * $Bjm{'3'}{$m}) / ($alpha_ref->{$m_max} * $Bjm{'3'}{$m_max}) ;
		}
		$Dj{'3'} *= $alpha_ref->{$m_max} * $Bjm{'3'}{$m_max} ;
	}
	# Compute D = \sum_j \pi_j D_j
	my $j_max = 1 ;
	for (my $j=1 ; $j<4; $j++) {
		if (($pi_ref->{$j} * $Dj{$j}) > ($pi_ref->{$j_max} * $Dj{$j_max})) {
			$j_max = $j ;
		}
	}
	for (my $j=1 ; $j<4; $j++) {
		$D += ($pi_ref->{$j} * $Dj{$j}) / ($pi_ref->{$j_max} * $Dj{$j_max}) ;
	}
	$D *= $pi_ref->{$j_max} * $Dj{$j_max} ;

	# Compute TMG_k_t
	foreach my $m (keys %{$alpha_ref}) {
		# j = 1
		$TMG_k_t_1{$m} = $alpha_ref->{$m} * $Bjm{'1'}{$m} / $Dj{'1'} ;
		# j = 2
		if ($pi_ref->{2}!= 0) {
			$TMG_k_t_2{$m} = $alpha_ref->{$m} * $Bjm{'2'}{$m} / $Dj{'2'} ;
		}
		# j = 3
		if ($pi_ref->{3}!= 0) {
			$TMG_k_t_3{$m} = $alpha_ref->{$m} * $Bjm{'3'}{$m} / $Dj{'3'} ;
		}
	}

	# Compute TFG_k_t
	# j=1
	foreach my $n (keys %{$alpha_ref}) {
		$TFG_k_t_1{$n} = $alpha_ref->{$n} * $Ajn{'1'}{$n} / $Dj{'1'} ;
	}
	# j=2
	if ($pi_ref->{2}!= 0) {
		foreach my $n (keys %{$beta_2_ref}) {
			$TFG_k_t_2{$n} = $beta_2_ref->{$n} * $Ajn{'2'}{$n} / $Dj{'2'} ;
		}
	}
	# j=3
	if ($pi_ref->{3}!= 0) {
		foreach my $n (keys %{$beta_3_ref}) {
			$TFG_k_t_3{$n} = $beta_3_ref->{$n} * $Ajn{'3'}{$n} / $Dj{'3'} ;
		}
	}

	# Compute S_k_t
	$S_k_t{1} = $pi_ref->{'1'} * $Dj{'1'} / $D ;
	if ($pi_ref->{2}!= 0) {
		$S_k_t{2} = $pi_ref->{'2'} * $Dj{'2'} / $D ;
	}
	if ($pi_ref->{3}!= 0) {
		$S_k_t{3} = $pi_ref->{'3'} * $Dj{'3'} / $D ;
	}

	return (\%S_k_t, \%TMG_k_t_1, \%TMG_k_t_2, \%TMG_k_t_3, \%TFG_k_t_1, \%TFG_k_t_2, \%TFG_k_t_3) ;
}


#---------------
# Compute %Q_k_t
#---------------
sub compute_Q_k_t {

	# Recovering parameters
	my $Parameters_ref = shift ;
	my $S_k_t_ref = shift ;
	my $TMG_k_t_1_ref = shift ;
	my $TMG_k_t_2_ref = shift ;
	my $TMG_k_t_3_ref = shift ;
	my $TFG_k_t_1_ref = shift ;
	my $TFG_k_t_2_ref = shift ;
	my $TFG_k_t_3_ref = shift ;
	my $GE_k_t_1_ref = shift ;
	my $GE_k_t_2_ref = shift ;
	my $GE_k_t_3_ref = shift ;
	my $YGE_k_t_2_ref = shift ;
	my $alpha_ref = $Parameters_ref->{'alpha'} ;
	my $beta_2_ref = $Parameters_ref->{'beta_2'} ;
	my $beta_3_ref = $Parameters_ref->{'beta_3'} ;
	my $pi_ref = $Parameters_ref->{'pi'} ;
	my $heterogametic_genotypes_ref = $Parameters_ref->{'heterogametic_genotypes'} ;
	my $homogametic_genotypes_ref = $Parameters_ref->{'homogametic_genotypes'} ;
	my $lambda_1_ref = $Parameters_ref->{'lambda_1'} ;
	my $lambda_1_E_ref = $Parameters_ref->{'lambda_1_E'} ;
	my $lambda_2_ref = $Parameters_ref->{'lambda_2'} ;
	my $lambda_2_E_ref = $Parameters_ref->{'lambda_2_E'} ;
	my $lambda_2_p_ref = $Parameters_ref->{'lambda_2_p'} ;
	my $lambda_2_E_p_ref = $Parameters_ref->{'lambda_2_E_p'} ;
	my $lambda_3_ref = $Parameters_ref->{'lambda_3'} ;
	my $lambda_3_E_ref = $Parameters_ref->{'lambda_3_E'} ;
	my $mu_2_het_ref = $Parameters_ref->{'mu_2_het'} ;
	my $mu_2_p_het_ref = $Parameters_ref->{'mu_2_p_het'} ;
	my $mu_2_E_het_ref = $Parameters_ref->{'mu_2_E_het'} ;
	my $mu_2_E_p_het_ref = $Parameters_ref->{'mu_2_E_p_het'} ;
	my $mu_3_het_ref = $Parameters_ref->{'mu_3_het'} ;
	my $mu_3_E_het_ref = $Parameters_ref->{'mu_3_E_het'} ;
	my $mu_hom_ref = $Parameters_ref->{'mu_hom'} ;
	my $mu_E_hom_ref = $Parameters_ref->{'mu_E_hom'} ;
	my $individuals_ref = $Parameters_ref->{'individuals'} ;
	my $het_individuals_ref = $Parameters_ref->{'het_individuals'} ;

	my $Q_k_t = 0 ;

	# observed variables
	# j=1
	if ($S_k_t_ref->{1}!= 0) {
		foreach my $n (keys %{$alpha_ref}) {
			foreach my $m (keys %{$alpha_ref}) {
				my $sum = 0 ;
				# heterogametic progeny
				foreach my $i (keys %{$heterogametic_genotypes_ref}) {
					my $l = $heterogametic_genotypes_ref->{$i} ;
					if ($l ne 'NN') {
						if (($GE_k_t_1_ref->{$n.','.$m}->{$i} == 1)&&($lambda_1_ref->{$n.','.$m}->{$l} == 0)) {
							# 0*log(0) = 0
						} else {
							$sum += (1-$GE_k_t_1_ref->{$n.','.$m}->{$i})*log($lambda_1_ref->{$n.','.$m}->{$l}) ;
						}
						if (($GE_k_t_1_ref->{$n.','.$m}->{$i} == 0)&&($lambda_1_E_ref->{$n.','.$m}->{$l} == 0)) {
							# 0*log(0) = 0
						} else {
							$sum += $GE_k_t_1_ref->{$n.','.$m}->{$i}*log($lambda_1_E_ref->{$n.','.$m}->{$l}) ;
						}
					}
				}
				# homogametic progeny
				foreach my $i (keys %{$homogametic_genotypes_ref}) {
					my $l = $homogametic_genotypes_ref->{$i} ;
					if ($l ne 'NN') {
						if (($GE_k_t_1_ref->{$n.','.$m}->{$i} == 1)&&($lambda_1_ref->{$n.','.$m}->{$l} == 0)) {
							# 0*log(0) = 0
						} else {
							$sum +=(1-$GE_k_t_1_ref->{$n.','.$m}->{$i})*log($lambda_1_ref->{$n.','.$m}->{$l}) ;
						}
						if (($GE_k_t_1_ref->{$n.','.$m}->{$i} == 0)&&($lambda_1_E_ref->{$n.','.$m}->{$l} == 0)) {
							# 0*log(0) = 0
						} else {
							$sum += $GE_k_t_1_ref->{$n.','.$m}->{$i}*log($lambda_1_E_ref->{$n.','.$m}->{$l}) ;
						}
					}
				}
				# heterogametic parent
				my $l = $Parameters_ref->{'heterogametic_parent_genotype'} ;
				my $i = $Parameters_ref->{'heterogametic_parent_name'} ;
				if (($GE_k_t_1_ref->{$n.','.$m}->{$i} == 1)&&($mu_hom_ref->{$n}->{$l} == 0)) {
					# 0*log(0) = 0
				} else {
					$sum += (1-$GE_k_t_1_ref->{$n.','.$m}->{$i})*log($mu_hom_ref->{$n}->{$l}) ;
				}
				if (($GE_k_t_1_ref->{$n.','.$m}->{$i} == 0)&&($mu_E_hom_ref->{$n}->{$l} == 0)) {
					# 0*log(0) = 0
				} else {
					$sum += $GE_k_t_1_ref->{$n.','.$m}->{$i}*log($mu_E_hom_ref->{$n}->{$l}) ;
				}
				# homogametic parent
				$l = $Parameters_ref->{'homogametic_parent_genotype'} ;
				$i = $Parameters_ref->{'homogametic_parent_name'} ;
				if (($GE_k_t_1_ref->{$n.','.$m}->{$i} == 1)&&($mu_hom_ref->{$m}->{$l} == 0)) {
					# 0*log(0) = 0
				} else {
					$sum += (1-$GE_k_t_1_ref->{$n.','.$m}->{$i})*log($mu_hom_ref->{$m}->{$l}) ;
				}
				if (($GE_k_t_1_ref->{$n.','.$m}->{$i} == 0)&&($mu_E_hom_ref->{$m}->{$l} == 0)) {
					# 0*log(0) = 0
				} else {
					$sum += $GE_k_t_1_ref->{$n.','.$m}->{$i}*log($mu_E_hom_ref->{$m}->{$l});
				}
				# incrementing Q_k_t for current n and m
				$Q_k_t += $S_k_t_ref->{1} * $TFG_k_t_1_ref->{$n} * $TMG_k_t_1_ref->{$m} * $sum ;
			}
		}
	}
	# j = 2
	if ($S_k_t_ref->{2}!= 0) {
		foreach my $n (keys %{$beta_2_ref}) {
			foreach my $m (keys %{$alpha_ref}) {
				my $sum = 0 ;
				# heterogametic progeny
				foreach my $i (keys %{$heterogametic_genotypes_ref}) {
					my $l = $heterogametic_genotypes_ref->{$i} ;
					if ($l ne 'NN') {
						if ((($GE_k_t_2_ref->{$n.','.$m}->{$i} == 1)||($YGE_k_t_2_ref->{$n.','.$m}->{$i} == 1))&&($lambda_2_ref->{$n.','.$m}->{'het,'.$l} == 0)) {
							# 0*log(0) = 0
						} else {
							$sum += (1-$GE_k_t_2_ref->{$n.','.$m}->{$i})*(1-$YGE_k_t_2_ref->{$n.','.$m}->{$i})*log($lambda_2_ref->{$n.','.$m}->{'het,'.$l}) ;
						}
						if ((($GE_k_t_2_ref->{$n.','.$m}->{$i} == 0)||($YGE_k_t_2_ref->{$n.','.$m}->{$i} == 1))&&($lambda_2_E_ref->{$n.','.$m}->{'het,'.$l} == 0)) {
							# 0*log(0) = 0
						} else {
							$sum += $GE_k_t_2_ref->{$n.','.$m}->{$i}*(1-$YGE_k_t_2_ref->{$n.','.$m}->{$i})*log($lambda_2_E_ref->{$n.','.$m}->{'het,'.$l}) ;
						}
						if ((($GE_k_t_2_ref->{$n.','.$m}->{$i} == 1)||($YGE_k_t_2_ref->{$n.','.$m}->{$i} == 0))&&($lambda_2_p_ref->{$n.','.$m}->{'het,'.$l} == 0)) {
							# 0*log(0) = 0
						} else {
							$sum += (1-$GE_k_t_2_ref->{$n.','.$m}->{$i})*$YGE_k_t_2_ref->{$n.','.$m}->{$i}*log($lambda_2_p_ref->{$n.','.$m}->{'het,'.$l}) ;
						}
						if ((($GE_k_t_2_ref->{$n.','.$m}->{$i} == 0)||($YGE_k_t_2_ref->{$n.','.$m}->{$i} == 0))&&($lambda_2_E_p_ref->{$n.','.$m}->{'het,'.$l} == 0)) {
							# 0*log(0) = 0
						} else {
							$sum += $GE_k_t_2_ref->{$n.','.$m}->{$i}*$YGE_k_t_2_ref->{$n.','.$m}->{$i}*log($lambda_2_E_p_ref->{$n.','.$m}->{'het,'.$l}) ;
						}
					}
				}
				# homogametic progeny
				foreach my $i (keys %{$homogametic_genotypes_ref}) {
					my $l = $homogametic_genotypes_ref->{$i} ;
					if ($l ne 'NN') {
						if (($GE_k_t_2_ref->{$n.','.$m}->{$i} == 1)&&($lambda_2_ref->{$n.','.$m}->{'hom,'.$l} == 0)) {
							# 0*log(0) = 0
						} else {
							$sum += (1-$GE_k_t_2_ref->{$n.','.$m}->{$i})*log($lambda_2_ref->{$n.','.$m}->{'hom,'.$l}) ;
						}
						if (($GE_k_t_2_ref->{$n.','.$m}->{$i} == 0)&&($lambda_2_E_ref->{$n.','.$m}->{'hom,'.$l} == 0)) {
							# 0*log(0) = 0
						} else {
							$sum += $GE_k_t_2_ref->{$n.','.$m}->{$i}*log($lambda_2_E_ref->{$n.','.$m}->{'hom,'.$l}) ;
						}
					}
				}
				# heterogametic parent
				my $l = $Parameters_ref->{'heterogametic_parent_genotype'} ;
				my $i = $Parameters_ref->{'heterogametic_parent_name'} ;
				if ((($GE_k_t_2_ref->{$n.','.$m}->{$i} == 1)||($YGE_k_t_2_ref->{$n.','.$m}->{$i} == 1))&&($mu_2_het_ref->{$n}->{$l} == 0)) {
					# 0*log(0) = 0
				} else {
					$sum += (1-$GE_k_t_2_ref->{$n.','.$m}->{$i})*(1-$YGE_k_t_2_ref->{$n.','.$m}->{$i})*log($mu_2_het_ref->{$n}->{$l}) ;
				}
				if ((($GE_k_t_2_ref->{$n.','.$m}->{$i} == 0)||($YGE_k_t_2_ref->{$n.','.$m}->{$i} == 1))&&($mu_2_E_het_ref->{$n}->{$l} == 0)) {
					# 0*log(0) = 0
				} else {
					$sum += $GE_k_t_2_ref->{$n.','.$m}->{$i}*(1-$YGE_k_t_2_ref->{$n.','.$m}->{$i})*log($mu_2_E_het_ref->{$n}->{$l}) ;
				}
				if ((($GE_k_t_2_ref->{$n.','.$m}->{$i} == 1)||($YGE_k_t_2_ref->{$n.','.$m}->{$i} == 0))&&($mu_2_p_het_ref->{$n}->{$l} == 0)) {
					# 0*log(0) = 0
				} else {
					$sum += (1-$GE_k_t_2_ref->{$n.','.$m}->{$i})*$YGE_k_t_2_ref->{$n.','.$m}->{$i}*log($mu_2_p_het_ref->{$n}->{$l}) ;
				}
				if ((($GE_k_t_2_ref->{$n.','.$m}->{$i} == 0)||($YGE_k_t_2_ref->{$n.','.$m}->{$i} == 0))&&($mu_2_E_p_het_ref->{$n}->{$l} == 0)) {
					# 0*log(0) = 0
				} else {
					$sum += $GE_k_t_2_ref->{$n.','.$m}->{$i}*$YGE_k_t_2_ref->{$n.','.$m}->{$i}*log($mu_2_E_p_het_ref->{$n}->{$l}) ;
				}
				# homogametic parent
				$l = $Parameters_ref->{'homogametic_parent_genotype'} ;
				$i = $Parameters_ref->{'homogametic_parent_name'} ;
				if (($GE_k_t_2_ref->{$n.','.$m}->{$i} == 1)&&($mu_hom_ref->{$m}->{$l} == 0)) {
					# 0*log(0) = 0
				} else {
					$sum += (1-$GE_k_t_2_ref->{$n.','.$m}->{$i})*log($mu_hom_ref->{$m}->{$l}) ;
				}
				if (($GE_k_t_2_ref->{$n.','.$m}->{$i} == 0)&&($mu_E_hom_ref->{$m}->{$l} == 0)) {
					# 0*log(0) = 0
				} else {
					$sum += $GE_k_t_2_ref->{$n.','.$m}->{$i}*log($mu_E_hom_ref->{$m}->{$l}) ;
				}
				# incrementing Q_k_t for current n and m
				$Q_k_t += $S_k_t_ref->{2} * $TFG_k_t_2_ref->{$n} * $TMG_k_t_2_ref->{$m} * $sum ;
			}
		}
	}
	# j=3
	if ($S_k_t_ref->{3}!= 0) {
		foreach my $n (keys %{$beta_3_ref}) {
			foreach my $m (keys %{$alpha_ref}) {
				my $sum = 0 ;
				# heterogametic progeny
				foreach my $i (keys %{$heterogametic_genotypes_ref}) {
					my $l = $heterogametic_genotypes_ref->{$i} ;
					if ($l ne 'NN') {
						if (($GE_k_t_3_ref->{$n.','.$m}->{$i} == 1)&&($lambda_3_ref->{$n.','.$m}->{'het,'.$l} == 0)) {
							# 0*log(0) = 0
						} else {
							$sum += (1-$GE_k_t_3_ref->{$n.','.$m}->{$i})*log($lambda_3_ref->{$n.','.$m}->{'het,'.$l}) ;
						}
						if (($GE_k_t_3_ref->{$n.','.$m}->{$i} == 0)&&($lambda_3_E_ref->{$n.','.$m}->{'het,'.$l} == 0)) {
							# 0*log(0) = 0
						} else {
							$sum += $GE_k_t_3_ref->{$n.','.$m}->{$i}*log($lambda_3_E_ref->{$n.','.$m}->{'het,'.$l}) ;
						}
					}
				}
				# homogametic progeny
				foreach my $i (keys %{$homogametic_genotypes_ref}) {
					my $l = $homogametic_genotypes_ref->{$i} ;
					if ($l ne 'NN') {
						if (($GE_k_t_3_ref->{$n.','.$m}->{$i} == 1)&&($lambda_3_ref->{$n.','.$m}->{'hom,'.$l} == 0)) {
							# 0*log(0) = 0
						} else {
							$sum += (1-$GE_k_t_3_ref->{$n.','.$m}->{$i})*log($lambda_3_ref->{$n.','.$m}->{'hom,'.$l}) ;
						}
						if (($GE_k_t_3_ref->{$n.','.$m}->{$i} == 0)&&($lambda_3_E_ref->{$n.','.$m}->{'hom,'.$l} == 0)) {
							# 0*log(0) = 0
						} else {
							$sum += $GE_k_t_3_ref->{$n.','.$m}->{$i}*log($lambda_3_E_ref->{$n.','.$m}->{'hom,'.$l}) ;
						}
					}
				}
				# heterogametic parent
				my $l = $Parameters_ref->{'heterogametic_parent_genotype'} ;
				my $i = $Parameters_ref->{'heterogametic_parent_name'} ;
				if (($GE_k_t_3_ref->{$n.','.$m}->{$i} == 1)&&($mu_3_het_ref->{$n}->{$l} == 0)) {
					# 0*log(0) = 0
				} else {
					$sum += (1-$GE_k_t_3_ref->{$n.','.$m}->{$i})*log($mu_3_het_ref->{$n}->{$l}) ;
				}
				if (($GE_k_t_3_ref->{$n.','.$m}->{$i} == 0)&&($mu_3_E_het_ref->{$n}->{$l} == 0)) {
					# 0*log(0) = 0
				} else {
					$sum += $GE_k_t_3_ref->{$n.','.$m}->{$i}*log($mu_3_E_het_ref->{$n}->{$l}) ;
				}
				# homogametic parent
				$l = $Parameters_ref->{'homogametic_parent_genotype'} ;
				$i = $Parameters_ref->{'homogametic_parent_name'} ;
				if (($GE_k_t_3_ref->{$n.','.$m}->{$i} == 1)&&($mu_hom_ref->{$m}->{$l} == 0)) {
					# 0*log(0) = 0
				} else {
					$sum += (1-$GE_k_t_3_ref->{$n.','.$m}->{$i})*log($mu_hom_ref->{$m}->{$l}) ;
				}
				if (($GE_k_t_3_ref->{$n.','.$m}->{$i} == 0)&&($mu_E_hom_ref->{$m}->{$l} == 0)) {
					# 0*log(0) = 0
				} else {
					$sum += $GE_k_t_3_ref->{$n.','.$m}->{$i}*log($mu_E_hom_ref->{$m}->{$l}) ;
				}
				# incrementing Q_k_t for current n and m
				$Q_k_t += $S_k_t_ref->{3} * $TFG_k_t_3_ref->{$n} * $TMG_k_t_3_ref->{$m} * $sum ;
			}
		}
	}

	# hidden variables
	foreach my $j (1,2,3) {
		if ($S_k_t_ref->{$j}!= 0) {
			$Q_k_t += $S_k_t_ref->{$j}*log($pi_ref->{$j}) ;
		}
	}
	foreach my $m (keys %{$alpha_ref}) {
		if ($alpha_ref->{$m}!= 0) {
			if ($S_k_t_ref->{1}!= 0) {
				$Q_k_t += $S_k_t_ref->{1} * $TMG_k_t_1_ref->{$m}*log($alpha_ref->{$m}) ;
			}
			if ($S_k_t_ref->{2}!= 0) {
				$Q_k_t += $S_k_t_ref->{2} * $TMG_k_t_2_ref->{$m}*log($alpha_ref->{$m}) ;
			}
			if ($S_k_t_ref->{3}!= 0) {
				$Q_k_t += $S_k_t_ref->{3} * $TMG_k_t_3_ref->{$m}*log($alpha_ref->{$m}) ;
			}
		}
	}
	if ($S_k_t_ref->{1}!= 0) {
		foreach my $n (keys %{$TFG_k_t_1_ref}) {
			if ($alpha_ref->{$n}!= 0) {
				$Q_k_t += $S_k_t_ref->{1}*$TFG_k_t_1_ref->{$n}*log($alpha_ref->{$n}) ;
			}
		}
	}
	if ($S_k_t_ref->{2}!= 0) {
		foreach my $n (keys %{$TFG_k_t_2_ref}) {
			if ($beta_2_ref->{$n}!= 0) {
				$Q_k_t += $S_k_t_ref->{2}*$TFG_k_t_2_ref->{$n}*log($beta_2_ref->{$n}) ;
			}
		}
	}
	if ($S_k_t_ref->{3}!= 0) {
		foreach my $n (keys %{$TFG_k_t_3_ref}) {
			if ($beta_3_ref->{$n}!= 0) {
				$Q_k_t += $S_k_t_ref->{3}*$TFG_k_t_3_ref->{$n}*log($beta_3_ref->{$n}) ;
			}
		}
	}
	foreach my $i (keys %{$individuals_ref}) {
		foreach my $m (keys %{$alpha_ref}) {
			if ($S_k_t_ref->{1}!= 0) {
				foreach my $n (keys %{$TFG_k_t_1_ref}) {
					$Q_k_t += $S_k_t_ref->{1} * $TMG_k_t_1_ref->{$m} * $TFG_k_t_1_ref->{$n} * ( $GE_k_t_1_ref->{$n.','.$m}->{$i} * log($Parameters_ref->{'E'}) + (1-$GE_k_t_1_ref->{$n.','.$m}->{$i})*log(1-$Parameters_ref->{'E'}) ) ;
				}
			}
			if ($S_k_t_ref->{2}!= 0) {
				foreach my $n (keys %{$TFG_k_t_2_ref}) {
					$Q_k_t += $S_k_t_ref->{2} * $TMG_k_t_2_ref->{$m} * $TFG_k_t_2_ref->{$n} * ( $GE_k_t_2_ref->{$n.','.$m}->{$i} * log($Parameters_ref->{'E'}) + (1-$GE_k_t_2_ref->{$n.','.$m}->{$i})*log(1-$Parameters_ref->{'E'}) ) ;
				}
			}
			if ($S_k_t_ref->{3}!= 0) {
				foreach my $n (keys %{$TFG_k_t_3_ref}) {
					$Q_k_t += $S_k_t_ref->{3} * $TMG_k_t_3_ref->{$m} * $TFG_k_t_3_ref->{$n} * ( $GE_k_t_3_ref->{$n.','.$m}->{$i} * log($Parameters_ref->{'E'}) + (1-$GE_k_t_3_ref->{$n.','.$m}->{$i})*log(1-$Parameters_ref->{'E'}) ) ;
				}
			}
		}
	}
	if ($S_k_t_ref->{2}!= 0) {
		foreach my $i (keys %{$het_individuals_ref}) {
			foreach my $m (keys %{$alpha_ref}) {
				foreach my $n (keys %{$TFG_k_t_2_ref}) {
					$Q_k_t += $S_k_t_ref->{2} * $TMG_k_t_2_ref->{$m} * $TFG_k_t_2_ref->{$n} * ( $YGE_k_t_2_ref->{$n.','.$m}->{$i} * log($Parameters_ref->{'p'}) + (1-$YGE_k_t_2_ref->{$n.','.$m}->{$i})*log(1-$Parameters_ref->{'p'}) ) ;
				}
			}
		}
	}

	return $Q_k_t;
}


#---------------
# Compute %H_k_t
#---------------
sub compute_H_k_t {

	# Recovering parameters
	my $Parameters_ref = shift ;
	my $S_k_t_ref = shift ;
	my $TMG_k_t_1_ref = shift ;
	my $TMG_k_t_2_ref = shift ;
	my $TMG_k_t_3_ref = shift ;
	my $TFG_k_t_1_ref = shift ;
	my $TFG_k_t_2_ref = shift ;
	my $TFG_k_t_3_ref = shift ;
	my $GE_k_t_1_ref = shift ;
	my $GE_k_t_2_ref = shift ;
	my $GE_k_t_3_ref = shift ;
	my $YGE_k_t_2_ref = shift ;
	my $pi_ref = $Parameters_ref->{'pi'} ;
	my $individuals_ref = $Parameters_ref->{'individuals'} ;
	my $het_individuals_ref = $Parameters_ref->{'het_individuals'} ;

	my $H_k_t = 0 ;

	foreach my $j (1,2,3) {
		if ($S_k_t_ref->{$j}!= 0) {
			if (($S_k_t_ref->{$j}!= 1)&&($S_k_t_ref->{$j}!= 0)) {
				$H_k_t += $S_k_t_ref->{$j}*log($S_k_t_ref->{$j}) ;
			}
		}
	}
	foreach my $m (keys %{$TMG_k_t_1_ref}) {
		if ($S_k_t_ref->{1}!= 0) {
			if (($TMG_k_t_1_ref->{$m}!= 1)&&($TMG_k_t_1_ref->{$m}!= 0)) {
				$H_k_t += $S_k_t_ref->{1} * $TMG_k_t_1_ref->{$m}*log($TMG_k_t_1_ref->{$m}) ;
			}
		}
		if ($S_k_t_ref->{2}!= 0) {
			if (($TMG_k_t_2_ref->{$m}!= 1)&&($TMG_k_t_2_ref->{$m}!= 0)) {
				$H_k_t += $S_k_t_ref->{2} * $TMG_k_t_2_ref->{$m}*log($TMG_k_t_2_ref->{$m}) ;
			}
		}
		if ($S_k_t_ref->{3}!= 0) {
			if (($TMG_k_t_3_ref->{$m}!= 1)&&($TMG_k_t_3_ref->{$m}!= 0)) {
				$H_k_t += $S_k_t_ref->{3} * $TMG_k_t_3_ref->{$m}*log($TMG_k_t_3_ref->{$m}) ;
			}
		}
	}
	if ($S_k_t_ref->{1}!= 0) {
		foreach my $n (keys %{$TFG_k_t_1_ref}) {
			if (($TFG_k_t_1_ref->{$n}!= 1)&&($TFG_k_t_1_ref->{$n}!= 0)) {
				$H_k_t += $S_k_t_ref->{1}*$TFG_k_t_1_ref->{$n}*log($TFG_k_t_1_ref->{$n}) ;
			}
		}
	}
	if ($S_k_t_ref->{2}!= 0) {
		foreach my $n (keys %{$TFG_k_t_2_ref}) {
			if (($TFG_k_t_2_ref->{$n}!= 1)&&($TFG_k_t_2_ref->{$n}!= 0)) {
				$H_k_t += $S_k_t_ref->{2}*$TFG_k_t_2_ref->{$n}*log($TFG_k_t_2_ref->{$n}) ;
			}
		}
	}
	if ($S_k_t_ref->{3}!= 0) {
		foreach my $n (keys %{$TFG_k_t_3_ref}) {
			if (($TFG_k_t_3_ref->{$n}!= 1)&&($TFG_k_t_3_ref->{$n}!= 0)) {
				$H_k_t += $S_k_t_ref->{3}*$TFG_k_t_3_ref->{$n}*log($TFG_k_t_3_ref->{$n}) ;
			}
		}
	}
	foreach my $i (keys %{$individuals_ref}) {
		foreach my $m (keys %{$TMG_k_t_1_ref}) {
			if ($S_k_t_ref->{1}!= 0) {
				foreach my $n (keys %{$TFG_k_t_1_ref}) {
					if ($GE_k_t_1_ref->{$n.','.$m}->{$i}!= 0) {
						$H_k_t += $S_k_t_ref->{1} * $TMG_k_t_1_ref->{$m} * $TFG_k_t_1_ref->{$n} * $GE_k_t_1_ref->{$n.','.$m}->{$i}*log($GE_k_t_1_ref->{$n.','.$m}->{$i}) ;
					}
					if ($GE_k_t_1_ref->{$n.','.$m}->{$i}!= 1) {
						$H_k_t += $S_k_t_ref->{1} * $TMG_k_t_1_ref->{$m} * $TFG_k_t_1_ref->{$n} * (1-$GE_k_t_1_ref->{$n.','.$m}->{$i})*log(1-$GE_k_t_1_ref->{$n.','.$m}->{$i}) ;
					}
				}
			}
			if ($S_k_t_ref->{2}!= 0) {
				foreach my $n (keys %{$TFG_k_t_2_ref}) {
					if ($GE_k_t_2_ref->{$n.','.$m}->{$i}!= 0) {
						$H_k_t += $S_k_t_ref->{2} * $TMG_k_t_2_ref->{$m} * $TFG_k_t_2_ref->{$n} * $GE_k_t_2_ref->{$n.','.$m}->{$i}*log($GE_k_t_2_ref->{$n.','.$m}->{$i}) ;
					}
					if ($GE_k_t_2_ref->{$n.','.$m}->{$i}!= 1) {
						$H_k_t += $S_k_t_ref->{2} * $TMG_k_t_2_ref->{$m} * $TFG_k_t_2_ref->{$n} * (1-$GE_k_t_2_ref->{$n.','.$m}->{$i})*log(1-$GE_k_t_2_ref->{$n.','.$m}->{$i}) ;
					}
				}
			}
			if ($S_k_t_ref->{3}!= 0) {
				foreach my $n (keys %{$TFG_k_t_3_ref}) {
					if ($GE_k_t_3_ref->{$n.','.$m}->{$i}!= 0) {
						$H_k_t += $S_k_t_ref->{3} * $TMG_k_t_3_ref->{$m} * $TFG_k_t_3_ref->{$n} * $GE_k_t_3_ref->{$n.','.$m}->{$i}*log($GE_k_t_3_ref->{$n.','.$m}->{$i}) ;
					}
					if ($GE_k_t_3_ref->{$n.','.$m}->{$i}!= 1) {
						$H_k_t += $S_k_t_ref->{3} * $TMG_k_t_3_ref->{$m} * $TFG_k_t_3_ref->{$n} * (1-$GE_k_t_3_ref->{$n.','.$m}->{$i})*log(1-$GE_k_t_3_ref->{$n.','.$m}->{$i}) ;
					}
				}
			}
		}
	}
	if ($S_k_t_ref->{2}!= 0) {
		foreach my $i (keys %{$het_individuals_ref}) {
			foreach my $m (keys %{$TMG_k_t_1_ref}) {
				foreach my $n (keys %{$TFG_k_t_2_ref}) {
					if (($YGE_k_t_2_ref->{$n.','.$m}->{$i}!= 1)&&($YGE_k_t_2_ref->{$n.','.$m}->{$i}!= 0)) {
						$H_k_t += $S_k_t_ref->{2} * $TMG_k_t_2_ref->{$m} * $TFG_k_t_2_ref->{$n} * ( $YGE_k_t_2_ref->{$n.','.$m}->{$i}*log($YGE_k_t_2_ref->{$n.','.$m}->{$i}) + (1-$YGE_k_t_2_ref->{$n.','.$m}->{$i})*log(1-$YGE_k_t_2_ref->{$n.','.$m}->{$i}) ) ;
					}
				}
			}
		}
	}

	return $H_k_t;
}


#----------------------------
# prepare sex-linked details
#----------------------------
sub prepare_sex_linked_details {

	# Recovering parameters
	my ($Parameters_ref, $S_k_t_ref, $higher_2_type, $higher_3_type, $higher_TFG_k_t_2, $higher_TMG_k_t_2, $higher_TFG_k_t_3, $higher_TMG_k_t_3, $homogametic_expression_ref, $heterogametic_expression_ref, $homogametic_genotypes_ref, $heterogametic_genotypes_ref, $number_error_1, $number_error_2, $number_error_3, $number_Y_error) ;
	($Parameters_ref, $S_k_t_ref, $higher_2_type, $higher_3_type, $higher_TFG_k_t_2, $higher_TMG_k_t_2, $higher_TFG_k_t_3, $higher_TMG_k_t_3, $homogametic_expression_ref, $heterogametic_expression_ref, $homogametic_genotypes_ref, $heterogametic_genotypes_ref, $number_error_1, $number_error_2, $number_error_3, $number_Y_error) = @_ ;
	my %heterogametic_expression = %{$heterogametic_expression_ref} ;
	my %homogametic_expression = %{$homogametic_expression_ref} ;


	if (($S_k_t_ref->{2} > $S_k_t_ref->{1})&&($S_k_t_ref->{2} > $S_k_t_ref->{3})&&($higher_2_type ne 'monomorphic')&&($higher_2_type ne 'not_informative')) {
		# XY SNP
		my $Y = substr($higher_TFG_k_t_2, 1, 1) ;
		my $X1 = substr($higher_TFG_k_t_2, 0, 1) ;
		my $X3 = substr($higher_TMG_k_t_2, 1, 1) ;
		my $X2 = substr($higher_TMG_k_t_2, 0, 1) ;
		my $X2_expr = retrieve_base_expression($X2, $Parameters_ref->{'homogametic_parent_expression'}) ;
		my $X3_expr = retrieve_base_expression($X3, $Parameters_ref->{'homogametic_parent_expression'}) ;
		my $X2_plus_X3_expr = $X3_expr + $X2_expr ;
		my $Y_expr = retrieve_base_expression($Y, $Parameters_ref->{'heterogametic_parent_expression'}) ;
		my $X1_expr = retrieve_base_expression($X1, $Parameters_ref->{'heterogametic_parent_expression'}) ;
		my $Y_plus_X1_expr = $Y_expr + $X1_expr ;
		if ($X2 ne $X3) {
			# 2 possible genotypes in each type of progeny
			$Parameters_ref->{'detail_sex-linked_SNPs'} .= $Parameters_ref->{'contig_name'}."\t".$Parameters_ref->{'current_contig_position'}."\t".$S_k_t_ref->{2}."\t".$higher_2_type."\t".$number_error_1."\t".$number_error_2."\t".$number_error_3."\t".$number_Y_error."\t".$Parameters_ref->{'number_individuals_with_aberrant_reads'}."\t".$X2."\t".$X2_expr."\t".$X3."\t".$X3_expr."\t".'NA'."\t".$X2_plus_X3_expr ;
			if ($Y ne $X1) {
				$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".$X1."\t".$X1_expr."\t".$Y."\t".$Y_expr."\t".'NA'."\t".$Y_plus_X1_expr ;
			} else {
				$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".'NA'."\t".'NA'."\t".'NA'."\t".'NA'."\t".$X1."\t".$X1_expr ;
			}
			foreach my $hom_ind (@homogametic) {
				if ($homogametic_genotypes_ref->{$hom_ind} eq "NN") {
					# no genotype for this individual
					$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t"."N"."\t"."NA"."\t"."N"."\t"."NA"."\t"."N"."\t"."NA" ;
				} else {
					$X1_expr = retrieve_base_expression($X1, $homogametic_expression{$hom_ind}) ;
					$X2_expr = retrieve_base_expression($X2, $homogametic_expression{$hom_ind}) ;
					$X3_expr = retrieve_base_expression($X3, $homogametic_expression{$hom_ind}) ;
					my $X1_plus_X2_expr = $X1_expr + $X2_expr ;
					my $X1_plus_X3_expr = $X1_expr + $X3_expr ;
					if (($homogametic_genotypes_ref->{$hom_ind} eq "$X1$X2")||($homogametic_genotypes_ref->{$hom_ind} eq "$X2$X1")) {
						if ($X1 ne $X2) {
							$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".$X1."\t".$X1_expr."\t".$X2."\t".$X2_expr."\t".'NA'."\t".$X1_plus_X2_expr ;
						} else {
							$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".'NA'."\t".'NA'."\t".'NA'."\t".'NA'."\t".$X2."\t".$X2_expr ;
						}
					} elsif (($homogametic_genotypes_ref->{$hom_ind} eq "$X1$X3")||($homogametic_genotypes_ref->{$hom_ind} eq "$X3$X1")) {
						if ($X1 ne $X3) {
							$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".$X1."\t".$X1_expr."\t".$X3."\t".$X3_expr."\t".'NA'."\t".$X1_plus_X3_expr ;
						} else {
							$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".'NA'."\t".'NA'."\t".'NA'."\t".'NA'."\t".$X1."\t".$X1_expr ;
						}
					} else {
						# a genotyping error occured
						my $total_expr = retrieve_total_expression($homogametic_expression{$hom_ind}) ;
						if (($X1_expr / $total_expr) > 0.02) {
							if ((($X2_expr / $total_expr) > 0.02)&&(($X3_expr / $total_expr) > 0.02)) {
								# both X3 and X2 alleles expressed in individual, this is imposssible!
								$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t"."?"."\t"."NA"."\t"."?"."\t"."NA"."\t"."?"."\t"."NA" ;
							} elsif (($X2_expr / $total_expr) > 0.02) {
								# X1X2 genotype
								if ($X1 ne $X2) {
									$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".$X1.'*'."\t".$X1_expr."\t".$X2.'*'."\t".$X2_expr."\t".'NA'."\t".$X1_plus_X2_expr ;
								} else {
									$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".'NA'."\t".'NA'."\t".'NA'."\t".'NA'."\t".$X2.'*'."\t".$X2_expr ;
								}
							} elsif (($X3_expr / $total_expr) > 0.02) {
								# X1X3 genotype
								if ($X1 ne $X3) {
									$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".$X1.'*'."\t".$X1_expr."\t".$X3.'*'."\t".$X3_expr."\t".'NA'."\t".$X1_plus_X3_expr ;
								} else {
									$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".'NA'."\t".'NA'."\t".'NA'."\t".'NA'."\t".$X1.'*'."\t".$X1_expr ;
								}
							} else {
								# neither X3 nor X2 alleles expressed in individual, this is imposssible!
								$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t"."?"."\t"."NA"."\t"."?"."\t"."NA"."\t"."?"."\t"."NA" ;
							} 
						} else {
							# no X1 allele expressed in individual, this is imposssible!
							$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t"."?"."\t"."NA"."\t"."?"."\t"."NA"."\t"."?"."\t"."NA" ;
						}
					}
				}
			}
			foreach my $het_ind (@heterogametic) {
				if ($heterogametic_genotypes_ref->{$het_ind} eq "NN") {
					# no genotype for this individual
					$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t"."N"."\t"."NA"."\t"."N"."\t"."NA"."\t"."N"."\t"."NA" ;
				} else {
					$Y_expr = retrieve_base_expression($Y, $heterogametic_expression{$het_ind}) ;
					$X2_expr = retrieve_base_expression($X2, $heterogametic_expression{$het_ind}) ;
					$X3_expr = retrieve_base_expression($X3, $heterogametic_expression{$het_ind}) ;
					my $Y_plus_X2_expr = $Y_expr + $X2_expr ;
					my $Y_plus_X3_expr = $Y_expr + $X3_expr ;
					if (($heterogametic_genotypes_ref->{$het_ind} eq "$X2$Y")||($heterogametic_genotypes_ref->{$het_ind} eq "$Y$X2")||($heterogametic_genotypes_ref->{$het_ind} eq "$X2$X2")) {
						if ($Y ne $X2) {
							$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".$X2."\t".$X2_expr."\t".$Y."\t".$Y_expr."\t".'NA'."\t".$Y_plus_X2_expr ;
						} else {
							$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".'NA'."\t".'NA'."\t".'NA'."\t".'NA'."\t".$X2."\t".$X2_expr ;
						} 
					} elsif (($heterogametic_genotypes_ref->{$het_ind} eq "$X3$Y")||($heterogametic_genotypes_ref->{$het_ind} eq "$Y$X3")||($heterogametic_genotypes_ref->{$het_ind} eq "$X3$X3")) {
						if ($Y ne $X3) {
							$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".$X3."\t".$X3_expr."\t".$Y."\t".$Y_expr."\t".'NA'."\t".$Y_plus_X3_expr ;
						} else {
							$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".'NA'."\t".'NA'."\t".'NA'."\t".'NA'."\t".$X3."\t".$X3_expr ;
						} 
					} else {
						# a genotyping error occured
						my $total_expr = retrieve_total_expression($heterogametic_expression{$het_ind}) ;
						if ((($X2_expr / $total_expr) > 0.02)&&(($X3_expr / $total_expr) > 0.02)) {
							# both X2 and X3 alleles expressed in individual, this is imposssible!
							$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t"."?"."\t"."NA"."\t"."?"."\t"."NA"."\t"."?"."\t"."NA" ;
						} elsif (($X2_expr / $total_expr) > 0.02) {
							# YX2 genotype
							if ($Y ne $X2) {
								$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".$X2.'*'."\t".$X2_expr."\t".$Y.'*'."\t".$Y_expr."\t".'NA'."\t".$Y_plus_X2_expr ;
							} else {
								$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".'NA'."\t".'NA'."\t".'NA'."\t".'NA'."\t".$X2.'*'."\t".$X2_expr ;
							} 
						} elsif (($X3_expr / $total_expr) > 0.02) {
							# YX3 genotype
							if ($Y ne $X2) {
								$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".$X3.'*'."\t".$X3_expr."\t".$Y.'*'."\t".$Y_expr."\t".'NA'."\t".$Y_plus_X3_expr ;
							} else {
								$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".'NA'."\t".'NA'."\t".'NA'."\t".'NA'."\t".$X3.'*'."\t".$X3_expr ;
							} 
						} else {
							# neither X3 nor X2 alleles expressed in individual, this is imposssible!
							$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t"."?"."\t"."NA"."\t"."?"."\t"."NA"."\t"."?"."\t"."NA" ;
						} 
					}
				}
			}
		} else {
			# 1 possible genotype in each type of progeny
			$Parameters_ref->{'detail_sex-linked_SNPs'} .= $Parameters_ref->{'contig_name'}."\t".$Parameters_ref->{'current_contig_position'}."\t".$S_k_t_ref->{2}."\t".$higher_2_type."\t".$number_error_1."\t".$number_error_2."\t".$number_error_3."\t".$number_Y_error."\t".$Parameters_ref->{'number_individuals_with_aberrant_reads'}."\t".'NA'."\t".'NA'."\t".'NA'."\t".'NA'."\t".$X2."\t".$X2_expr ;
			if ($Y ne $X1) {
				$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".$X1."\t".$X1_expr."\t".$Y."\t".$Y_expr."\t".'NA'."\t".$Y_plus_X1_expr ;
			} else {
				$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".'NA'."\t".'NA'."\t".'NA'."\t".'NA'."\t".$X1."\t".$X1_expr ;
			}
			foreach my $hom_ind (@homogametic) {
				if ($homogametic_genotypes_ref->{$hom_ind} eq "NN") {
					# no genotype for this individual
					$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t"."N"."\t"."NA"."\t"."N"."\t"."NA"."\t"."N"."\t"."NA" ;
				} else {
					$X1_expr = retrieve_base_expression($X1, $homogametic_expression{$hom_ind}) ;
					$X2_expr = retrieve_base_expression($X2, $homogametic_expression{$hom_ind}) ;
					my $X1_plus_X2_expr = $X1_expr + $X2_expr ;
					if (($homogametic_genotypes_ref->{$hom_ind} eq "$X1$X2")||($homogametic_genotypes_ref->{$hom_ind} eq "$X2$X1")) {
						if ($X1 ne $X2) {
							$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".$X1."\t".$X1_expr."\t".$X2."\t".$X2_expr."\t".'NA'."\t".$X1_plus_X2_expr ;
						} else {
							$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".'NA'."\t".'NA'."\t".'NA'."\t".'NA'."\t".$X2."\t".$X2_expr ;
						}
					} else {
						# a genotyping error occured
						my $total_expr = retrieve_total_expression($homogametic_expression{$hom_ind}) ;
						if ((($X1_expr / $total_expr) > 0.02)&&(($X2_expr / $total_expr) > 0.02)) {
							# X1X2 genotype
							if ($X1 ne $X2) {
								$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".$X1.'*'."\t".$X1_expr."\t".$X2.'*'."\t".$X2_expr."\t".'NA'."\t".$X1_plus_X2_expr ;
							} else {
								$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".'NA'."\t".'NA'."\t".'NA'."\t".'NA'."\t".$X2.'*'."\t".$X2_expr ;
							}
						} else {
							# either X1 or X2 allele not expressed in individual, this is imposssible!
							$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t"."?"."\t"."NA"."\t"."?"."\t"."NA"."\t"."?"."\t"."NA" ;
						}
					}
				}
			}
			foreach my $het_ind (@heterogametic) {
				if ($heterogametic_genotypes_ref->{$het_ind} eq "NN") {
					# no genotype for this individual
					$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t"."N"."\t"."NA"."\t"."N"."\t"."NA"."\t"."N"."\t"."NA" ;
				} else {
					$Y_expr = retrieve_base_expression($Y, $heterogametic_expression{$het_ind}) ;
					$X2_expr = retrieve_base_expression($X2, $heterogametic_expression{$het_ind}) ;
					my $Y_plus_X2_expr = $Y_expr + $X2_expr ;
					if (($heterogametic_genotypes_ref->{$het_ind} eq "$X2$Y")||($heterogametic_genotypes_ref->{$het_ind} eq "$Y$X2")||($heterogametic_genotypes_ref->{$het_ind} eq "$X2$X2")) {
						if ($Y ne $X2) {
							$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".$X2."\t".$X2_expr."\t".$Y."\t".$Y_expr."\t".'NA'."\t".$Y_plus_X2_expr ;
						} else {
							$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".'NA'."\t".'NA'."\t".'NA'."\t".'NA'."\t".$X2."\t".$X2_expr ;
						} 
					} else {
						# a genotyping error occured
						my $total_expr = retrieve_total_expression($heterogametic_expression{$het_ind}) ;
						if (($X2_expr / $total_expr) > 0.02) {
							# YX2 genotype
							if ($Y ne $X2) {
								$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".$X2.'*'."\t".$X2_expr."\t".$Y.'*'."\t".$Y_expr."\t".'NA'."\t".$Y_plus_X2_expr ;
							} else {
								$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".'NA'."\t".'NA'."\t".'NA'."\t".'NA'."\t".$X2.'*'."\t".$X2_expr ;
							} 
						} else {
							# no X2 allele expressed in individual, this is imposssible!
							$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t"."?"."\t"."NA"."\t"."?"."\t"."NA"."\t"."?"."\t"."NA" ;
						} 
					}
				}
			}
		}
		$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\n" ;
	} elsif (($S_k_t_ref->{3} > $S_k_t_ref->{1})&&($S_k_t_ref->{3} > $S_k_t_ref->{2})&&($higher_3_type ne 'monomorphic')&&($higher_3_type ne 'not_informative')) {
		# X hemizygous SNP
		my $X1 = substr($higher_TFG_k_t_3, 0, 1) ;
		my $X3 = substr($higher_TMG_k_t_3, 1, 1) ;
		my $X2 = substr($higher_TMG_k_t_3, 0, 1) ;
		my $X2_expr = retrieve_base_expression($X2, $Parameters_ref->{'homogametic_parent_expression'}) ;
		my $X3_expr = retrieve_base_expression($X3, $Parameters_ref->{'homogametic_parent_expression'}) ;
		my $X2_plus_X3_expr = $X3_expr + $X2_expr ;
		my $X1_expr = retrieve_base_expression($X1, $Parameters_ref->{'heterogametic_parent_expression'}) ;
		if ($X2 ne $X3) {
			# 2 possible genotypes in each type of progeny
			$Parameters_ref->{'detail_sex-linked_SNPs'} .= $Parameters_ref->{'contig_name'}."\t".$Parameters_ref->{'current_contig_position'}."\t".$S_k_t_ref->{3}."\t".$higher_3_type."\t".$number_error_1."\t".$number_error_2."\t".$number_error_3."\t".$number_Y_error."\t".$Parameters_ref->{'number_individuals_with_aberrant_reads'}."\t".$X2."\t".$X2_expr."\t".$X3."\t".$X3_expr."\t".'NA'."\t".$X2_plus_X3_expr ;
			$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".$X1."\t".$X1_expr."\t".'N'."\t".'0'."\t".'NA'."\t".$X1_expr ;
			foreach my $hom_ind (@homogametic) {
				if ($homogametic_genotypes_ref->{$hom_ind} eq "NN") {
					# no genotype for this individual
					$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t"."N"."\t"."NA"."\t"."N"."\t"."NA"."\t"."N"."\t"."NA" ;
				} else {
					$X1_expr = retrieve_base_expression($X1, $homogametic_expression{$hom_ind}) ;
					$X2_expr = retrieve_base_expression($X2, $homogametic_expression{$hom_ind}) ;
					$X3_expr = retrieve_base_expression($X3, $homogametic_expression{$hom_ind}) ;
					my $X1_plus_X2_expr = $X1_expr + $X2_expr ;
					my $X1_plus_X3_expr = $X1_expr + $X3_expr ;
					if (($homogametic_genotypes_ref->{$hom_ind} eq "$X1$X2")||($homogametic_genotypes_ref->{$hom_ind} eq "$X2$X1")) {
						if ($X1 ne $X2) {
							$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".$X1."\t".$X1_expr."\t".$X2."\t".$X2_expr."\t".'NA'."\t".$X1_plus_X2_expr ;
						} else {
							$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".'NA'."\t".'NA'."\t".'NA'."\t".'NA'."\t".$X2."\t".$X2_expr ;
						}
					} elsif (($homogametic_genotypes_ref->{$hom_ind} eq "$X1$X3")||($homogametic_genotypes_ref->{$hom_ind} eq "$X3$X1")) {
						if ($X1 ne $X3) {
							$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".$X1."\t".$X1_expr."\t".$X3."\t".$X3_expr."\t".'NA'."\t".$X1_plus_X3_expr ;
						} else {
							$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".'NA'."\t".'NA'."\t".'NA'."\t".'NA'."\t".$X1."\t".$X1_expr ;
						}
					} else {
						# a genotyping error occured
						my $total_expr = retrieve_total_expression($homogametic_expression{$hom_ind}) ;
						if (($X1_expr / $total_expr) > 0.02) {
							if ((($X2_expr / $total_expr) > 0.02)&&(($X3_expr / $total_expr) > 0.02)) {
								# both X3 and X2 alleles expressed in individual, this is imposssible!
								$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t"."?"."\t"."NA"."\t"."?"."\t"."NA"."\t"."?"."\t"."NA" ;
							} elsif (($X2_expr / $total_expr) > 0.02) {
								# X1X2 genotype
								if ($X1 ne $X2) {
									$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".$X1.'*'."\t".$X1_expr."\t".$X2.'*'."\t".$X2_expr."\t".'NA'."\t".$X1_plus_X2_expr ;
								} else {
									$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".'NA'."\t".'NA'."\t".'NA'."\t".'NA'."\t".$X2.'*'."\t".$X2_expr ;
								}
							} elsif (($X3_expr / $total_expr) > 0.02) {
								# X1X3 genotype
								if ($X1 ne $X3) {
									$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".$X1.'*'."\t".$X1_expr."\t".$X3.'*'."\t".$X3_expr."\t".'NA'."\t".$X1_plus_X3_expr ;
								} else {
									$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".'NA'."\t".'NA'."\t".'NA'."\t".'NA'."\t".$X1.'*'."\t".$X1_expr ;
								}
							} else {
								# neither X3 nor X2 alleles expressed in individual, this is imposssible!
								$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t"."?"."\t"."NA"."\t"."?"."\t"."NA"."\t"."?"."\t"."NA" ;
							} 
						} else {
							# no X1 allele expressed in individual, this is imposssible!
							$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t"."?"."\t"."NA"."\t"."?"."\t"."NA"."\t"."?"."\t"."NA" ;
						}
					}
				}
			}
			foreach my $het_ind (@heterogametic) {
				if ($heterogametic_genotypes_ref->{$het_ind} eq "NN") {
					# no genotype for this individual
					$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t"."N"."\t"."NA"."\t"."N"."\t"."NA"."\t"."N"."\t"."NA" ;
				} else {
					$X2_expr = retrieve_base_expression($X2, $heterogametic_expression{$het_ind}) ;
					$X3_expr = retrieve_base_expression($X3, $heterogametic_expression{$het_ind}) ;
					if ($heterogametic_genotypes_ref->{$het_ind} eq "$X2$X2") {
						$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".$X2."\t".$X2_expr."\t".'NA'."\t".'0'."\t".'NA'."\t".$X2_expr ;
					} elsif ($heterogametic_genotypes_ref->{$het_ind} eq "$X3$X3") {
						$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".$X3."\t".$X3_expr."\t".'NA'."\t".'0'."\t".'NA'."\t".$X3_expr ;
					} else {
						# a genotyping error occured
						my $total_expr = retrieve_total_expression($heterogametic_expression{$het_ind}) ;
						if ((($X2_expr / $total_expr) > 0.02)&&(($X3_expr / $total_expr) > 0.02)) {
							# both X2 and X3 alleles expressed in individual, this is imposssible!
							$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t"."?"."\t"."NA"."\t"."?"."\t"."NA"."\t"."?"."\t"."NA" ;
						} elsif (($X2_expr / $total_expr) > 0.02) {
							# X2 genotype
							$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".$X2.'*'."\t".$X2_expr."\t".'NA'."\t".'0'."\t".'NA'."\t".$X2_expr ;
						} elsif (($X3_expr / $total_expr) > 0.02) {
							# X3 genotype
							$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".$X3.'*'."\t".$X3_expr."\t".'NA'."\t".'0'."\t".'NA'."\t".$X3_expr ;
						} else {
							# neither X3 nor X2 alleles expressed in individual, this is imposssible!
							$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t"."?"."\t"."NA"."\t"."?"."\t"."NA"."\t"."?"."\t"."NA" ;
						} 
					}
				}
			}
		} else {
			# 1 possible genotype in each type of progeny
			$Parameters_ref->{'detail_sex-linked_SNPs'} .= $Parameters_ref->{'contig_name'}."\t".$Parameters_ref->{'current_contig_position'}."\t".$S_k_t_ref->{3}."\t".$higher_3_type."\t".$number_error_1."\t".$number_error_2."\t".$number_error_3."\t".$number_Y_error."\t".$Parameters_ref->{'number_individuals_with_aberrant_reads'}."\t".'NA'."\t".'NA'."\t".'NA'."\t".'NA'."\t".$X2."\t".$X2_expr ;
			$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".$X1."\t".$X1_expr."\t".'N'."\t".'0'."\t".'NA'."\t".$X1_expr ;
			foreach my $hom_ind (@homogametic) {
				if ($homogametic_genotypes_ref->{$hom_ind} eq "NN") {
					# no genotype for this individual
					$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t"."N"."\t"."NA"."\t"."N"."\t"."NA"."\t"."N"."\t"."NA" ;
				} else {
					$X1_expr = retrieve_base_expression($X1, $homogametic_expression{$hom_ind}) ;
					$X2_expr = retrieve_base_expression($X2, $homogametic_expression{$hom_ind}) ;
					my $X1_plus_X2_expr = $X1_expr + $X2_expr ;
					if (($homogametic_genotypes_ref->{$hom_ind} eq "$X1$X2")||($homogametic_genotypes_ref->{$hom_ind} eq "$X2$X1")) {
						if ($X1 ne $X2) {
							$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".$X1."\t".$X1_expr."\t".$X2."\t".$X2_expr."\t".'NA'."\t".$X1_plus_X2_expr ;
						} else {
							$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".'NA'."\t".'NA'."\t".'NA'."\t".'NA'."\t".$X2."\t".$X2_expr ;
						}
					} else {
						# a genotyping error occured
						my $total_expr = retrieve_total_expression($homogametic_expression{$hom_ind}) ;
						if ((($X1_expr / $total_expr) > 0.02)&&(($X2_expr / $total_expr) > 0.02)) {
							# X1X2 genotype
							if ($X1 ne $X2) {
								$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".$X1.'*'."\t".$X1_expr."\t".$X2.'*'."\t".$X2_expr."\t".'NA'."\t".$X1_plus_X2_expr ;
							} else {
								$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".'NA'."\t".'NA'."\t".'NA'."\t".'NA'."\t".$X2.'*'."\t".$X2_expr ;
							}
						} else {
							# either X1 or X2 allele not expressed in individual, this is imposssible!
							$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t"."?"."\t"."NA"."\t"."?"."\t"."NA"."\t"."?"."\t"."NA" ;
						}
					}
				}
			}
			foreach my $het_ind (@heterogametic) {
				if ($heterogametic_genotypes_ref->{$het_ind} eq "NN") {
					# no genotype for this individual
					$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t"."N"."\t"."NA"."\t"."N"."\t"."NA"."\t"."N"."\t"."NA" ;
				} else {
					$X2_expr = retrieve_base_expression($X2, $heterogametic_expression{$het_ind}) ;
					if ($heterogametic_genotypes_ref->{$het_ind} eq "$X2$X2") {
						$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".$X2."\t".$X2_expr."\t".'NA'."\t".'0'."\t".'NA'."\t".$X2_expr ;
					} else {
						# a genotyping error occured
						my $total_expr = retrieve_total_expression($heterogametic_expression{$het_ind}) ;
						if (($X2_expr / $total_expr) > 0.02) {
							# YX2 genotype
							$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t".$X2.'*'."\t".$X2_expr."\t".'NA'."\t".'0'."\t".'NA'."\t".$X2_expr ;
						} else {
							# no X2 allele expressed in individual, this is imposssible!
							$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\t"."?"."\t"."NA"."\t"."?"."\t"."NA"."\t"."?"."\t"."NA" ;
						} 
					}
				}
			}
		}
		$Parameters_ref->{'detail_sex-linked_SNPs'} .= "\n" ;
	}
}


#----------------------------------
# Analyze all results from a contig
#----------------------------------
sub analyze_contig_results {
	
	# Recovering parameters
	my $Parameters_ref = shift ;

	# calcul nombre SNP sans erreur
	my $number_SNP_without_error = $Parameters_ref->{'autosomal_without_error'} + $Parameters_ref->{'XY_without_error'} + $Parameters_ref->{'hemizygous_without_error'} ;

	# bring sum of probabilities to 1
	if (($Parameters_ref->{'S_k_1'}!= 0)||($Parameters_ref->{'S_k_2'}!= 0)||($Parameters_ref->{'S_k_3'}!= 0)) {
		$Parameters_ref->{'S_k_1_end'} = $Parameters_ref->{'S_k_1'} / ( $Parameters_ref->{'S_k_1'} + $Parameters_ref->{'S_k_2'} + $Parameters_ref->{'S_k_3'} ) ;
		$Parameters_ref->{'S_k_2_end'} = $Parameters_ref->{'S_k_2'} / ( $Parameters_ref->{'S_k_1'} + $Parameters_ref->{'S_k_2'} + $Parameters_ref->{'S_k_3'} ) ;
		$Parameters_ref->{'S_k_3_end'} = $Parameters_ref->{'S_k_3'} / ( $Parameters_ref->{'S_k_1'} + $Parameters_ref->{'S_k_2'} + $Parameters_ref->{'S_k_3'} ) ;
	} else {
		$Parameters_ref->{'S_k_1_end'} = $Parameters_ref->{'S_k_1'} ;
		$Parameters_ref->{'S_k_2_end'} = $Parameters_ref->{'S_k_2'} ;
		$Parameters_ref->{'S_k_3_end'} = $Parameters_ref->{'S_k_3'} ;
	}

	# print output line for contig
	open(OUTPUT, '>>'.$Parameters_ref->{'output_file_name'}) ;
	print(OUTPUT $Parameters_ref->{'contig_name'}."\t".$Parameters_ref->{'S_k_1_end'}."\t".$Parameters_ref->{'S_k_2_end'}."\t".$Parameters_ref->{'S_k_3_end'}."\t") ;
	if (($Parameters_ref->{'S_k_1'} == 0)&&($Parameters_ref->{'S_k_2'} == 0)&&($Parameters_ref->{'S_k_3'} == 0)) {
		# not enough information (expression or individuals)
		print(OUTPUT 'lack-information'."\t") ;
	} else {
		if ($number_SNP_without_error > 0) {
			if (($Parameters_ref->{'S_k_1_end'}>$Parameters_ref->{'threshold'})&&($Parameters_ref->{'S_k_1_end'}>($Parameters_ref->{'S_k_2_end'}+$Parameters_ref->{'S_k_3_end'}))&&($Parameters_ref->{'autosomal_without_error'}>0)) {
				print(OUTPUT 'autosomal'."\t") ;
			} elsif ((($Parameters_ref->{'S_k_2_end'}+$Parameters_ref->{'S_k_3_end'})>$Parameters_ref->{'threshold'})&&($Parameters_ref->{'S_k_1_end'}<($Parameters_ref->{'S_k_2_end'}+$Parameters_ref->{'S_k_3_end'}))&&(($Parameters_ref->{'number_clean_XY_SNPs_without_error'}>0)||($Parameters_ref->{'number_clean_hemizygous_SNPs_without_error'}>0))) {
				print(OUTPUT 'sex-linked'."\t") ;
				if (defined($Parameters_ref->{'detail-sex-linked'})) {
					open(SEX_DETAIL, ">>".$Parameters_ref->{'sex-linked_detail_file_name'}) ;
					print(SEX_DETAIL $Parameters_ref->{'detail_sex-linked_SNPs'}) ;
					close(SEX_DETAIL) ;
				}
				if (defined($Parameters_ref->{'sequences'})) {
					open(SEQ, ">>".$Parameters_ref->{'sequences_file_name'}) ;
					print(SEQ '>'.$Parameters_ref->{'contig_name'}."_X1\n".$Parameters_ref->{'current_contig_sequence_X1'}."\n") ;
					if ($Parameters_ref->{'current_contig_sequence_X2'} ne $Parameters_ref->{'current_contig_sequence_X1'}) {
						print(SEQ '>'.$Parameters_ref->{'contig_name'}."_X2\n".$Parameters_ref->{'current_contig_sequence_X2'}."\n") ;
					}
					if (($Parameters_ref->{'current_contig_sequence_X3'} ne $Parameters_ref->{'current_contig_sequence_X1'})&&($Parameters_ref->{'current_contig_sequence_X3'} ne $Parameters_ref->{'current_contig_sequence_X2'})) {
						print(SEQ '>'.$Parameters_ref->{'contig_name'}."_X3\n".$Parameters_ref->{'current_contig_sequence_X3'}."\n") ;
					}
					if (($Parameters_ref->{'XY_without_error'}+$Parameters_ref->{'XY_with_error'}) > 0) {
						# the gene is not hemizygous
						print(SEQ '>'.$Parameters_ref->{'contig_name'}."_Y\n".$Parameters_ref->{'current_contig_sequence_Y'}."\n") ;
					}
					close(SEQ) ;
				}
			} else {
				print(OUTPUT 'lack-information'."\t") ;
			}
		} else {
			print(OUTPUT 'lack-information'."\t") ;
		}
	}
	print(OUTPUT $number_SNP_without_error."\t".$Parameters_ref->{'autosomal_without_error'}."\t".$Parameters_ref->{'autosomal_with_error'}."\t".$Parameters_ref->{'XY_without_error'}."\t".$Parameters_ref->{'XY_with_error'}."\t".$Parameters_ref->{'hemizygous_without_error'}."\t".$Parameters_ref->{'hemizygous_with_error'}) ;
	if (defined($Parameters_ref->{'ALR_file'})) {
		print(OUTPUT "\t".$Parameters_ref->{'number_clean_XY_SNPs_without_error'}."\t".$Parameters_ref->{'number_clean_hemizygous_SNPs_without_error'}."\n") ;
	} else {
		print(OUTPUT "\n") ;
	}
	close(OUTPUT) ;
	
	return 0;
}


#------------------------------
# Compute new parameters values
#------------------------------
sub Maximization_Step {
	
	# Recovering parameters
	my $Parameters_ref = shift ;
	my $pi_ref = $Parameters_ref->{'pi'} ;
	my $alpha_ref = $Parameters_ref->{'alpha'} ;
	my $beta_2_ref = $Parameters_ref->{'beta_2'} ;
	my $beta_3_ref = $Parameters_ref->{'beta_3'} ;
	my $sum_k_t_S_k_t_ref = $Parameters_ref->{'sum_k_t_S_k_t'} ;
	my $alpha_new_ref = $Parameters_ref->{'alpha_new'} ;
	my $beta_2_new_ref = $Parameters_ref->{'beta_2_new'} ;
	my $beta_3_new_ref = $Parameters_ref->{'beta_3_new'} ;
	my $TMG_1_counts_ref = $Parameters_ref->{'TMG_1_counts'} ;
	my $TMG_2_counts_ref = $Parameters_ref->{'TMG_2_counts'} ;
	my $TMG_3_counts_ref = $Parameters_ref->{'TMG_3_counts'} ;
	my $TFG_1_counts_ref = $Parameters_ref->{'TFG_1_counts'} ;
	my $TFG_2_counts_ref = $Parameters_ref->{'TFG_2_counts'} ;
	my $TFG_3_counts_ref = $Parameters_ref->{'TFG_3_counts'} ;
	my $S_counts_ref = $Parameters_ref->{'S_counts'} ;

	# registering old values of parameters and likelihoods in order to compare them with new values
	my %alpha_old = %{$alpha_ref} ;
	my %beta_2_old = %{$beta_2_ref} ;
	my %beta_3_old = %{$beta_3_ref} ;
	my %pi_old = %{$pi_ref} ;
	$Parameters_ref->{'pi_old'} = \%pi_old ;
	$Parameters_ref->{'alpha_old'} = \%alpha_old ;
	$Parameters_ref->{'beta_2_old'} = \%beta_2_old ;
	$Parameters_ref->{'beta_3_old'} = \%beta_3_old ;
	$Parameters_ref->{'p_old'} = $Parameters{'p'} ;
	$Parameters_ref->{'E_old'} = $Parameters{'E'} ;

	# compute number of free parameters
	my $number_free_param = 27 - $Parameters_ref->{'non_optimized_parameters'} ;

	# New pi
	if ($S_counts_ref->{1} < 1) {
		print("Error! Not a single position was attributed to autosomal segregation type\n") ;
		exit(1);
	} else {
		$pi_ref->{1} = $sum_k_t_S_k_t_ref->{1} / $Parameters_ref->{'sum_k_t'} ;
	}
	if ($S_counts_ref->{2} < 1) {
		# Not a single position was attributed to hemizygous segregation type
		$pi_ref->{2} = 0 ;
	} else {
		$pi_ref->{2} = $sum_k_t_S_k_t_ref->{2} / $Parameters_ref->{'sum_k_t'} ;
	}
	if ($S_counts_ref->{3} < 1) {
		# Not a single position was attributed to hemizygous segregation type
		$pi_ref->{3} = 0 ;
	} else {
		$pi_ref->{3} = $sum_k_t_S_k_t_ref->{3} / $Parameters_ref->{'sum_k_t'} ;
	}

	# New P parameter
	if ($S_counts_ref->{2} > 1) {
		$Parameters_ref->{'p'} = $Parameters_ref->{'p_new_numerator'} / $Parameters_ref->{'p_new_denominator'} ;
		# Contrainte sur p
		if ($Parameters_ref->{'p'} > 0.1) {
			$Parameters_ref->{'p'} = 0.1 ;
		}
	} else {
		$Parameters_ref->{'p'} = 0 ;
	}
	if ($Parameters_ref->{'p_fixed_to_zero'} eq "yes") {
		$Parameters_ref->{'p'} = 0 ;
	}

	# New E parameter
	$Parameters_ref->{'E'} = $Parameters_ref->{'E_new_numerator'} / $Parameters_ref->{'E_new_denominator'} ;

	# New %alpha and %beta_j, adding 1 to convergence count if convergence was obtained
	my @diff = () ;
	for my $m (keys %{$alpha_ref}) {
		if (($TMG_1_counts_ref->{$m} < 1)&&($TMG_2_counts_ref->{$m} < 1)&&($TMG_3_counts_ref->{$m} < 1)&&($TFG_1_counts_ref->{$m} < 1)) {
			# Not a single position has $m as parental genotype
			$alpha_ref->{$m} = 0 ;
		} else {
			$alpha_ref->{$m} = $alpha_new_ref->{$m} / (2*$sum_k_t_S_k_t_ref->{1} + $sum_k_t_S_k_t_ref->{2} + $sum_k_t_S_k_t_ref->{3}) ;
		}
		if ($alpha_ref->{$m}!= 0) {
			push(@diff, (abs($alpha_old{$m}-$alpha_ref->{$m})/$alpha_ref->{$m})) ;
		}
	}
	for my $n (keys %{$beta_2_ref}) {
		if ($TFG_2_counts_ref->{$n} < 1) {
			# Not a single position has $n as parental genotype for segregation type 2
			$beta_2_ref->{$n} = 0 ;
		} else {
			$beta_2_ref->{$n} = $beta_2_new_ref->{$n} / $sum_k_t_S_k_t_ref->{2} ;
		}
		if ($beta_2_ref->{$n}!= 0) {
			push(@diff, (abs($beta_2_old{$n}-$beta_2_ref->{$n})/$beta_2_ref->{$n})) ;
		}
	}
	for my $n (keys %{$beta_3_ref}) {
		if ($TFG_3_counts_ref->{$n} < 1) {
			# Not a single position has $n as parental genotype for segregation type 3
			$beta_3_ref->{$n} = 0 ;
		} else {
			$beta_3_ref->{$n} = $beta_3_new_ref->{$n} / $sum_k_t_S_k_t_ref->{3} ;
		}
		if ($beta_3_ref->{$n}!= 0) {
			push(@diff, (abs($beta_3_old{$n}-$beta_3_ref->{$n})/$beta_3_ref->{$n})) ;
		}
	}
	
	# Testing convergence
	my $convergence = "no" ;
	push(@diff, (abs($Parameters_ref->{'E_old'}-$Parameters_ref->{'E'})/$Parameters_ref->{'E'})) ;
	if (($S_counts_ref->{2} > 1)&&($Parameters_ref->{'p'}!= 0)) {
		push(@diff, (abs($Parameters_ref->{'p_old'}-$Parameters_ref->{'p'})/$Parameters_ref->{'p'})) ;
	}
	push(@diff, (abs($pi_old{1}-$pi_ref->{1})/$pi_ref->{1})) ;
	if ($S_counts_ref->{2} > 1) {
		push(@diff, (abs($pi_old{2}-$pi_ref->{2})/$pi_ref->{2})) ;
	}
	if ($S_counts_ref->{3} > 1) {
		push(@diff, (abs($pi_old{3}-$pi_ref->{3})/$pi_ref->{3})) ;
	}
	if (max(@diff)<0.01) {
		$convergence = "yes" ;
	}

	# printing new parameters
	open(PARAM, '>>'.$Parameters_ref->{'output_parameters'}) ;
	if ((defined($Parameters_ref->{'SEM'}))&&($Parameters_ref->{'run'} == 1)) {
		print(PARAM "\nS+M Iteration ".$Parameters_ref->{'run'}."\t") ;
	} elsif ((defined($Parameters_ref->{'SEM'}))&&($Parameters_ref->{'run'} > 1)&&($Parameters_ref->{'run'} < 11)) {
		print(PARAM "\nS+E+M Iteration ".$Parameters_ref->{'run'}."\t") ;
	} else {
		print(PARAM "\nE+M Iteration ".$Parameters_ref->{'run'}."\t") ;
	}
	print(PARAM $pi_ref->{1}."\t".$pi_ref->{2}."\t".$pi_ref->{3}."\t".$Parameters_ref->{'p'}."\t".$Parameters_ref->{'E'}."\t".$alpha_ref->{'AA'}."\t".$alpha_ref->{'AC'}."\t".$alpha_ref->{'AG'}."\t".$alpha_ref->{'AT'}."\t".$alpha_ref->{'CC'}."\t".$alpha_ref->{'CG'}."\t".$alpha_ref->{'CT'}."\t".$alpha_ref->{'GG'}."\t".$alpha_ref->{'GT'}."\t".$alpha_ref->{'TT'}."\t".$beta_2_ref->{'AC'}."\t".$beta_2_ref->{'CA'}."\t".$beta_2_ref->{'AG'}."\t".$beta_2_ref->{'GA'}."\t".$beta_2_ref->{'AT'}."\t".$beta_2_ref->{'TA'}."\t".$beta_2_ref->{'CG'}."\t".$beta_2_ref->{'GC'}."\t".$beta_2_ref->{'CT'}."\t".$beta_2_ref->{'TC'}."\t".$beta_2_ref->{'GT'}."\t".$beta_2_ref->{'TG'}."\t".$beta_3_ref->{'A'}."\t".$beta_3_ref->{'C'}."\t".$beta_3_ref->{'G'}."\t".$beta_3_ref->{'T'}) ;
	# printing likelihoods if parameters estimations over
	if ((defined($Parameters{'debug'})&&($Parameters_ref->{'run'}>1))||(($Parameters{'parameters_estimations_over'} eq "yes")&&(defined($Parameters{'compute_Likelihood'})))) {
		$Parameters_ref->{'L_h'} = $Parameters_ref->{'Q'} - $Parameters_ref->{'H'} ;
		print(PARAM "\t".$Parameters_ref->{'Q'}."\t".$Parameters_ref->{'H'}."\t".$Parameters_ref->{'L_h'}) ;
		if ($Parameters{'parameters_estimations_over'} eq "yes") {
			my $BIC = - 2 * $Parameters_ref->{'L_h'} + $number_free_param * log($Parameters_ref->{'sample_size'}) ;
			print(PARAM "\t".$number_free_param."\t".$Parameters_ref->{'sample_size'}."\t".$BIC) ;
		}
	}
	close(PARAM) ;

	# incrementing run number
	$Parameters_ref->{'run'}++ ;

	return $convergence ;
}


#---------------------
# Display Help message
#---------------------
sub displayHelpMessage {
	print("\n");
	print('#################################################' . "\n");
	print('#           SEX-DETector - Help section         #' . "\n"); 
	print('#################################################' . "\n");
	print('-h : print this help message' . "\n");
	print('-alr : alr input file, necessary for expression levels and X/Y sequences prediction' . "\n");
	print('-alr_gen : alr_gen input file' . "\n");
	print('-alr_gen_sum : alr_gen_sum input file' . "\n");
	print('-out : outputs base name' . "\n");
	print('-system : heterogametic system : xy or zw (only used for output headers)' . "\n");
	print('-hom : name of one of the homogametic sex individual' . "\n");
	print('-het : name of one of the heterogametic sex individual' . "\n");
	print('-hom_par : name of the homogametic parent' . "\n");
	print('-het_par : name of the heterogametic parent' . "\n");
	print('-seq : output X and Y (or W and Z or U and V) sequences for sex-linked contigs, requires alr file' . "\n");
	print('-detail : output all the details of all SNPs for all contigs' . "\n");
	print('-detail-sex-linked : output sex-linked SNPs alleles and expression levels for each individual, for sex-linked contigs' . "\n");
	print('-p : probability p of genotyping error when the Y/W is not expressed enough, default value 0.1'."\n") ;
	print('-E : probability E of genotyping error, default value 0.01'."\n") ;
	print('-pi_1 : probability of segregation type 1 (autosomal)'."\n") ;
	print('-pi_2 : probability of segregation type 2 (X/Y)'."\n") ;
	print('-pi_3 : probability of segregation type 3 (hemizygous)'."\n") ;
	print('-param : control file with input parameters'."\n") ;
	print('-L : ask for Likelihood and BIC computation once parameters were estimated'."\n") ;
	print('-skip_opt : no parameter optimization, input parameters must be provided in the control file'."\n") ;
	print('-debug : Likelihood is computed at avery iteration'."\n") ;
	print('-no_sex_chr : the model only includes autosomal segregation'."\n") ;
	print('-SEM : Stochastic Expectation Maximization algorithm for the first 10 iterations'."\n") ;
	print("\n");
	print('#################################################' . "\n");
	print("\n");
	return 0;
}


#-----------------
# Check parameters
#-----------------
sub read_control_file {

	# Recovering parameters
	my $Parameters_ref = shift;
	my $control_file_name = $Parameters_ref->{'control_file'} ;

	open (CONTROL, '<' . $control_file_name) or die ('Error: Cannot open/read file: ' . $control_file_name . "\n");
	my @control_file = <CONTROL> ;
	my @split_control_file_first_line = split(/\t/, $control_file[0]) ;
	$Parameters_ref->{'pi_1'} = $split_control_file_first_line[1] ;
	$Parameters_ref->{'pi_2'} = $split_control_file_first_line[2] ;
	$Parameters_ref->{'pi_3'} = $split_control_file_first_line[3] ;
	$Parameters_ref->{'p'} = $split_control_file_first_line[4] ;
	$Parameters_ref->{'E'} = $split_control_file_first_line[5] ;
	$Parameters_ref->{'alpha_AA'} = $split_control_file_first_line[6] ;
	$Parameters_ref->{'alpha_AC'} = $split_control_file_first_line[7] ;
	$Parameters_ref->{'alpha_AG'} = $split_control_file_first_line[8] ;
	$Parameters_ref->{'alpha_AT'} = $split_control_file_first_line[9] ;
	$Parameters_ref->{'alpha_CC'} = $split_control_file_first_line[10] ;
	$Parameters_ref->{'alpha_CG'} = $split_control_file_first_line[11] ;
	$Parameters_ref->{'alpha_CT'} = $split_control_file_first_line[12] ;
	$Parameters_ref->{'alpha_GG'} = $split_control_file_first_line[13] ;
	$Parameters_ref->{'alpha_GT'} = $split_control_file_first_line[14] ;
	$Parameters_ref->{'alpha_TT'} = $split_control_file_first_line[15] ;
	$Parameters_ref->{'beta_2_AC'} = $split_control_file_first_line[16] ;
	$Parameters_ref->{'beta_2_CA'} = $split_control_file_first_line[17] ;
	$Parameters_ref->{'beta_2_AG'} = $split_control_file_first_line[18] ;
	$Parameters_ref->{'beta_2_GA'} = $split_control_file_first_line[19] ;
	$Parameters_ref->{'beta_2_AT'} = $split_control_file_first_line[20] ;
	$Parameters_ref->{'beta_2_TA'} = $split_control_file_first_line[21] ;
	$Parameters_ref->{'beta_2_CG'} = $split_control_file_first_line[22] ;
	$Parameters_ref->{'beta_2_GC'} = $split_control_file_first_line[23] ;
	$Parameters_ref->{'beta_2_CT'} = $split_control_file_first_line[24] ;
	$Parameters_ref->{'beta_2_TC'} = $split_control_file_first_line[25] ;
	$Parameters_ref->{'beta_2_GT'} = $split_control_file_first_line[26] ;
	$Parameters_ref->{'beta_2_TG'} = $split_control_file_first_line[27] ;
	$Parameters_ref->{'beta_3_A'} = $split_control_file_first_line[28] ;
	$Parameters_ref->{'beta_3_C'} = $split_control_file_first_line[29] ;
	$Parameters_ref->{'beta_3_G'} = $split_control_file_first_line[30] ;
	$Parameters_ref->{'beta_3_T'} = $split_control_file_first_line[31] ;
	chomp($Parameters_ref->{'beta_3_T'}) ;
	close(CONTROL) ;

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
	# Display error messages for missing input files or necessary parameters
	my $error = "ok" ;
	if ((!defined($Parameters_ref->{'control_file'}))&&(defined($Parameters_ref->{'skip_optimization'}))) {
		print 'Error: No control file provided for parameters values and skipping of parameters optimization asked!' . "\n";
		$error = "bad" ;
	}
	if ((!defined($Parameters_ref->{'ALR_file'}))&&(defined($Parameters_ref->{'detail-sex-linked'}))) {
		print 'Error: No .alr input file defined when it is needed for detail-sex-linked output!' . "\n";
		$error = "bad" ;
	}
	if (!defined($Parameters_ref->{'ALR_file'})) {
		#no need for alr file apparently
	} else {
		checkInputFile($Parameters_ref->{'ALR_file'}, 'alr', $Parameters_ref);
	}
	if (!defined($Parameters_ref->{'ALR_gen_file'})) {
		print 'Error: No .alr_gen input file defined!' . "\n";
		$error = "bad" ;
	} else {
		checkInputFile($Parameters_ref->{'ALR_gen_file'}, 'alr_gen', $Parameters_ref);
	}
	if (!defined($Parameters_ref->{'ALR_gen_summary_file'})) {
		print 'Error: No .alr_gen_summary input file defined!' . "\n";
		$error = "bad" ;
	} else {
		checkInputFile($Parameters_ref->{'ALR_gen_summary_file'}, 'alr_gen_sum', $Parameters_ref);
	}
	if (!defined($Parameters_ref->{'output_file'})) {
		print 'Error: No output file defined!' . "\n";
		$error = "bad" ;
	}
 	if (!defined($Parameters_ref->{'system'})) {
		print 'Error: no heterogametic system defined (xy or zw)!' . "\n";
		$error = "bad" ;
	}
	if ((!@{$homogametic_ref})||(!@{$heterogametic_ref})) {
		print 'Error: No sex for individuals defined!' . "\n";
		$error = "bad" ;
	}
	# Checking that parameters are coherent with the dataset
	if (!defined($Parameters_ref->{'threshold'})) {
		$Parameters_ref->{'threshold'} = 0.8 ;
	}
	# probability p of genotyping error when the Y/W is not expressed enough
	if (!defined($Parameters_ref->{'p'})) {
		$Parameters_ref->{'p'} = 0.1 ;
	} elsif ($Parameters_ref->{'p'} == 0) {
		$Parameters_ref->{'non_optimized_parameters'} += 1 ;
	}
	# probability E of genotyping error for all genotypes
	if (!defined($Parameters_ref->{'E'})) {
		$Parameters_ref->{'E'} = 0.01 ;
	} elsif ($Parameters_ref->{'E'} == 0) {
		$Parameters_ref->{'non_optimized_parameters'} += 1 ;
	}
	# Segregation probabilities (default values 1/3)
	if (defined($Parameters_ref->{'no_sex_chr'})) {
		$Parameters_ref->{'compute_Likelihood'} = "yes" ;
		$Parameters_ref->{'pi'} = {1=>1,2=>0,3=>0} ;
		$Parameters_ref->{'non_optimized_parameters'} += 16 ;
		if ($Parameters_ref->{'p'} != 0) {
			$Parameters_ref->{'non_optimized_parameters'} += 1 ;
		}
	} elsif ((!defined($Parameters_ref->{'pi_1'}))&&(!defined($Parameters_ref->{'pi_2'}))&&(!defined($Parameters_ref->{'pi_3'}))) {
		$Parameters_ref->{'pi'} = {1=>1/3,2=>1/3,3=>1/3} ;
	} elsif ((defined($Parameters_ref->{'pi_1'}))&&(defined($Parameters_ref->{'pi_2'}))&&(defined($Parameters_ref->{'pi_3'}))&&(0<$Parameters_ref->{'pi_1'})&&($Parameters_ref->{'pi_1'}<=1)&&(0<=$Parameters_ref->{'pi_2'})&&($Parameters_ref->{'pi_2'}<1)&&(0<=$Parameters_ref->{'pi_3'})&&($Parameters_ref->{'pi_3'}<1)) {
		if ((($Parameters_ref->{'pi_1'} + $Parameters_ref->{'pi_2'} + $Parameters_ref->{'pi_3'}) - 1) < 0.01) {
			$Parameters_ref->{'pi'} = {1=>$Parameters_ref->{'pi_1'},2=>$Parameters_ref->{'pi_2'},3=>$Parameters_ref->{'pi_3'}} ;
		} else {
			print 'Error: sum of pi values must be one!' . "\n";
			$error = "bad" ;
		}
	} elsif ((defined($Parameters_ref->{'pi_1'}))&&(defined($Parameters_ref->{'pi_2'}))&&(0<$Parameters_ref->{'pi_1'})&&($Parameters_ref->{'pi_1'}<=1)&&(0<=$Parameters_ref->{'pi_2'})&&($Parameters_ref->{'pi_2'}<1)) {
		$Parameters_ref->{'pi'} = {1=>$Parameters_ref->{'pi_1'},2=>$Parameters_ref->{'pi_2'},3=>(1-$Parameters_ref->{'pi_1'}-$Parameters_ref->{'pi_2'})} ;
	} elsif ((defined($Parameters_ref->{'pi_1'}))&&(defined($Parameters_ref->{'pi_3'}))&&(0<$Parameters_ref->{'pi_1'})&&($Parameters_ref->{'pi_1'}<=1)&&(0<=$Parameters_ref->{'pi_3'})&&($Parameters_ref->{'pi_3'}<1)) {
		$Parameters_ref->{'pi'} = {1=>$Parameters_ref->{'pi_1'},2=>(1-$Parameters_ref->{'pi_1'}-$Parameters_ref->{'pi_3'}),3=>$Parameters_ref->{'pi_3'}} ;
	} elsif ((defined($Parameters_ref->{'pi_2'}))&&(defined($Parameters_ref->{'pi_3'}))&&(0<=$Parameters_ref->{'pi_2'})&&($Parameters_ref->{'pi_2'}<1)&&(0<=$Parameters_ref->{'pi_3'})&&($Parameters_ref->{'pi_3'}<1)) {
		$Parameters_ref->{'pi'} = {1=>(1-$Parameters_ref->{'pi_2'}-$Parameters_ref->{'pi_3'}),2=>$Parameters_ref->{'pi_2'},3=>$Parameters_ref->{'pi_3'}} ;
	} else {
		print("ERROR! Please define at least two segregation probabilities between 0 and 1 (both excluded) or none\n") ;
		$error = "bad" ;
	}
	# Check that if sequences are required, alr file was provided
	if ((defined($Parameters_ref->{'sequences'}))&&(!defined($Parameters_ref->{'ALR_file'}))) {
		print("ERROR! Please provide alr file to obtain sex-linked sequences (option -seq).\n") ;
		$error = "bad" ;
	}

	if ($error eq "bad") {
		print 'rerun with -h to see the help section' . "\n";
		exit(1);
	}
	return 0;
}

#-----------------------------
# Check that input files exist
#-----------------------------
sub checkInputFile {
	# Recovers parameters
	my ($Input_file, $Type, $Parameters_ref) = @_ ;
	
	# Check input file existence
	if ((! -e $Input_file)||(-z $Input_file)) {
		print 'Error: the ' . $Type . ' file specified in the command line does not exist or is empty! Please check file name, path and content.' . "\n";
		print 'Selected file: ' . basename($Input_file) . "\n";
		print 'directory: ' . dirname($Input_file) . "\n";
		exit(1);
	}
	
	# Check input file format
	checkInputFileFormat($Input_file, $Type, $Parameters_ref);
	return 0;
}

#-------------------------
# Check input Files format
#-------------------------
sub checkInputFileFormat {
	# Recovers parameters
	my ($Input_file, $Type, $Parameters_ref) = @_ ;
	my $heterogametic_ref = $Parameters_ref->{'heterogametic'} ;
	my $homogametic_ref = $Parameters_ref->{'homogametic'} ;
	my $alr_columns_ref = $Parameters_ref->{'alr_columns'} ;
	my $alr_gen_columns_ref = $Parameters_ref->{'alr_gen_columns'} ;

	# Initializations
	my $Counter = 0;

	# Cheks that the first lines of the alr file are well formated
	if ($Type eq 'alr') {
		open (INPUT, '<' . $Input_file) or die ('Error: Cannot open/read file: ' . $Input_file . "\n");
		while (my $Alr_line = <INPUT>){	
			if ($Counter > 2) {
				last;			
			} elsif ($Counter == 0 && $Alr_line!~ /^>([\w\.]+)/) {
				# Checks the sequence name on the first line
				print 'Error: The format of the first line of the ' . $Type . ' file specified in the command line is not valid!' . "\n";
				print 'Selected file: ' . basename($Input_file) . "\n";
				print 'Corresponding directory: ' . dirname($Input_file) . "\n";
				print 'Examples for correct first lines: >Contig23735 or >gi_291382774_ref_XM_002708152.1' . "\n";
				exit(1);
			} elsif ($Counter == 1) {
				# header line
			} elsif ($Counter == 2 && $Alr_line!~ /^[ACGT]{1}\t[MP]{1}\t[\d\[\/\]]+/) {
				# Checks the format of the third line (first position)
				print 'Error: The format of the third line of the ' . $Type . ' file specified in the command line is not a valid!' . "\n";
				print 'Selected file: ' . basename($Input_file) . "\n";
				print 'Corresponding directory: ' . dirname($Input_file) . "\n";
				print 'Example for correct line - Monomorphic position (4 individuals): T	M	3	4	3	7' . "\n";
				print 'Example for correct line - Polymorphic position (3 individuals): A	P	3[0/3/0/0]	4[3/1/0/0]	2[0/2/0/0]' . "\n";
				print 'Your invalid line: ' . $Alr_line . "\n";
				exit(1);
			}
			$Counter++;
		}
		close (INPUT);
		# Don't forget the too short files : Checks the number of line.
		if ($Counter < 3) {
			print 'Error: The ' . $Type . ' file specified in the command line contains less than 3 lines!' . "\n";
			print 'Selected file: ' . basename($Input_file) . "\n";
			print 'Corresponding directory: ' . dirname($Input_file) . "\n";
			exit(1);
		}
	}

	# Cheks that the first lines of the alr_gen file are well formated
	if ($Type eq 'alr_gen') {
		open (INPUT, '<' . $Input_file) or die ('Error: Cannot open/read file: ' . $Input_file . "\n");
		while (my $Alr_gen_line = <INPUT>){	
			if ($Counter > 2) {
				last;			
			} elsif ($Counter == 0 && $Alr_gen_line!~ /^>([\w\.]+)/) {
				# Checks the sequence name on the first line
				print 'Error: The format of the first line of the ' . $Type . ' file specified in the command line is not valid!' . "\n";
				print 'Selected file: ' . basename($Input_file) . "\n";
				print 'Corresponding directory: ' . dirname($Input_file) . "\n";
				print 'Examples for correct first lines: >Contig23735 or >gi_291382774_ref_XM_002708152.1' . "\n";
				exit(1);
				
			} elsif ($Counter == 1) {
				# header line
			} elsif ($Counter == 2 && $Alr_gen_line!~ /^.+(\t[ATGCN]{2}\|\d{1}(\.\d+)?)+/) {
				# Checks the format of the third line (first position)
				print 'Error: The format of the third line of the ' . $Type . ' file specified in the command line is not a valid!' . "\n";
				print 'Selected file: ' . basename($Input_file) . "\n";
				print 'Corresponding directory: ' . dirname($Input_file) . "\n";
				print 'Example for correct line (2 individuals) : 41	NN|0	TT|0.995604' . "\n";
				print 'Your invalid line: ' . $Alr_gen_line . "\n";
				exit(1);
			}
			$Counter++;
		}
		close (INPUT);
		# Don't forget the too short files : Checks the number of line.
		if ($Counter < 3) {
			print 'Error: The ' . $Type . ' file specified in the command line contains less than 3 lines!' . "\n";
			print 'Selected file: ' . basename($Input_file) . "\n";
			print 'Corresponding directory: ' . dirname($Input_file) . "\n";
			exit(1);
		}
	}
	# Cheks that the first lines of the alr_gen_summary file are well formated
	if ($Type eq 'alr_gen_sum') {
		open (INPUT, '<' . $Input_file) or die ('Error: Cannot open/read file: ' . $Input_file . "\n");
		while (my $Alr_gen_line = <INPUT>){	
			if ($Counter > 2) {
				last;			
			} elsif ($Counter == 0 && $Alr_gen_line!~ /^>([\w\.]+)/) {
				# Checks the sequence name on the first line
				print 'Error: The format of the first line of the ' . $Type . ' file specified in the command line is not valid!' . "\n";
				print 'Selected file: ' . basename($Input_file) . "\n";
				print 'Corresponding directory: ' . dirname($Input_file) . "\n";
				print 'Examples for correct first lines: >Contig23735 or >gi_291382774_ref_XM_002708152.1' . "\n";
				exit(1);
				
			} elsif ($Counter == 1) {
				# header line
			} elsif ($Counter == 2 && $Alr_gen_line!~ /^[\d]+(\t[ATGCN]{2})+/) {
				# Checks the format of the third line (first position)
				print 'Error: The format of the third line of the ' . $Type . ' file specified in the command line is not a valid!' . "\n";
				print 'Selected file: ' . basename($Input_file) . "\n";
				print 'Corresponding directory: ' . dirname($Input_file) . "\n";
				print 'Example for correct line (2 individuals) : 41	NN	TT' . "\n";
				print 'Your invalid line: ' . $Alr_gen_line . "\n";
				exit(1);
			}
			$Counter++;
		}
		close (INPUT);
		# Don't forget the too short files : Checks the number of line.
		if ($Counter < 3) {
			print 'Error: The ' . $Type . ' file specified in the command line contains less than 3 lines!' . "\n";
			print 'Selected file: ' . basename($Input_file) . "\n";
			print 'Corresponding directory: ' . dirname($Input_file) . "\n";
			exit(1);
		}
	}
	return 0;
}

#---------------------------
# Generate output file names
#---------------------------
sub generateOutputfileNames {
	# Recovering parameters
	my $Parameters_ref = shift;
	my $Output_file = $Parameters_ref->{'output_file'} ;

	# Creating main output file name
	$Parameters_ref->{'output_file_name'} = $Output_file . '_assignment.txt';
	# writing the header
	open(OUTPUT, '>'.$Parameters_ref->{'output_file_name'}) ;
	if ($Parameters_ref->{'system'} eq 'xy') {
		print(OUTPUT "contig\tprobability_autosomal\tprobability_sex-linked\tprobability_hemizygous\tassignment\tnumber_SNP_without_error\tnumber_autosomal_SNPs_without_error\tnumber_autosomal_SNPs_with_error\tnumber_X_Y_SNPs_without_error\tnumber_X_Y_SNPs_with_error\tnumber_hemizygous_SNPs_without_error\tnumber_hemizygous_SNPs_with_error") ;
		if (defined($Parameters_ref->{'ALR_file'})) {
			print(OUTPUT "\t".'number_clean_XY_SNPs_without_error'."\t".'number_clean_hemizygous_SNPs_without_error'."\n") ;
		} else {
			print(OUTPUT "\n") ;
		}
	} elsif ($Parameters_ref->{'system'} eq 'zw') {
		print(OUTPUT "contig\tprobability_autosomal\tprobability_sex-linked\tprobability_hemizygous\tassignment\tnumber_SNP_without_error\tnumber_autosomal_SNPs_without_error\tnumber_autosomal_SNPs_with_error\tnumber_Z_W_SNPs_without_error\tnumber_Z_W_SNPs_with_error\tnumber_hemizygous_SNPs_without_error\tnumber_hemizygous_SNPs_with_error") ;
		if (defined($Parameters_ref->{'ALR_file'})) {
			print(OUTPUT "\t".'number_clean_ZW_SNPs_without_error'."\t".'number_clean_hemizygous_SNPs_without_error'."\n") ;
		} else {
			print(OUTPUT "\n") ;
		}
	}
	close(OUTPUT) ;

	# Creating output for parameters values
	$Parameters_ref->{'output_parameters'} = $Output_file . '_parameters.txt';
	# writing the header
	open(PARAM, '>'.$Parameters_ref->{'output_parameters'}) ;
	print(PARAM "iteration number\tautosomal_probability\tXY_or_ZW_probability\themizygous_probability\tp\tE\talpha_AA\talpha_AC\talpha_AG\talpha_AT\talpha_CC\talpha_CG\talpha_CT\talpha_GG\talpha_GT\talpha_TT\tbeta_2_AC\tbeta_2_CA\tbeta_2_AG\tbeta_2_GA\tbeta_2_AT\tbeta_2_TA\tbeta_2_CG\tbeta_2_GC\tbeta_2_CT\tbeta_2_TC\tbeta_2_GT\tbeta_2_TG\tbeta_3_A\tbeta_3_C\tbeta_3_G\tbeta_3_T") ;
	if ((defined($Parameters{'debug'}))||(defined($Parameters{'compute_Likelihood'}))) {
		print(PARAM "\tQ\tH\tL\tnumber_free_parameters\tsample_size\tBIC") ;
	}
	close(PARAM) ;

	# Creating optional output file names that were asked for
	if (defined($Parameters_ref->{'sequences'})) {
		$Parameters_ref->{'sequences_file_name'} = $Output_file . '_sex-linked_sequences.fasta';
		# clearing eventual previous outputs from file
		open(SEQ, ">".$Parameters_ref->{'sequences_file_name'}) ;
		close(SEQ) ;
	}
	if (defined($Parameters_ref->{'detail'})) {
		$Parameters_ref->{'detail_file_name'} = $Output_file . '_SNPs_detail.txt';
		# writing the header
		open(DETAIL, ">".$Parameters_ref->{'detail_file_name'}) ;
		my $homogametic_ref = $Parameters{'homogametic'} ;
		my @homogametic = @{$homogametic_ref} ;
		my $heterogametic_ref = $Parameters{'heterogametic'} ;
		my @heterogametic = @{$heterogametic_ref} ;
		if ($Parameters_ref->{'system'} eq 'xy') {
			if (defined($Parameters_ref->{'ALR_file'})) {
				print(DETAIL "contig_name\tposition\tautosomal_proba\tXY_proba\themizygous_proba\tinferred_het_par_gen,hom_par_gen_autosomal\tautosomal_type\tinferred_het_par_genXY,hom_par_gen_XX\tXY_type\tinferred_het_par_gen,hom_par_gen_hemizygous\themizygous_type\t#_autosomal_errors\t#_XY_errors\t#_hemizygous_errors\t#_Y_errors\t#_individuals_with_aberrant_reads\t".$Parameters_ref->{'heterogametic_parent_name'}."_obs_gen\t".$Parameters_ref->{'heterogametic_parent_name'}."_expr\t".$Parameters_ref->{'homogametic_parent_name'}."_obs_gen\t".$Parameters_ref->{'homogametic_parent_name'}."_expr") ;
				foreach my $hom_ind (@homogametic) {
					print(DETAIL "\t".$hom_ind."_obs_gen\t".$hom_ind."_expr") ;
				}
				foreach my $het_ind (@heterogametic) {
					print(DETAIL "\t".$het_ind."_obs_gen\t".$het_ind."_expr") ;
				}
				print(DETAIL "\n") ;
			} else {
				print(DETAIL "contig_name\tposition\tautosomal_proba\tXY_proba\themizygous_proba\tinferred_het_par_gen,hom_par_gen_autosomal\tautosomal_type\tinferred_het_par_genXY,hom_par_gen_XX\tXY_type\tinferred_het_par_gen,hom_par_gen_hemizygous\themizygous_type\t#_autosomal_errors\t#_XY_errors\t#_hemizygous_errors\t#_Y_errors\t".$Parameters_ref->{'heterogametic_parent_name'}."_obs_gen\t".$Parameters_ref->{'homogametic_parent_name'}."_obs_gen") ;
				foreach my $hom_ind (@homogametic) {
					print(DETAIL "\t".$hom_ind."_obs_gen") ;
				}
				foreach my $het_ind (@heterogametic) {
					print(DETAIL "\t".$het_ind."_obs_gen") ;
				}
				print(DETAIL "\n") ;
			}
		} elsif ($Parameters_ref->{'system'} eq 'zw') {
			if (defined($Parameters_ref->{'ALR_file'})) {
				print(DETAIL "contig_name\tposition\tautosomal_proba\tZW_proba\themizygous_proba\tinferred_het_par_gen,hom_par_gen_autosomal\tautosomal_type\tinferred_het_par_genZW,hom_par_gen_ZZ\tZW_type\tinferred_het_par_gen,hom_par_gen_hemizygous\themizygous_type\t#_autosomal_errors\t#_ZW_errors\t#_hemizygous_errors\t#_W_errors\t#_individuals_with_aberrant_reads\t".$Parameters_ref->{'heterogametic_parent_name'}."_obs_gen\t".$Parameters_ref->{'heterogametic_parent_name'}."_expr\t".$Parameters_ref->{'homogametic_parent_name'}."_obs_gen\t".$Parameters_ref->{'homogametic_parent_name'}."_expr") ;
				foreach my $hom_ind (@homogametic) {
					print(DETAIL "\t".$hom_ind."_obs_gen\t".$hom_ind."_expr") ;
				}
				foreach my $het_ind (@heterogametic) {
					print(DETAIL "\t".$het_ind."_obs_gen\t".$het_ind."_expr") ;
				}
				print(DETAIL "\n") ;
			} else {
				print(DETAIL "contig_name\tposition\tautosomal_proba\tZW_proba\themizygous_proba\tinferred_het_par_gen,hom_par_gen_autosomal\tautosomal_type\tinferred_het_par_genZW,hom_par_gen_ZZ\tZW_type\tinferred_het_par_gen,hom_par_gen_hemizygous\themizygous_type\t#_autosomal_errors\t#_ZW_errors\t#_hemizygous_errors\t#_W_errors\t".$Parameters_ref->{'heterogametic_parent_name'}."_obs_gen\t".$Parameters_ref->{'homogametic_parent_name'}."_obs_gen") ;
				foreach my $hom_ind (@homogametic) {
					print(DETAIL "\t".$hom_ind."_obs_gen") ;
				}
				foreach my $het_ind (@heterogametic) {
					print(DETAIL "\t".$het_ind."_obs_gen") ;
				}
				print(DETAIL "\n") ;
			}
		}
		close(DETAIL) ;
	}
	if (defined($Parameters_ref->{'detail-sex-linked'})) {
		$Parameters_ref->{'sex-linked_detail_file_name'} = $Output_file . '_sex-linked_detail.txt';
		# writing the header
		open(SEX_DETAIL, ">".$Parameters_ref->{'sex-linked_detail_file_name'}) ;
		my $homogametic_ref = $Parameters{'homogametic'} ;
		my @homogametic = @{$homogametic_ref} ;
		my $heterogametic_ref = $Parameters{'heterogametic'} ;
		my @heterogametic = @{$heterogametic_ref} ;
		if ($Parameters_ref->{'system'} eq 'xy') {
			print(SEX_DETAIL "contig_name\tposition\tsex-linked_probability\tSNP_type\tnumber_error_1\tnumber_error_2\tnumber_error_3\tnumber_Y_error\tindividuals_with_aberrant_reads\t".$Parameters_ref->{'homogametic_parent_name'}.'_X1_allele'."\t".$Parameters_ref->{'homogametic_parent_name'}.'_X1_expr'."\t".$Parameters_ref->{'homogametic_parent_name'}.'_X2_allele'."\t".$Parameters_ref->{'homogametic_parent_name'}.'_X2_expr'."\t".$Parameters_ref->{'homogametic_parent_name'}.'_homozygous_allele'."\t".$Parameters_ref->{'homogametic_parent_name'}.'_total_expr'."\t".$Parameters_ref->{'heterogametic_parent_name'}.'_X_allele'."\t".$Parameters_ref->{'heterogametic_parent_name'}.'_X_expr'."\t".$Parameters_ref->{'heterogametic_parent_name'}.'_Y_allele'."\t".$Parameters_ref->{'heterogametic_parent_name'}.'_Y_expr'."\t".$Parameters_ref->{'heterogametic_parent_name'}.'_homozygous_allele'."\t".$Parameters_ref->{'heterogametic_parent_name'}.'_total_expr'."\t") ;
			foreach my $hom_ind (@homogametic) {
				print(SEX_DETAIL "\t".$hom_ind.'_X1_allele'."\t".$hom_ind.'_X1_expr'."\t".$hom_ind.'_X2_allele'."\t".$hom_ind.'_X2_expr'."\t".$hom_ind.'_homozygous_allele'."\t".$hom_ind.'_total_expr') ;
			}
			foreach my $het_ind (@heterogametic) {
				print(SEX_DETAIL "\t".$het_ind.'_X_allele'."\t".$het_ind.'_X_expr'."\t".$het_ind.'_Y_allele'."\t".$het_ind.'_Y_expr'."\t".$het_ind.'_homozygous_allele'."\t".$het_ind.'_total_expr') ;
			}
			print(SEX_DETAIL "\n") ;
		} elsif ($Parameters_ref->{'system'} eq 'zw') {
			print(SEX_DETAIL "contig_name\tposition\tsex-linked_probability\tSNP_type\tnumber_error_1\tnumber_error_2\tnumber_error_3\tnumber_Y_error\tindividuals_with_aberrant_reads\t".$Parameters_ref->{'homogametic_parent_name'}.'_Z1_allele'."\t".$Parameters_ref->{'homogametic_parent_name'}.'_Z1_expr'."\t".$Parameters_ref->{'homogametic_parent_name'}.'_Z2_allele'."\t".$Parameters_ref->{'homogametic_parent_name'}.'_Z2_expr'."\t".$Parameters_ref->{'homogametic_parent_name'}.'_homozygous_allele'."\t".$Parameters_ref->{'homogametic_parent_name'}.'_total_expr'."\t".$Parameters_ref->{'heterogametic_parent_name'}.'_Z_allele'."\t".$Parameters_ref->{'heterogametic_parent_name'}.'_Z_expr'."\t".$Parameters_ref->{'heterogametic_parent_name'}.'_W_allele'."\t".$Parameters_ref->{'heterogametic_parent_name'}.'_W_expr'."\t".$Parameters_ref->{'heterogametic_parent_name'}.'_homozygous_allele'."\t".$Parameters_ref->{'heterogametic_parent_name'}.'_total_expr'."\t") ;
			foreach my $hom_ind (@homogametic) {
				print(SEX_DETAIL "\t".$hom_ind.'_Z1_allele'."\t".$hom_ind.'_Z1_expr'."\t".$hom_ind.'_Z2_allele'."\t".$hom_ind.'_Z2_expr'."\t".$hom_ind.'_homozygous_allele'."\t".$hom_ind.'_total_expr') ;
			}
			foreach my $het_ind (@heterogametic) {
				print(SEX_DETAIL "\t".$het_ind.'_Z_allele'."\t".$het_ind.'_Z_expr'."\t".$het_ind.'_W_allele'."\t".$het_ind.'_W_expr'."\t".$het_ind.'_homozygous_allele'."\t".$het_ind.'_total_expr') ;
			}
			print(SEX_DETAIL "\n") ;
		}
		close(SEX_DETAIL) ;
	}
	return 0;
}


#--------------------------------------------------------
# initialize parents genotypes frequencies alpha and beta
#--------------------------------------------------------
sub initialize_alpha_and_beta {

	# Recovering parameters
	my $Parameters_ref = shift ;
	my $parents_genotypes_numbers_dataset_ref = $Parameters_ref->{'parents_genotypes_numbers_dataset'} ;
	my $alpha_ref = $Parameters_ref->{'alpha'} ;
	my $beta_2_ref = $Parameters_ref->{'beta_2'} ;
	my $beta_3_ref = $Parameters_ref->{'beta_3'} ;

	# Initializations
	my @DNA_bases = ('A', 'C', 'G', 'T') ;
	my @possible_genotypes = ('AA', 'AC', 'AG', 'AT', 'CC', 'CG', 'CT', 'GG', 'GT', 'TT') ;

	# computing total number of genotypes observed in dataset, idem for homozygous genotypes and heterozygous genotypes
	my $total_number_genotypes = 0; 
	foreach my $genotype (keys %{$parents_genotypes_numbers_dataset_ref}) { 
		$total_number_genotypes += $parents_genotypes_numbers_dataset_ref->{$genotype} ;
	}
	my $total_number_homozygous_genotypes = 0; 
	my @homozygous_genotype = ('AA', 'CC', 'GG', 'TT') ;
	foreach my $hom_gen (@homozygous_genotype) { 
		$total_number_homozygous_genotypes += $parents_genotypes_numbers_dataset_ref->{$hom_gen} ; 
	}
	my $total_number_heterozygous_genotypes = 0; 
	my @heterozygous_genotype = ('AC', 'AG', 'AT', 'CG', 'CT', 'GT') ;
	foreach my $het_gen (@heterozygous_genotype) {
		$total_number_heterozygous_genotypes += $parents_genotypes_numbers_dataset_ref->{$het_gen} ; 
	}

	# creating beta_2
	foreach my $X1 (@DNA_bases) {
		foreach my $Y (@DNA_bases) {
			# computing the probability of the real heterogametic parents genotype
			if ($X1 ne $Y) {
				if (grep(/$X1$Y/, @possible_genotypes) > 0) {
					$beta_2_ref->{$X1.$Y} = 0.5 * $parents_genotypes_numbers_dataset_ref->{$X1.$Y} / $total_number_heterozygous_genotypes ;
				} else {
					$beta_2_ref->{$X1.$Y} = 0.5 * $parents_genotypes_numbers_dataset_ref->{$Y.$X1} / $total_number_heterozygous_genotypes ;
				}
			}
		}
	}
	# creating alpha (=beta_1)
	foreach my $X1X2 (@possible_genotypes) {
		# computing the probability of the real homogametic parent genotype (or heterogametic parent in autosomal segregation)
		$alpha_ref->{$X1X2} = $parents_genotypes_numbers_dataset_ref->{$X1X2} / $total_number_genotypes ;
	}
	# creating beta_3
	foreach my $X1 (@DNA_bases) {
		foreach my $X2X3 (@possible_genotypes) {
			# computing the probability of the real parents' genotypes
			$beta_3_ref->{$X1} = $parents_genotypes_numbers_dataset_ref->{$X1.$X1} / $total_number_homozygous_genotypes ;
		}
	}

	return 0;
}


#---------------------------------------------------------
# create segregation tables for parents and tagged progeny
#---------------------------------------------------------
sub create_segregation_table_tagged_progeny {

	# Recovering parameters
	my $Parameters_ref = shift ;
	my $p = $Parameters_ref->{'p'} ;
	my $E = $Parameters_ref->{'E'} ;
	my $parents_genotypes_numbers_dataset_ref = $Parameters_ref->{'parents_genotypes_numbers_dataset'} ;

	# Initializations
	my %lambda_1;
	$Parameters_ref->{'lambda_1'} = \%lambda_1 ;
	my %lambda_1_E;
	$Parameters_ref->{'lambda_1_E'} = \%lambda_1_E ;
	my %lambda_2;
	$Parameters_ref->{'lambda_2'} = \%lambda_2 ;
	my %lambda_2_p;
	$Parameters_ref->{'lambda_2_p'} = \%lambda_2_p ;
	my %lambda_2_E;
	$Parameters_ref->{'lambda_2_E'} = \%lambda_2_E ;
	my %lambda_2_E_p;
	$Parameters_ref->{'lambda_2_E_p'} = \%lambda_2_E_p ;
	my %lambda_3;
	$Parameters_ref->{'lambda_3'} = \%lambda_3 ;
	my %lambda_3_E;
	$Parameters_ref->{'lambda_3_E'} = \%lambda_3_E ;

	my %mu_2_het;
	$Parameters_ref->{'mu_2_het'} = \%mu_2_het ;
	my %mu_2_p_het;
	$Parameters_ref->{'mu_2_p_het'} = \%mu_2_p_het ;
	my %mu_2_E_het;
	$Parameters_ref->{'mu_2_E_het'} = \%mu_2_E_het ;
	my %mu_2_E_p_het;
	$Parameters_ref->{'mu_2_E_p_het'} = \%mu_2_E_p_het ;
	my %mu_3_het;
	$Parameters_ref->{'mu_3_het'} = \%mu_3_het ;
	my %mu_3_E_het;
	$Parameters_ref->{'mu_3_E_het'} = \%mu_3_E_het ;

	my %mu_hom;
	$Parameters_ref->{'mu_hom'} = \%mu_hom ;
	my %mu_E_hom;
	$Parameters_ref->{'mu_E_hom'} = \%mu_E_hom ;

	my @DNA_bases = ('A', 'C', 'G', 'T') ;
	my @possible_genotypes = ('AA', 'AC', 'AG', 'AT', 'CC', 'CG', 'CT', 'GG', 'GT', 'TT') ;

	
	# for loops to fill the tables, lines of the tables correspond to parents' genotypes, Columns correspond to progeny genotypes
	# creating lambda_2 and mu_hom and mu_2_het
	foreach my $X1 (@DNA_bases) {
		foreach my $Y (@DNA_bases) {
			# initilizing to 0 all genotypes probabilities
			my %line_mu_2_het = ('AA'=>0, 'AC'=>0, 'AG'=>0, 'AT'=>0, 'CC'=>0, 'CG'=>0, 'CT'=>0, 'GG'=>0, 'GT'=>0, 'TT'=>0);
			my %line_mu_2_p_het = ('AA'=>0, 'AC'=>0, 'AG'=>0, 'AT'=>0, 'CC'=>0, 'CG'=>0, 'CT'=>0, 'GG'=>0, 'GT'=>0, 'TT'=>0);
			my %line_mu_2_E_het = ('AA'=>0, 'AC'=>0, 'AG'=>0, 'AT'=>0, 'CC'=>0, 'CG'=>0, 'CT'=>0, 'GG'=>0, 'GT'=>0, 'TT'=>0);
			my %line_mu_2_E_p_het = ('AA'=>0, 'AC'=>0, 'AG'=>0, 'AT'=>0, 'CC'=>0, 'CG'=>0, 'CT'=>0, 'GG'=>0, 'GT'=>0, 'TT'=>0);
			# incrementing probabilities given real parents genotype X1Y and X2X3
			foreach my $genotype (@possible_genotypes) {
				if (($genotype ne "$X1$Y")&&($genotype ne "$Y$X1")&&($genotype ne "$X1$X1")) {
					$line_mu_2_E_het{$genotype} += 1/9 ;
					$line_mu_2_E_p_het{$genotype} += 1/9 ;
				} elsif ($genotype eq "$X1$X1") {
					$line_mu_2_E_het{$genotype} += 1/9 ;
				} else {
					$line_mu_2_het{$genotype} += 1 ;
					if ($X1 ne $Y) {
						$line_mu_2_E_p_het{$genotype} += 1/9 ;
					}
				}
			}
			$line_mu_2_p_het{$X1.$X1} += 1 ;
			# Registering line of table
			$mu_2_het{$X1.$Y} = \%line_mu_2_het ;
			$mu_2_p_het{$X1.$Y} = \%line_mu_2_p_het ;
			$mu_2_E_het{$X1.$Y} = \%line_mu_2_E_het ;
			$mu_2_E_p_het{$X1.$Y} = \%line_mu_2_E_p_het ;
		}
	}
	foreach my $X2X3 (@possible_genotypes) {
		# retrieving X1 and X2
		my $X2 = substr($X2X3,0,1) ;
		my $X3 = substr($X2X3,1,1) ;
		# initilizing to 0 all genotypes probabilities
		my %line_mu_hom = ('AA'=>0, 'AC'=>0, 'AG'=>0, 'AT'=>0, 'CC'=>0, 'CG'=>0, 'CT'=>0, 'GG'=>0, 'GT'=>0, 'TT'=>0);
		my %line_mu_E_hom = ('AA'=>0, 'AC'=>0, 'AG'=>0, 'AT'=>0, 'CC'=>0, 'CG'=>0, 'CT'=>0, 'GG'=>0, 'GT'=>0, 'TT'=>0);
		# incrementing probabilities given real parents genotype X1Y and X2X3
		foreach my $genotype (@possible_genotypes) {
			if (($genotype ne "$X2$X3")&&($genotype ne "$X3$X2")) {
				$line_mu_E_hom{$genotype} += 1/9 ;
			} else {
				$line_mu_hom{$genotype} += 1 ;
			}
		}
		# Registering line of table
		$mu_hom{$X2.$X3} = \%line_mu_hom ;
		$mu_E_hom{$X2.$X3} = \%line_mu_E_hom ;
	}
	foreach my $X1 (@DNA_bases) {
		foreach my $Y (@DNA_bases) {
			foreach my $X2X3 (@possible_genotypes) {
				# retrieving X1 and X2
				my $X2 = substr($X2X3,0,1) ;
				my $X3 = substr($X2X3,1,1) ;
				# initilizing to 0 all genotypes probabilities, het for heterogametic and hom for homogametic
				my %line_lambda_2 = ('het,AA'=>0, 'het,AC'=>0, 'het,AG'=>0, 'het,AT'=>0, 'het,CC'=>0, 'het,CG'=>0, 'het,CT'=>0, 'het,GG'=>0, 'het,GT'=>0, 'het,TT'=>0, 'hom,AA'=>0, 'hom,AC'=>0, 'hom,AG'=>0, 'hom,AT'=>0, 'hom,CC'=>0, 'hom,CG'=>0, 'hom,CT'=>0, 'hom,GG'=>0, 'hom,GT'=>0, 'hom,TT'=>0);
				my %line_lambda_2_E = ('het,AA'=>0, 'het,AC'=>0, 'het,AG'=>0, 'het,AT'=>0, 'het,CC'=>0, 'het,CG'=>0, 'het,CT'=>0, 'het,GG'=>0, 'het,GT'=>0, 'het,TT'=>0, 'hom,AA'=>0, 'hom,AC'=>0, 'hom,AG'=>0, 'hom,AT'=>0, 'hom,CC'=>0, 'hom,CG'=>0, 'hom,CT'=>0, 'hom,GG'=>0, 'hom,GT'=>0, 'hom,TT'=>0);
				my %line_lambda_2_p = ('het,AA'=>0, 'het,AC'=>0, 'het,AG'=>0, 'het,AT'=>0, 'het,CC'=>0, 'het,CG'=>0, 'het,CT'=>0, 'het,GG'=>0, 'het,GT'=>0, 'het,TT'=>0);
				my %line_lambda_2_E_p = ('het,AA'=>0, 'het,AC'=>0, 'het,AG'=>0, 'het,AT'=>0, 'het,CC'=>0, 'het,CG'=>0, 'het,CT'=>0, 'het,GG'=>0, 'het,GT'=>0, 'het,TT'=>0);
				# registering X and Y alleles
				$line_lambda_2{'Y_allele'} = $Y;
				my @Xalleles = ($X1, $X2, $X3) ;
				$line_lambda_2{'X_alleles'} = \@Xalleles;
				# Registering SNP type
				if (($X1 eq $Y)&&($X1 eq $X2)&&($X1 eq $X3)) {
					$line_lambda_2{'SNP_type'} = 'monomorphic' ;
				} elsif ($X1 eq $Y) {
					#######################################################
					# ajouter ici if (($X2 ne $X1)&&($X3 ne $X1)) { #informative XX if ($X3 eq $X2), XXX if ($X2 ne $X3)
					# else not informative :
					$line_lambda_2{'SNP_type'} = 'not_informative' ;
				} elsif (($X1.$Y eq $X2.$X3)||($X1.$Y eq $X3.$X2)||($Y.$X1 eq $X2.$X3)||($Y.$X1 eq $X3.$X2)) {
					$line_lambda_2{'SNP_type'} = 'not_informative' ;
				} elsif (($X1 ne $Y)&&($X2 eq $X3)&&($X1 eq $X2)) {
					$line_lambda_2{'SNP_type'} = 'XY' ;
				} elsif (($X1 ne $Y)&&($X2 eq $X3)&&($Y eq $X2)) {
					$line_lambda_2{'SNP_type'} = 'XX' ;
				} elsif (($X1 ne $Y)&&($X2 eq $X3)&&($Y ne $X2)&&($X1 ne $X2)) {
					$line_lambda_2{'SNP_type'} = 'XXY' ;
				} elsif (($X1 ne $Y)&&($X3 ne $X2)&&(($X1 eq $X3)||($X1 eq $X2))) {
					$line_lambda_2{'SNP_type'} = 'XXY' ;
				} elsif (($X1 ne $Y)&&($X3 ne $X2)&&(($Y eq $X3)||($Y eq $X2))) {
					$line_lambda_2{'SNP_type'} = 'XXX' ;
				} elsif (($X1 ne $Y)&&($X3 ne $X2)&&($X1 ne $X3)&&($X1 ne $X2)) {
					$line_lambda_2{'SNP_type'} = 'XXXY' ;
				}
				# incrementing probabilities given real parents genotype X1Y and X2X3
				# homogametic progeny
				if ($X2 eq $X3) {
					# only one genotype in homogametic progeny X1X2
					foreach my $genotype (@possible_genotypes) {
						if (($genotype ne $X1.$X2)&&($genotype ne $X2.$X1)) {
							$line_lambda_2_E{'hom,'.$genotype} += 1/9 ;
						} else {
							$line_lambda_2{'hom,'.$genotype} += 1 ;
						}
					}
				} else {
					# two genotypes in homogametic progeny X1X2 and X1X3
					foreach my $genotype (@possible_genotypes) {
						if (($genotype ne $X1.$X2)&&($genotype ne $X2.$X1)&&($genotype ne $X1.$X3)&&($genotype ne $X3.$X1)) {
							$line_lambda_2_E{'hom,'.$genotype} += 1/9 ;
						} else {
							$line_lambda_2{'hom,'.$genotype} += 0.5 ;
							$line_lambda_2_E{'hom,'.$genotype} += 0.5*1/9 ;
						}
					}
				}
				# heterogametic progeny
				if (($X2 eq $X3)&&($Y eq $X3)) {
					# only one genotype in heterogametic progeny YX2=YX3=X2X2=X3X3
					foreach my $genotype (@possible_genotypes) {
						if ($genotype ne $Y.$X2) {
							$line_lambda_2_E{'het,'.$genotype} += 1/9 ;
							$line_lambda_2_E_p{'het,'.$genotype} += 1/9 ;
						} else {
							$line_lambda_2{'het,'.$genotype} += 1 ;
							$line_lambda_2_p{'het,'.$genotype} += 1 ;
						}
					}
				} elsif ($X2 eq $X3) {
					# two genotypes in heterogametic progeny YX2=YX3 and X2X2=X3X3
					foreach my $genotype (@possible_genotypes) {
						if (($genotype ne $Y.$X2)&&($genotype ne $X2.$Y)&&($genotype ne $X2.$X2)) {
							$line_lambda_2_E{'het,'.$genotype} += 1/9 ;
							$line_lambda_2_E_p{'het,'.$genotype} += 1/9 ;
						} elsif ($genotype eq $X2.$X2) {
							$line_lambda_2_E{'het,'.$genotype} += 1/9 ;
							$line_lambda_2_p{'het,'.$genotype} += 1 ;
						} else {
							$line_lambda_2{'het,'.$genotype} += 1 ;
							$line_lambda_2_E_p{'het,'.$genotype} += 1/9 ;
						}
					}
				} elsif (($Y eq $X2)||($Y eq $X3)) {
					# two genotypes in heterogametic progeny
					if ($Y eq $X2) {
						# three genotypes in heterogametic progeny : YX2=X2X2, YX3, X3X3
						foreach my $genotype (@possible_genotypes) {
							if (($genotype ne $X2.$X2)&&($genotype ne $X3.$Y)&&($genotype ne $Y.$X3)&&($genotype ne $X3.$X3)) {
								$line_lambda_2_E{'het,'.$genotype} += 1/9 ;
								$line_lambda_2_E_p{'het,'.$genotype} += 1/9 ;
							} elsif ($genotype eq $X3.$X3) {
								$line_lambda_2_E{'het,'.$genotype} += 1/9 ;
								$line_lambda_2_p{'het,'.$genotype} += 0.5 ;
								$line_lambda_2_E_p{'het,'.$genotype} += 0.5*1/9 ;
							} elsif ($genotype eq $X2.$X2) {
								$line_lambda_2{'het,'.$genotype} += 0.5 ;
								$line_lambda_2_p{'het,'.$genotype} += 0.5 ;
								$line_lambda_2_E{'het,'.$genotype} += 0.5*1/9 ;
								$line_lambda_2_E_p{'het,'.$genotype} += 0.5*1/9 ;
							} else { #YX3
								$line_lambda_2{'het,'.$genotype} += 0.5 ;
								$line_lambda_2_E{'het,'.$genotype} += 0.5*1/9 ;
								$line_lambda_2_E_p{'het,'.$genotype} += 1/9 ;
							}
						}
					}
					if ($Y eq $X3) {
						# three genotypes in heterogametic progeny : YX3=X3X3, YX2, X2X2
						foreach my $genotype (@possible_genotypes) {
							if (($genotype ne $X3.$X3)&&($genotype ne $X2.$Y)&&($genotype ne $Y.$X2)&&($genotype ne $X2.$X2)) {
								$line_lambda_2_E{'het,'.$genotype} += 1/9 ;
								$line_lambda_2_E_p{'het,'.$genotype} += 1/9 ;
							} elsif ($genotype eq $X2.$X2) {
								$line_lambda_2_p{'het,'.$genotype} += 0.5 ;
								$line_lambda_2_E{'het,'.$genotype} += 1/9 ;
								$line_lambda_2_E_p{'het,'.$genotype} += 0.5*1/9 ;
							} elsif ($genotype eq $X3.$X3) {
								$line_lambda_2_p{'het,'.$genotype} += 0.5 ;
								$line_lambda_2{'het,'.$genotype} += 0.5 ;
								$line_lambda_2_E{'het,'.$genotype} += 0.5*1/9 ;
								$line_lambda_2_E_p{'het,'.$genotype} += 0.5*1/9 ;
							} else {
								$line_lambda_2{'het,'.$genotype} += 0.5 ;
								$line_lambda_2_E{'het,'.$genotype} += 0.5*1/9 ;
								$line_lambda_2_E_p{'het,'.$genotype} += 1/9 ;
							}
						}
					}
				} else {
					# four genotypes in heterogametic progeny YX2, YX3, X2X2, X3X3
					foreach my $genotype (@possible_genotypes) {
						if (($genotype ne $Y.$X2)&&($genotype ne $X2.$Y)&&($genotype ne $Y.$X3)&&($genotype ne $X3.$Y)&&($genotype ne $X2.$X2)&&($genotype ne $X3.$X3)) {
							$line_lambda_2_E{'het,'.$genotype} += 1/9 ;
							$line_lambda_2_E_p{'het,'.$genotype} += 1/9 ;
						} elsif (($genotype eq $X2.$X2)||($genotype eq $X3.$X3)) {
							$line_lambda_2_E{'het,'.$genotype} += 1/9 ;
							$line_lambda_2_p{'het,'.$genotype} += 0.5 ;
							$line_lambda_2_E_p{'het,'.$genotype} += 0.5*1/9 ;
						} else {
							$line_lambda_2{'het,'.$genotype} += 0.5 ;
							$line_lambda_2_E{'het,'.$genotype} += 0.5*1/9 ;
							$line_lambda_2_E_p{'het,'.$genotype} += 1/9 ;
						}
					}
				}
				# Registering line of table
				$lambda_2{$X1.$Y.','.$X2.$X3} = \%line_lambda_2 ;
				$lambda_2_p{$X1.$Y.','.$X2.$X3} = \%line_lambda_2_p ;
				$lambda_2_E{$X1.$Y.','.$X2.$X3} = \%line_lambda_2_E ;
				$lambda_2_E_p{$X1.$Y.','.$X2.$X3} = \%line_lambda_2_E_p ;
			}
		}
	}

	# creating lambda_1
	foreach my $X1X2 (@possible_genotypes) {
		foreach my $X3X4 (@possible_genotypes) {
			# retrieving X1, X2, X3 and X4
			my $X1 = substr($X1X2,0,1) ;
			my $X2 = substr($X1X2,1,1) ;
			my $X3 = substr($X3X4,0,1) ;
			my $X4 = substr($X3X4,1,1) ;
			# initilizing to 0 all parents (par) and progeny's (pro) genotypes probabilities, het for heterogametic and hom for homogametic
			my %line_lambda_1 = ('AA'=>0, 'AC'=>0, 'AG'=>0, 'AT'=>0, 'CC'=>0, 'CG'=>0, 'CT'=>0, 'GG'=>0, 'GT'=>0, 'TT'=>0) ;
			my %line_lambda_1_E = ('AA'=>0, 'AC'=>0, 'AG'=>0, 'AT'=>0, 'CC'=>0, 'CG'=>0, 'CT'=>0, 'GG'=>0, 'GT'=>0, 'TT'=>0) ;
			# Registering SNP type
			if (($X1 eq $X2)&&($X1 eq $X4)&&($X1 eq $X3)) {
				$line_lambda_1{'SNP_type'} = "monomorphic" ;
			} elsif ($X1 eq $X2) {
				$line_lambda_1{'SNP_type'} = "not_informative" ;
			} elsif (($X1.$X2 eq $X3.$X4)||($X1.$X2 eq $X4.$X3)||($X2.$X1 eq $X3.$X4)||($X2.$X1 eq $X4.$X3)) {
				$line_lambda_1{'SNP_type'} = "not_informative" ;
			} else {
				$line_lambda_1{'SNP_type'} = "informative" ;
			}
			# computing probabilities given real parents genotype X1X2 and X3X4
			# progeny
			if (($X1 eq $X2)&&($X3 eq $X4)) {
				# only one genotype in progeny X1X3 = X1X4 = X2X3 = X2X4
				foreach my $genotype (@possible_genotypes) {
					if (($genotype ne $X1.$X3)&&($genotype ne $X3.$X1)) {
						$line_lambda_1_E{$genotype} += 1/9 ;
					} else {
						$line_lambda_1{$genotype} += 1 ;
					}
				}
			} elsif (($X1 eq $X2)||($X3 eq $X4)) {
				# 2 genotypes in progeny 
				if ($X1 eq $X2) {
					# X1X3 = X2X3 and X1X4 = X2X4
					foreach my $genotype (@possible_genotypes) {
						if (($genotype ne $X1.$X3)&&($genotype ne $X3.$X1)&&($genotype ne $X1.$X4)&&($genotype ne $X4.$X1)) {
							$line_lambda_1_E{$genotype} += 1/9 ;
						} else {
							$line_lambda_1{$genotype} += 0.5 ;
							$line_lambda_1_E{$genotype} += 0.5*1/9 ;
						}
					}
				} else {
					# X1X3 = X1X4 and X2X3 = X2X4
					foreach my $genotype (@possible_genotypes) {
						if (($genotype ne $X1.$X3)&&($genotype ne $X3.$X1)&&($genotype ne $X2.$X3)&&($genotype ne $X3.$X2)) {
							$line_lambda_1_E{$genotype} += 1/9 ;
						} else {
							$line_lambda_1{$genotype} += 0.5 ;
							$line_lambda_1_E{$genotype} += 0.5*1/9 ;
						}
					}
				}
			} elsif ((($X1 eq $X3)&&($X2 eq $X4))||(($X1 eq $X4)&&($X2 eq $X3))) {
				# 3 genotypes in progeny X1X1 (1/4) X1X2 (1/2) X2X2 (1/4)
				foreach my $genotype (@possible_genotypes) {
					if (($genotype ne $X1.$X1)&&($genotype ne $X1.$X2)&&($genotype ne $X2.$X1)&&($genotype ne $X2.$X2)) {
						$line_lambda_1_E{$genotype} += 1/9 ;
					} elsif (($genotype eq $X1.$X2)||($genotype eq $X2.$X1)) {
						$line_lambda_1{$genotype} += 0.5 ;
						$line_lambda_1_E{$genotype} += 2*0.25*1/9 ;
					} else {
						$line_lambda_1{$genotype} += 0.25 ;
						$line_lambda_1_E{$genotype} += 0.25*1/9 + 0.5*1/9 ;
					}
				}
			} else {
				# 4 genotypes in progeny X1X3, X1X4, X2X3 and X2X4
				foreach my $genotype (@possible_genotypes) {
					if (($genotype ne $X1.$X3)&&($genotype ne $X3.$X1)&&($genotype ne $X1.$X4)&&($genotype ne $X4.$X1)&&($genotype ne $X2.$X3)&&($genotype ne $X3.$X2)&&($genotype ne $X2.$X4)&&($genotype ne $X4.$X2)) {
						$line_lambda_1_E{$genotype} += 1/9 ;
					} else {
						$line_lambda_1{$genotype} += 0.25 ;
						$line_lambda_1_E{$genotype} += 0.25*3*1/9 ;
					}
				}
			}
			# Registering line of table
			$lambda_1{$X1.$X2.','.$X3.$X4} = \%line_lambda_1 ;
			$lambda_1_E{$X1.$X2.','.$X3.$X4} = \%line_lambda_1_E ;
		}
	}

	# creating lambda_3 and mu_3
	foreach my $X1 (@DNA_bases) {
		# initilizing to 0 all genotypes probabilities
		my %line_mu_3_het = ('AA'=>0, 'AC'=>0, 'AG'=>0, 'AT'=>0, 'CC'=>0, 'CG'=>0, 'CT'=>0, 'GG'=>0, 'GT'=>0, 'TT'=>0);
		my %line_mu_3_E_het = ('AA'=>0, 'AC'=>0, 'AG'=>0, 'AT'=>0, 'CC'=>0, 'CG'=>0, 'CT'=>0, 'GG'=>0, 'GT'=>0, 'TT'=>0);
		# incrementing probabilities given real parents genotype X1 and X2X3
		foreach my $genotype (@possible_genotypes) {
			if ($genotype ne "$X1$X1") {
				$line_mu_3_E_het{$genotype} += 1/9 ;
			} else {
				$line_mu_3_het{$genotype} += 1 ;
			}
		}
		# Registering line of table
		$mu_3_het{$X1} = \%line_mu_3_het ;
		$mu_3_E_het{$X1} = \%line_mu_3_E_het ;
	}
	foreach my $X1 (@DNA_bases) {
		foreach my $X2X3 (@possible_genotypes) {
			# retrieving X2 and X3
			my $X2 = substr($X2X3,0,1) ;
			my $X3 = substr($X2X3,1,1) ;
			# initilizing to 0 all genotypes probabilities, het for heterogametic and hom for homogametic
			my %line_lambda_3 = ('het,AA'=>0, 'het,AC'=>0, 'het,AG'=>0, 'het,AT'=>0, 'het,CC'=>0, 'het,CG'=>0, 'het,CT'=>0, 'het,GG'=>0, 'het,GT'=>0, 'het,TT'=>0, 'hom,AA'=>0, 'hom,AC'=>0, 'hom,AG'=>0, 'hom,AT'=>0, 'hom,CC'=>0, 'hom,CG'=>0, 'hom,CT'=>0, 'hom,GG'=>0, 'hom,GT'=>0, 'hom,TT'=>0);
			my %line_lambda_3_E = ('het,AA'=>0, 'het,AC'=>0, 'het,AG'=>0, 'het,AT'=>0, 'het,CC'=>0, 'het,CG'=>0, 'het,CT'=>0, 'het,GG'=>0, 'het,GT'=>0, 'het,TT'=>0, 'hom,AA'=>0, 'hom,AC'=>0, 'hom,AG'=>0, 'hom,AT'=>0, 'hom,CC'=>0, 'hom,CG'=>0, 'hom,CT'=>0, 'hom,GG'=>0, 'hom,GT'=>0, 'hom,TT'=>0);
			# registering X alleles
			my @Xalleles ;
			if (($X1 ne $X2)&&($X1 ne $X3)&&($X2 ne $X3)) {
				@Xalleles = ($X1,$X2,$X3) ;
			} elsif (($X1 ne $X2)&&($X1 eq $X3)) {
				@Xalleles = ($X1,$X2) ;
			} elsif (($X1 eq $X2)&&($X1 ne $X3)) {
				@Xalleles = ($X1,$X3) ;
			} elsif (($X1 ne $X2)&&($X2 eq $X3)) {
				@Xalleles = ($X1,$X2) ;
			} elsif (($X1 eq $X2)&&($X1 eq $X3)) {
				@Xalleles = ($X1) ;
			}
			@Xalleles = ($X1, $X2, $X3) ;
			$line_lambda_3{'X_alleles'} = \@Xalleles;
			# Registering SNP type
			if (($X1 eq $X2)&&($X1 eq $X3)) {
				$line_lambda_3{'SNP_type'} = "monomorphic" ;
			} elsif (($X2 eq $X3)&&($X1 ne $X2)) {
				$line_lambda_3{'SNP_type'} = "XX0" ;
			} elsif (($X1 eq $X3)||($X1 eq $X2)) {
				$line_lambda_3{'SNP_type'} = "XX0" ;
			} elsif (($X3 ne $X2)&&($X1 ne $X3)&&($X1 ne $X2)) {
				$line_lambda_3{'SNP_type'} = "XXX0" ;
			}
			# incrementing probabilities given real parents genotype X1 and X2X3
			# heterogametic progeny
			if ($X2 eq $X3) {
				# only one genotype in heterogametic progeny X2X2=X3X3
				foreach my $genotype (@possible_genotypes) {
					if ($genotype ne $X2.$X2) {
						$line_lambda_3_E{'het,'.$genotype} += 1/9 ;
					} else {
						$line_lambda_3{'het,'.$genotype} += 1 ;
					}
				}
			} else {
				# two genotypes in heterogametic progeny X2X2 and X3X3
				if (($X1 eq $X2)||($X1 eq $X3)) {
					# this is easily counfounded with an autosomal segregation pattern --> increase probability of the genotype observed only in hemizygous segregation to counter that
					if ($X1 eq $X2) {
						# genotype X3 in heterogametic sex can only be obtained with a hemizygous segregation
						foreach my $genotype (@possible_genotypes) {
							if (($genotype ne $X2.$X2)&&($genotype ne $X3.$X3)) {
								$line_lambda_3_E{'het,'.$genotype} += 1/9 ;
							} elsif ($genotype eq $X2.$X2) {
								$line_lambda_3{'het,'.$genotype} += 0.25 ;
								$line_lambda_3_E{'het,'.$genotype} += 0.75*1/9 ;
							} elsif ($genotype eq $X3.$X3) {
								$line_lambda_3{'het,'.$genotype} += 0.75 ;
								$line_lambda_3_E{'het,'.$genotype} += 0.25*1/9 ;
							}
						}
					} elsif ($X1 eq $X3) {
						# genotype X2 in heterogametic sex can only be obtained with a hemizygous segregation
						foreach my $genotype (@possible_genotypes) {
							if (($genotype ne $X2.$X2)&&($genotype ne $X3.$X3)) {
								$line_lambda_3_E{'het,'.$genotype} += 1/9 ;
							} elsif ($genotype eq $X3.$X3) {
								$line_lambda_3{'het,'.$genotype} += 0.25 ;
								$line_lambda_3_E{'het,'.$genotype} += 0.75*1/9 ;
							} elsif ($genotype eq $X2.$X2) {
								$line_lambda_3{'het,'.$genotype} += 0.75 ;
								$line_lambda_3_E{'het,'.$genotype} += 0.25*1/9 ;
							}
						}
					}
				} else {
					foreach my $genotype (@possible_genotypes) {
						if (($genotype ne $X2.$X2)&&($genotype ne $X3.$X3)) {
							$line_lambda_3_E{'het,'.$genotype} += 1/9 ;
						} else {
							$line_lambda_3{'het,'.$genotype} += 0.5 ;
							$line_lambda_3_E{'het,'.$genotype} += 0.5*1/9 ;
						}
					}
				}
			}
			# homogametic progeny
			if ($X2 eq $X3) {
				# only one genotype in homogametic progeny X1X2
				foreach my $genotype (@possible_genotypes) {
					if (($genotype ne $X1.$X2)&&($genotype ne $X2.$X1)) {
						$line_lambda_3_E{'hom,'.$genotype} += 1/9 ;
					} else {
						$line_lambda_3{'hom,'.$genotype} += 1 ;
					}
				}
			} else {
				# two genotypes in homogametic progeny X1X2 and X1X3
				foreach my $genotype (@possible_genotypes) {
					if (($genotype ne $X1.$X2)&&($genotype ne $X2.$X1)&&($genotype ne $X1.$X3)&&($genotype ne $X3.$X1)) {
						$line_lambda_3_E{'hom,'.$genotype} += 1/9 ;
					} else {
						$line_lambda_3{'hom,'.$genotype} += 0.5 ;
						$line_lambda_3_E{'hom,'.$genotype} += 0.5*1/9;
					}
				}
			}
			# Registering line of table
			$lambda_3{$X1.','.$X2.$X3} = \%line_lambda_3 ;
			$lambda_3_E{$X1.','.$X2.$X3} = \%line_lambda_3_E ;
		}
	}

	return 0;
}


