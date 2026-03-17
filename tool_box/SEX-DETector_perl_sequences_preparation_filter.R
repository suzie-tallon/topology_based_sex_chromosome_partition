#!/usr/bin/env Rscript

# developper : Aline Muyle & Suzie Tallon
# last update : November 2025
# This code reads SEX-DETector perl output files and create a file that can be used with SDpop wxyz_genotyper code to extract X/Y or Z/W sequences.
# command line :
# SEX-DETector_perl_sequences_preparation.R working_directory_path assignment_file SNP_detail_file output_file
# example: SEX-DETector_perl_sequences_preparation.R /home/muyle/Documents/enseignement/Physalia/sex_chromosomes/SEX-DETector_STAR_conica Silene_latifolia_Family_perl_assignment.txt Silene_latifolia_Family_perl_SNPs_detail.txt perl_sequence_extraction/Silene_latifolia_Family_perl_allele_freq.txt

# output file format :
# #>contig_name prob_auto prob_xhemi prob_xy
# #position alleles fx_mean fy_mean
# >AmTr_v6.0_c1.110.1 1.040570e-04
# 166 C,T 0.00000 1.000000
# the above line means C has a frequency 0 on the X and 1 on the Y (fixed), conversely the T allele is fixed on the X.
# 180 A,G 0.50000 1.000000
# the above line means A has a frequency 0.5 on the X (polymorphic A/G) and 1 on the Y (fixed), allele G allele also has frequency 0.5 on the X


args <- commandArgs(trailingOnly = TRUE)
DirectoryPath <- args[1] # DirectoryPath 
assignment_filename <- args[2] # assignment_filename <- "Silene_latifolia_Family_perl_assignment.txt"
snp_detail_filename <- args[3] # snp_detail_filename <- "Silene_latifolia_Family_perl_SNPs_detail.txt"
snp_to_remove <- args[4] #snp_to_remove <-"snp_to_remove.csv"
output_filename <- args[5] # output_filename <- "perl_sequence_extraction/Silene_latifolia_Family_perl_allele_freq.txt"


library(data.table)
setwd(DirectoryPath)
Assignments <- fread(assignment_filename, h=T)
SNPs <- fread(snp_detail_filename, h=T, sep="\t")
SNP_to_remove <- fread(snp_to_remove, h=T, sep=",")

# print header of the output file:
header <- "#>contig_name prob_auto prob_sex_linked\n#position alleles fx_mean fy_mean"
write(header, file=output_filename, append=FALSE)

# keep only sex-linked genes, discard others
Assignments <- Assignments[assignment == 'sex-linked'] 

# remove SNPs where individuals have aberrant read patterns
SNPs <- SNPs[`#_individuals_with_aberrant_reads` == 0]

# remove SNPs where population individuals do not fit with sex-linkage (i.e. females with a Y allele in their genotype)

SNPs <- merge(SNPs, SNP_to_remove, by = c("contig_name", "position"), all = TRUE)
SNPs <- SNPs[Yfemcount == 0 | is.na(Yfemcount)]
#SNPs <- SNPs[noYmalcount == 0 | is.na(noYmalcount)] (males with no Y allele in their genotype)
SNPs <- SNPs[, !c("Xm", "Ym", "Xf1","Xf2","Yfemcount","noYmalcount","noNcountfemale","noNcountmale"), with = FALSE]


# loop on all sex-linked genes
for (i in seq(1, nrow(Assignments))) {
  # gene info
  CurrentGene <- Assignments$contig[i]
  Current_probability_autosomal <- Assignments$probability_autosomal[i]
  Current_probability_sex_linked <- Assignments$`probability_sex-linked`[i]
  Current_probability_hemizygous <- Assignments$probability_hemizygous[i]
  Current_number_XY_SNPs_total <- Assignments$number_X_Y_SNPs_without_error[i] +  Assignments$number_X_Y_SNPs_with_error[i]
  
  # if gene has zero X/Y SNPs (with or without error) then discard it, it is a X-hemizygous gene, no Y sequence can be output
  if (Current_number_XY_SNPs_total > 0) {
    
    # extract gene SNPs
    GeneSNPs <- SNPs[contig_name == CurrentGene]
    
    if (nrow(GeneSNPs) > 0) {
      # find the maximum probability of segregation for each SNP
      GeneSNPs[, max:=pmax(autosomal_proba, XY_proba, hemizygous_proba)]
      GeneSNPs[, SNP_segregation := if (max == autosomal_proba) {'autosomal'} else if (max == XY_proba) {'XY'} else if (max == hemizygous_proba) {'X_hemi'}, by=1:nrow(GeneSNPs)]
      
      # remove X-hemizygous SNPs because they have no Y sequence to output, they will have an X in the sequence output
      GeneSNPs <- GeneSNPs[SNP_segregation != 'X_hemi']
      
      # count number of XY SNPs
      numberXY_SNPs <- nrow(GeneSNPs[(SNP_segregation == "XY")&(XY_type != 'not_informative')])
      
      # test if the contig has any XY SNPs, if not the gene is X-hemizygous, so no Y sequence to output, stop there
      if (numberXY_SNPs > 0) {
        # print contig info line in output
        line <- paste('>', CurrentGene, ' ', Current_probability_autosomal, ' ', Current_probability_hemizygous+Current_probability_sex_linked, sep="")
        write(line, file=output_filename, append=TRUE)
        
        # identify the alleles at each SNP position
        GeneSNPs[, parent_genotypes := if (SNP_segregation == 'autosomal') {`inferred_het_par_gen,hom_par_gen_autosomal`} else if (SNP_segregation == 'XY') {`inferred_het_par_genXY,hom_par_gen_XX`} else if (SNP_segregation == 'X_hemi') {`inferred_het_par_gen,hom_par_gen_hemizygous`}, by=1:nrow(GeneSNPs)]
        GeneSNPs[, alleles := paste(unique(strsplit(paste(strsplit(parent_genotypes, ',')[[1]], sep='', collapse=''), '')[[1]]), sep='', collapse = ','), by=1:nrow(GeneSNPs)]
        GeneSNPs[, alleleNumber := length(strsplit(alleles, ',')[[1]]), by=1:nrow(GeneSNPs)]
        GeneSNPs[, allele_1_father := strsplit(strsplit(parent_genotypes, ',')[[1]][1], '')[[1]][1], by=1:nrow(GeneSNPs)]
        GeneSNPs[, allele_2_father := strsplit(strsplit(parent_genotypes, ',')[[1]][1], '')[[1]][2], by=1:nrow(GeneSNPs)]
        GeneSNPs[, allele_1_mother := strsplit(strsplit(parent_genotypes, ',')[[1]][2], '')[[1]][1], by=1:nrow(GeneSNPs)]
        GeneSNPs[, allele_2_mother := strsplit(strsplit(parent_genotypes, ',')[[1]][2], '')[[1]][2], by=1:nrow(GeneSNPs)]
        
        # loop on each SNP
        for (s in seq(1, nrow(GeneSNPs))) {
          CurrentSNPposition <- GeneSNPs$position[s]
          Current_allele_Number <- GeneSNPs$alleleNumber[s]
          Current_allele_1_father <- GeneSNPs$allele_1_father[s]
          Current_allele_2_father <- GeneSNPs$allele_2_father[s]
          Current_allele_1_mother <- GeneSNPs$allele_1_mother[s]
          Current_allele_2_mother <- GeneSNPs$allele_2_mother[s]
          Current_alleles <- GeneSNPs$alleles[s]
          CurrentSNPSegregation <- GeneSNPs$SNP_segregation[s]
          Current_XY_type <- GeneSNPs$XY_type[s]
          
          if (Current_allele_Number == 1) {
            # this is a monomorphic position, skip SNP
          } else if (Current_allele_Number == 2) {
            if (Current_allele_1_father == Current_allele_2_father) {
              # the father is homomorphic, prepare fx_mean and fy_mean values
              fx_mean = 0
              fy_mean = 0
              if (Current_allele_1_father == strsplit(Current_alleles, ',')[[1]][1]) {
                fx_mean = 1
                fy_mean = 1
              }
              # write SNP output line
              line <- paste(CurrentSNPposition, Current_alleles, fx_mean, fy_mean, sep="\t")
              write(line, file=output_filename, append=TRUE)
            } else {
              # the father is heteromorphic
              if (CurrentSNPSegregation == 'autosomal') {
                # we don't know which allele is X and which is Y, so we put 0.5 frequency for both
                line <- paste(CurrentSNPposition, Current_alleles, 0.5, 0.5, sep="\t")
                write(line, file=output_filename, append=TRUE)
              } else if (CurrentSNPSegregation == 'XY') {
                if (Current_XY_type == 'XY') {
                  # this is an XY SNP with aparently true X-Y divergence (at least from what we can see from the mother and father)
                  # allele_1 is the X and allele_2 is the Y, prepare output
                  fx_mean = 0
                  fy_mean = 0
                  if (Current_allele_1_father == strsplit(Current_alleles, ',')[[1]][1]) {
                    fx_mean = 1
                  }
                  if (Current_allele_2_father == strsplit(Current_alleles, ',')[[1]][1]) {
                    fy_mean = 1
                  }
                  line <- paste(CurrentSNPposition, Current_alleles, fx_mean, fy_mean, sep="\t")
                  write(line, file=output_filename, append=TRUE)
                } else if (Current_XY_type == 'XX') {
                  # this is an XX SNP, where the father heterozygosity is not due to X-Y divergence but to X polymorphism between mother and father
                  # only output the X allele which is the same as the Y, if possible
                  if ((Current_allele_2_father == Current_allele_1_mother)|(Current_allele_2_father == Current_allele_2_mother)) {
                    # the Y allele is also an X allele in the mother
                    # output Current_allele_2_father as Y and X allele
                    fx_mean = 0
                    fy_mean = 0
                    if (Current_allele_2_father == strsplit(Current_alleles, ',')[[1]][1]) {
                      fx_mean = 1
                      fy_mean = 1
                    }
                    line <- paste(CurrentSNPposition, Current_alleles, fx_mean, fy_mean, sep="\t")
                    write(line, file=output_filename, append=TRUE)
                  }
                }
              }
            }
          } else if (Current_allele_Number > 2 ) {
            # there are more than 2 alleles at this SNP
            if (CurrentSNPSegregation == 'autosomal') {
              # we don't know which allele is X and which is Y, so we put 0.5 frequency for both
              # this position will appear as an X in the sequence output after running wxyz_genotyper
              line <- paste(CurrentSNPposition, paste(Current_allele_1_father, Current_allele_2_father, sep=','), 0.5, 0.5, sep="\t")
              write(line, file=output_filename, append=TRUE)
            } else if (CurrentSNPSegregation == 'XY') {
              # allele_1 is the X and allele_2 is the Y, prepare output
              # we check that the Y allele is not present in the mother, meaning it can also be an X allele
              if ((Current_allele_2_father == Current_allele_1_mother)|(Current_allele_2_father == Current_allele_2_mother)) {
                # this SNP is due to X polymorphism, because the father Y allele is an X allele in the mother
                # we output Current_allele_2_father as both X and Y allele, i.e we ignore this SNP for X and Y sequences
                Current_alleles <- paste(Current_allele_2_father, Current_allele_1_father, sep=",")
                line <- paste(CurrentSNPposition, Current_alleles, 1, 1, sep="\t")
                write(line, file=output_filename, append=TRUE)
              } else {
                # the SNP is due to true X-Y divergence
                line <- paste(CurrentSNPposition, paste(Current_allele_1_father, Current_allele_2_father, sep=','), 1, 0, sep="\t")
                write(line, file=output_filename, append=TRUE)
              }
            }
          }
        }
      }
    }
  }
}

warnings()
