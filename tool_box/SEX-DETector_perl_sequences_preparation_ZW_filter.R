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
DirectoryPath <- args[1] # DirectoryPath <- "/media/tallon@newcefe.newage.fr/DATAS/Alignements/A.sinica/merged$"
assignment_filename <- args[2] # assignment_filename <- "Silene_latifolia_Family_perl_assignment.txt"
snp_detail_filename <- args[3] # snp_detail_filename <- "Silene_latifolia_Family_perl_SNPs_detail.txt"
snp_to_remove<-args[4] ##snp_to_remove <-"snp_to_remove.csv"
output_filename <- args[5] # output_filename <- "perl_sequence_extraction/Silene_latifolia_Family_perl_allele_freq.txt"

library(data.table)
setwd(DirectoryPath)
Assignments <- fread(assignment_filename, h=T)
SNPs <- fread(snp_detail_filename, h=T, sep="\t")
SNP_to_remove <- fread(snp_to_remove, h=T, sep=",")

# print header of the output file:
header <- "#>contig_name prob_auto prob_sex_linked\n#position alleles fz_mean fw_mean"
write(header, file=output_filename, append=FALSE)

# keep only sex-linked genes, discard others
Assignments <- Assignments[assignment== 'sex-linked'] 

# remove SNPs where individuals have aberrant read patterns
SNPs <- SNPs[`#_individuals_with_aberrant_reads` == 0]

# remove SNPs where population individuals do not fit with sex-linkage (i.e. males with a W allele in their genotype)

SNPs <- merge(SNPs, SNP_to_remove, by = c("contig_name", "position"), all = TRUE)
SNPs <- SNPs[Wmalcount == 0 | is.na(Wmalcount)]
#SNPs <- SNPs[noWfemcount == 0 | is.na(noWfemcount)] (females with no W allele in their genotype)
SNPs <- SNPs[, !c("Zf", "Wf", "Zm1","Zm2","Wmalcount","noWfemcount","noNcountfemale","noNcountmale"), with = FALSE]

# loop on all sex-linked genes
for (i in seq(1, nrow(Assignments))) {
  # gene info
  CurrentGene <- Assignments$contig[i]
  Current_probability_autosomal <- Assignments$probability_autosomal[i]
  Current_probability_sex_linked <- Assignments$`probability_sex-linked`[i]
  Current_probability_hemizygous <- Assignments$probability_hemizygous[i]
  Current_number_ZW_SNPs_total <- Assignments$number_Z_W_SNPs_without_error[i] +  Assignments$number_Z_W_SNPs_with_error[i]
  
  # if gene has zero Z/W SNPs (with or without error) then discard it, it is a Z-hemizygous gene, no W sequence can be output
  if (Current_number_ZW_SNPs_total > 0) {
    
    # extract gene SNPs
    GeneSNPs <- SNPs[contig_name == CurrentGene]
    
    if (nrow(GeneSNPs) > 0) {
      # find the maximum probability of segregation for each SNP
      GeneSNPs[, max:=pmax(autosomal_proba, ZW_proba, hemizygous_proba)]
      GeneSNPs[, SNP_segregation := if (max == autosomal_proba) {'autosomal'} else if (max == ZW_proba) {'ZW'} else if (max == hemizygous_proba) {'Z_hemi'}, by=1:nrow(GeneSNPs)]
      
      # remove Z-hemizygous SNPs because they have no W sequence to output, they will have an Z in the sequence output
      GeneSNPs <- GeneSNPs[SNP_segregation != 'Z_hemi']
      
      # count number of ZW SNPs
      numberZW_SNPs <- nrow(GeneSNPs[(SNP_segregation == "ZW")&(ZW_type != 'not_informative')])
      
      # test if the contig has any ZW SNPs, if not the gene is Z-hemizygous, so no W sequence to output, stop there
      if (numberZW_SNPs > 0) {
        # print contig info line in output
        line <- paste('>', CurrentGene, ' ', Current_probability_autosomal, ' ', Current_probability_hemizygous+Current_probability_sex_linked, sep="")
        write(line, file=output_filename, append=TRUE)
        
        # identify the alleles at each SNP position
        GeneSNPs[, parent_genotypes := if (SNP_segregation == 'autosomal') {`inferred_het_par_gen,hom_par_gen_autosomal`} else if (SNP_segregation == 'ZW') {`inferred_het_par_genZW,hom_par_gen_ZZ`} else if (SNP_segregation == 'Z_hemi') {`inferred_het_par_gen,hom_par_gen_hemizygous`}, by=1:nrow(GeneSNPs)]
        GeneSNPs[, alleles := paste(unique(strsplit(paste(strsplit(parent_genotypes, ',')[[1]], sep='', collapse=''), '')[[1]]), sep='', collapse = ','), by=1:nrow(GeneSNPs)]
        GeneSNPs[, alleleNumber := length(strsplit(alleles, ',')[[1]]), by=1:nrow(GeneSNPs)]
        GeneSNPs[, allele_1_mother := strsplit(strsplit(parent_genotypes, ',')[[1]][1], '')[[1]][1], by=1:nrow(GeneSNPs)]
        GeneSNPs[, allele_2_mother := strsplit(strsplit(parent_genotypes, ',')[[1]][1], '')[[1]][2], by=1:nrow(GeneSNPs)]
        GeneSNPs[, allele_1_father := strsplit(strsplit(parent_genotypes, ',')[[1]][2], '')[[1]][1], by=1:nrow(GeneSNPs)]
        GeneSNPs[, allele_2_father := strsplit(strsplit(parent_genotypes, ',')[[1]][2], '')[[1]][2], by=1:nrow(GeneSNPs)]
        
        # loop on each SNP
        for (s in seq(1, nrow(GeneSNPs))) {
          CurrentSNPposition <- GeneSNPs$position[s]
          Current_allele_Number <- GeneSNPs$alleleNumber[s]
          Current_allele_1_mother <- GeneSNPs$allele_1_mother[s]
          Current_allele_2_mother <- GeneSNPs$allele_2_mother[s]
          Current_allele_1_father <- GeneSNPs$allele_1_father[s]
          Current_allele_2_father <- GeneSNPs$allele_2_father[s]
          Current_alleles <- GeneSNPs$alleles[s]
          CurrentSNPSegregation <- GeneSNPs$SNP_segregation[s]
          Current_ZW_type <- GeneSNPs$ZW_type[s]
          
          if (Current_allele_Number == 1) {
            # this is a monomorphic position, skip SNP
          } else if (Current_allele_Number == 2) {
            if (Current_allele_1_mother == Current_allele_2_mother) {
              # the mother is homomorphic, prepare fx_mean and fy_mean values
              fz_mean = 0
              fw_mean = 0
              if (Current_allele_1_mother == strsplit(Current_alleles, ',')[[1]][1]) {
                fz_mean = 1
                fw_mean = 1
              }
              # write SNP output line
              line <- paste(CurrentSNPposition, Current_alleles, fz_mean, fw_mean, sep="\t")
              write(line, file=output_filename, append=TRUE)
            } else {
              # the mother is heteromorphic
              if (CurrentSNPSegregation == 'autosomal') {
                # we don't know which allele is Z and which is W, so we put 0.5 frequency for both
                line <- paste(CurrentSNPposition, Current_alleles, 0.5, 0.5, sep="\t")
                write(line, file=output_filename, append=TRUE)
              } else if (CurrentSNPSegregation == 'ZW') {
                if (Current_ZW_type == 'ZW') {
                  # this is an ZW SNP with aparently true Z-W divergence (at least from what we can see from the father and mother)
                  # allele_1 is the Z and allele_2 is the W, prepare output
                  fz_mean = 0
                  fw_mean = 0
                  if (Current_allele_1_mother == strsplit(Current_alleles, ',')[[1]][1]) {
                    fz_mean = 1
                  }
                  if (Current_allele_2_mother == strsplit(Current_alleles, ',')[[1]][1]) {
                    fw_mean = 1
                  }
                  line <- paste(CurrentSNPposition, Current_alleles, fz_mean, fw_mean, sep="\t")
                  write(line, file=output_filename, append=TRUE)
                } else if (Current_ZW_type == 'ZZ') {
                  # this is an ZZ SNP, where the mother heterozygosity is not due to Z-W divergence but to Z polymorphism between father and mother
                  # only output the Z allele which is the same as the W, if possible
                  if ((Current_allele_2_mother == Current_allele_1_father)|(Current_allele_2_mother == Current_allele_2_father)) {
                    # the W allele is also an Z allele in the father
                    # output Current_allele_2_mother as W and Z allele
                    fz_mean = 0
                    fw_mean = 0
                    if (Current_allele_2_mother == strsplit(Current_alleles, ',')[[1]][1]) {
                      fz_mean = 1
                      fw_mean = 1
                    }
                    line <- paste(CurrentSNPposition, Current_alleles, fz_mean, fw_mean, sep="\t")
                    write(line, file=output_filename, append=TRUE)
                  }
                }
              }
            }
          } else if (Current_allele_Number > 2 ) {
            # there are more than 2 alleles at this SNP
            if (CurrentSNPSegregation == 'autosomal') {
              # we don't know which allele is Z and which is W, so we put 0.5 frequency for both
              # this position will appear as an Z in the sequence output after running wxyz_genotyper
              line <- paste(CurrentSNPposition, paste(Current_allele_1_mother, Current_allele_2_mother, sep=','), 0.5, 0.5, sep="\t")
              write(line, file=output_filename, append=TRUE)
            } else if (CurrentSNPSegregation == 'ZW') {
              # allele_1 is the Z and allele_2 is the W, prepare output
              # we check that the Z allele is not present in the father, meaning it can also be an W allele
              if ((Current_allele_2_mother == Current_allele_1_father)|(Current_allele_2_mother == Current_allele_2_father)) {
                # this SNP is due to Z polymorphism, because the mother W allele is an Z allele in the father
                # we output Current_allele_2_mother as both Z and W allele, i.e we ignore this SNP for Z and W sequences
                Current_alleles <- paste(Current_allele_2_mother, Current_allele_1_mother, sep=",")
                line <- paste(CurrentSNPposition, Current_alleles, 1, 1, sep="\t")
                write(line, file=output_filename, append=TRUE)
              } else {
                # the SNP is due to true Z-W divergence
                line <- paste(CurrentSNPposition, paste(Current_allele_1_mother, Current_allele_2_mother, sep=','), 1, 0, sep="\t")
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
