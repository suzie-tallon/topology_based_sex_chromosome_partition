#!/usr/bin/Rscript

# Author : Suzie Tallon

# Date: February 2026

# This script examines candidate W-linked SNPs and identifies those with population genotype patterns inconsistent with strict W-linkage.

library(tidyverse)
library(ggplot2)
library(readr)
library(tibble)
library(conflicted)

conflicts_prefer(dplyr::filter())

# Retrieve command line arguments: SNP detail file, ZW gene list, and directory containing .gen files
args <- commandArgs(trailingOnly = TRUE)
fichier_snp_detail <- args[1]
fichier_gene_ZW <- args[2]
directory_path <- args[3]

# Load SNP detail table and list of ZW contigs
donnees_snp_detail <- read_delim(file = fichier_snp_detail, delim = "\t")
donnees_gene_ZW <- read_csv(fichier_gene_ZW,  col_names = FALSE)

# Initialize output table storing SNPs with inconsistent sex-linked genotype patterns
snp_checked <- tibble( contig_name = character(), position = double(), Zf = character(), Wf = character(), Zm1 = character(), Zm2 = character(), Wmalcount = double(), noWfemcount = double() ,noNcountfemale=double(),noNcountmale=double())

# Process each candidate ZW contig
for (fichier in donnees_gene_ZW$X1){
  
  # Load genotype file for the current contig
  complete_path <- file.path(directory_path, paste0(fichier, ".gen"))
  donnees_gen <- read_delim(file = complete_path, delim = "\t")
  
  # Select SNPs inferred as ZW-linked with reliable genotype inference
  snp_to_check <- donnees_snp_detail %>% filter(contig_name==fichier) %>% 
    mutate(`inferred_het_par_genZW,hom_par_gen_ZZ` = gsub(",", "", `inferred_het_par_genZW,hom_par_gen_ZZ`))%>% 
    separate(col = `inferred_het_par_genZW,hom_par_gen_ZZ`, into = c("Zf", "Wf", "Zm1", "Zm2"), sep =  "(?!^)") %>%
    filter(ZW_proba == pmax(ZW_proba, autosomal_proba, hemizygous_proba))%>% # Retain SNPs best supported by the ZW segregation model
    filter(`#_individuals_with_aberrant_reads`==0) %>% # Remove SNPs with aberrant read patterns
    filter(ZW_type %in% c("ZW", "ZZW", "ZZZW")) %>% # Restrict to canonical ZW-compatible genotype configurations
    select(contig_name, position, Zf, Wf, Zm1, Zm2)
  
  if (nrow(snp_to_check) > 0) {
    # Examine segregation pattern in the genotype file
    for (i in 1:nrow(snp_to_check)) {
      Wmalcount=0
      noWfemcount=0
      noNcountfemale=0
      noNcountmale=0
      snp_data_gen <- donnees_gen %>% filter(position==snp_to_check$position[i])
      
      # Iterate over individuals and count W allele presence/absence by sex
      for (colonne in names(snp_data_gen)[-1]) {
        if (!grepl("N",snp_data_gen[1,][[colonne]]) & grepl("_male",colonne) ) { 
          noNcountmale=noNcountmale+1 # Male individuals without missing genotype
          if (grepl(snp_to_check$Wf[i],snp_data_gen[1,][[colonne]])) { 
            Wmalcount=Wmalcount+1 # Count males carrying the W allele 
          } 
        } else if (!grepl("N",snp_data_gen[1,][[colonne]]) & grepl("female",colonne) ) { 
          noNcountfemale=noNcountfemale+1 # Female individuals without missing genotype
          if  (!grepl(snp_to_check$Wf[i],snp_data_gen[1,][[colonne]])) { 
          noWfemcount=noWfemcount+1 # Count females lacking the W allele 
          }
        }
      }
      # Store results for the current SNP
      snp_checked <- bind_rows(snp_checked, tibble(
        contig_name = snp_to_check$contig_name[i],
        position = snp_to_check$position[i],
        Zf = snp_to_check$Zf[i],
        Wf = snp_to_check$Wf[i],
        Zm1 =snp_to_check$Zm1[i] , 
        Zm2 =snp_to_check$Zm2[i],
        Wmalcount=Wmalcount,
        noWfemcount=noWfemcount,
        noNcountfemale=noNcountfemale,
        noNcountmale=noNcountmale
      ))
    }
  }
}
# Extract output directory from SNP detail file path
output_directory <- dirname(fichier_snp_detail)

## Write list of SNPs showing inconsistent sex-linked segregation
output_file <- file.path(output_directory, "snp_to_remove.csv")
write_csv(snp_checked, output_file)



