#!/usr/bin/Rscript

# Author : Suzie Tallon

# Date: February 2026

# This script examines candidate Y-linked SNPs and identifies those with population genotype patterns inconsistent with strict Y-linkage.

library(tidyverse)
library(ggplot2)
library(readr)
library(tibble)
library(conflicted)

conflicts_prefer(dplyr::filter())

# Retrieve command line arguments: SNP detail file, XY gene list, and directory containing .gen files
args <- commandArgs(trailingOnly = TRUE)
fichier_snp_detail <- args[1]
fichier_gene_XY <- args[2]
directory_path <- args[3]

# Load SNP detail table and list of XY contigs
donnees_snp_detail <- read_delim(file = fichier_snp_detail, delim = "\t")
donnees_gene_XY <- read_csv(fichier_gene_XY,  col_names = FALSE)

# Initialize output table storing SNPs with inconsistent sex-linked genotype patterns
snp_checked <- tibble( contig_name = character(), position = double(), Xm = character(), Ym = character(), Xf1 = character(), Xf2 = character(), Yfemcount = double(), noYmalcount = double() ,noNcountfemale=double(),noNcountmale=double())

# Process each candidate XY contig
for (fichier in donnees_gene_XY$X1){
  
  # Load genotype file for the current contig
  complete_path <- file.path(directory_path, paste0(fichier, ".gen"))
  donnees_gen <- read_delim(file = complete_path, delim = "\t")
 
  # Select SNPs inferred as ZW-linked with reliable genotype inference
  snp_to_check <- donnees_snp_detail %>% filter(contig_name==fichier) %>% 
    mutate(`inferred_het_par_genXY,hom_par_gen_XX` = gsub(",", "", `inferred_het_par_genXY,hom_par_gen_XX`))%>% 
    separate(col = `inferred_het_par_genXY,hom_par_gen_XX`, into = c("Xm", "Ym", "Xf1", "Xf2"), sep =  "(?!^)") %>%
    filter(XY_proba == pmax(XY_proba, autosomal_proba, hemizygous_proba))%>% # Retain SNPs best supported by the ZW segregation model
    filter(`#_individuals_with_aberrant_reads`==0) %>%# Remove SNPs with aberrant read patterns
    filter(XY_type %in% c("XY", "XXY", "XXXY")) %>%# Restrict to canonical ZW-compatible genotype configurations
    select(contig_name, position, Xm, Ym, Xf1, Xf2)
  
  if (nrow(snp_to_check) > 0) {
    # Examine segregation pattern in the genotype file
    for (i in 1:nrow(snp_to_check)) {
      Yfemcount=0
      noYmalcount=0
      noNcountfemale=0
      noNcountmale=0
      snp_data_gen <- donnees_gen %>% filter(position==snp_to_check$position[i])
      # Iterate over individuals and count Y allele presence/absence by sex
      for (colonne in names(snp_data_gen)[-1]) {
        if (!grepl("N",snp_data_gen[1,][[colonne]]) & grepl("female",colonne) ) {
          noNcountfemale=noNcountfemale+1  # Female individuals without missing genotype
          if (grepl(snp_to_check$Ym[i],snp_data_gen[1,][[colonne]])) { 
            Yfemcount=Yfemcount+1 # Count females carrying the Y allele 
          } 
        } else if (!grepl("N",snp_data_gen[1,][[colonne]]) & grepl("_male",colonne) ) {
          noNcountmale=noNcountmale+1  # Male individuals without missing genotype
          if  (!grepl(snp_to_check$Ym[i],snp_data_gen[1,][[colonne]])) {
          noYmalcount=noYmalcount+1 # Count Males lacking the Y allele 
          }
        }
      }
      # Store results for the current SNP
      snp_checked <- bind_rows(snp_checked, tibble(
        contig_name = snp_to_check$contig_name[i],
        position = snp_to_check$position[i],
        Xm = snp_to_check$Xm[i],
        Ym = snp_to_check$Ym[i],
        Xf1 =snp_to_check$Xf1[i] , 
        Xf2 =snp_to_check$Xf2[i],
        Yfemcount=Yfemcount,
        noYmalcount=noYmalcount,
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



