#!/usr/bin/Rscript

# Author : Suzie Tallon

# Date: February 2026

# This script identifies SNPs likely originating from duplicated regions based on allele patterns

library(tidyverse)
library(readr)
library(tibble)

# Retrieve command line arguments: SNP detail file and gene list
args <- commandArgs(trailingOnly = TRUE)
fichier_snp_detail <- args[1]
fichier_gen <- args[2]

# Load gene list and format contig names to match SNP file
list_gen<-read.table(fichier_gen, header = FALSE, col.names = c("1"), sep = "\t", stringsAsFactors = FALSE, comment.char="") 
list_gen$X1 <- gsub("\\.fasta$", "", list_gen$X1)

# Load SNP detail table and standardize contig names
snp_to_exam <- read.table(fichier_snp_detail, header = TRUE,  sep = "\t", stringsAsFactors = FALSE, comment.char="") 
snp_to_exam$contig_name <- sub("_","",snp_to_exam$contig_name)

# Initialize containers for duplicated SNPs and per-gene proportions
snp_duplicate <- c()
prop_gen <- c()

# Process each contig independently
for (fichier in list_gen$X1) {
  fichier_snp <- fichier
  snp_to_exem_filtre<- filter(snp_to_exam, contig_name==fichier_snp) # Subset SNPs for the current contig
  snp_to_exem_filtre <- snp_to_exem_filtre %>% 
    select(contig_name,position,contains("obs")) # Keep only genotype observation columns
  nb_col<- ncol(snp_to_exem_filtre)
  nb_row <- nrow(snp_to_exem_filtre)
  
  snp_to_exem_filtre <- snp_to_exem_filtre %>%
    filter(.[[3]] == .[[4]]) %>% # Keep SNPs where both parents share the same genotype
    filter(substr(.[[3]], 1, 1) != substr(.[[3]], 2, 2))  %>% # Keep heterozygous parental genotypes
    mutate(proportion_identique = rowSums(.[, 5:nb_col] == .[[3]]) / (nb_col - 4)) %>% # Compute the proportion of offspring identical to the parental genotype
    filter(proportion_identique >= 0.7) %>% # Retain SNPs with excess identical genotypes, indicative of duplication
    select(contig_name, position, proportion_identique)
  snp_duplicate <- rbind(snp_duplicate,snp_to_exem_filtre)

  # Store proportion of suspicious SNPs per contig
  nb_gene_suspicious <-nrow(snp_to_exem_filtre)
  info_gene<- c(fichier,nb_gene_suspicious/nb_row)
  prop_gen <- rbind(prop_gen,info_gene)
}

prop_gen <- as.data.frame(prop_gen)
# Define output directory from input SNP file location
output_directory <- dirname(fichier_snp_detail)

# write output files
# Write duplicated SNP list
output_file <- file.path(output_directory, "snp_duplicate.csv")
write_csv(snp_duplicate, output_file)
# Write per-contig proportion of duplicated SNPs
output_file2 <- file.path(output_directory, "prop_snp_duplicate.csv")
write_csv(prop_gen, output_file2)

