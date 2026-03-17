#!/usr/bin/Rscript

# Author : Suzie Tallon

# Date: February 2026

# This script is used to compute and SNP-level phylogenetic topologies across two species and an outgroup and to summarize topologies per gene.

library(tidyverse)
library(readr)
library(tibble)

# Retrieve command-line arguments defining input directories, gene list, and output naming
args <- commandArgs(trailingOnly = TRUE)
dir_species_1 <- args[1]
dir_species_2 <- args[2]
name <- args[3]
fichier_gen <- args[4]
outputdirectory <- args[5]
outgroupname <- args[6]

# Load gene list and standardize identifiers to match sequence file names
list_gen<-read.table(fichier_gen, header = FALSE, col.names = c("1"), sep = "\t", stringsAsFactors = FALSE, comment.char="") 
list_gen$X1 <- gsub("\\.fasta$", "", list_gen$X1)
#list_gen$X1 <- sub("_", "", list_gen$X1) #to uncomment if working on Artemia (gtf file)


###### 1 : gene level topology ###### 

# Extract Z/X and W/Y haplotypes from fasta-derived tables and convert them into aligned allele matrices
split_function<-function(data,fichier_data,name){
  data <- read.table(fichier_data, header = FALSE,  sep = "\t", stringsAsFactors = FALSE, comment.char="")
  data <- subset(data, rownames(data) %in% c("2", "4"))
  data <- as.data.frame(do.call(rbind, strsplit(as.character(data$V1), "")))
  colnames(data) <- 1:ncol(data)
  rownames(data) <- c(paste0("Z_", name),paste0("W_", name))
  data <-as.data.frame(t(data))
  return(data)
}

# Retain only informative sites with exactly two alleles across outgroup and both species, excluding missing data
pair_filter_function <- function(data1,data2,outgroup) {
  data <- cbind(outgroup,data1,data2) %>%
    filter_all(all_vars(!grepl("N", .))) %>%
    rowwise() %>%
    mutate(n_uniques = length(unique(c_across(everything())))) %>%
    filter(n_uniques == 2) %>%
    select(-n_uniques)
  return(data)
}

# Convert alleles into binary states relative to the outgroup to standardize topology encoding
binerize_function <- function(data) {
  result <- data
  for (i in 2:ncol(data)) {
    result[, i] <- ifelse(data[, i] == data[, 1], 0, 1)
  }
  result[, 1] <- 0  # Replace outgroup state with reference state 0
  return(result)
}

# Enumerate all possible rooted binary topologies for four ingroup haplotypes relative to the outgroup
alltrees <- expand.grid(replicate(4, c(0, 1), simplify = FALSE))
alltrees <- cbind(0, alltrees)

# Count occurrences of each possible topology across all sites within a gene
treecount <- function(data) {
  sapply(1:nrow(alltrees), function(i) {
    sum(apply(data, 1, function(x) all(x == alltrees[i,])))
  })
}

# Compute genome-wide topology counts for each gene by comparing both species to the outgroup
topo_function <- function(fichier_data1,fichier_data2) {
  topologies<-c()
  for (fichier in list_gen$X1){
    fichier_outgroup <- paste0(outgroupname, fichier, ".csv")
    fichier_species1 <- paste0(fichier_data1, fichier, ".fasta")
    fichier_species2 <- paste0(fichier_data2, fichier, ".fasta")
    if (!file.exists(fichier_outgroup)) {
      next  # Skip genes without corresponding outgroup data
    }
    outgroup <- read_csv(fichier_outgroup,  col_names = FALSE, show_col_types = FALSE)
    sp1_gen<-split_function(sp1_gen,fichier_species1,"sp1")
    sp2_gen<-split_function(sp2_gen,fichier_species2,"sp2")
    count<-treecount(binerize_function(pair_filter_function(sp1_gen,sp2_gen,outgroup)))
    topologies <- rbind(topologies,count)
  }
  return(topologies)
}

# Generate summarized topology count table with topology codes as column names
topo_sp1_sp2 <- as.data.frame(topo_function(dir_species_1,dir_species_2))
names(topo_sp1_sp2) <- apply(alltrees, 1, function(x) paste(x, collapse = "")) # Assign topology identifiers

# Save topology count summary for downstream analyses
write_csv(as.data.frame(topo_sp1_sp2),paste0(outputdirectory, name , ".csv"))



# #### 2 : position specific topology#### 
# # optional, not used anymore
# 
# # Similar to pair_filter_function but retains genomic position to allow positional topology analysis
# pair_filter_function_full <- function(data1, data2, outgroup) {
#   position <- seq_len(nrow(data1))
#   data <- cbind(position,outgroup, data1, data2) %>%
#     filter_all(all_vars(!grepl("N", .))) %>%
#     rowwise() %>%
#     mutate(n_uniques = length(unique(c_across(2:6)))) %>%
#     filter(n_uniques == 2) %>%
#     select(-n_uniques)
#   return(data)
# }
# 
# # Binary encoding preserving genomic coordinates for breakpoint detection analyses
# binerize_function_full <- function(data) {
#   result <- data
#   for (i in 3:ncol(data)) {
#     result[, i] <- ifelse(data[, i] == data[, 2], 0, 1)
#   }
#   result[, 2] <- 0  # Remplacer la première colonne par 0
#   return(result)
# }
# 
# # Assign topology identity at each individual SNP position instead of aggregating across the gene
# treecount_full <- function(data) {
#   result <- lapply(1:nrow(alltrees), function(i) {
#     apply(select(data,-position),1, function(x) all(x == alltrees[i,]))
#   })
#   result <- do.call(cbind, result)
#   result <- ifelse(result, 1, 0)  # Convert logical matches to binary indicators
#   result <- cbind(data$position,result)
#   return(result)
# }
# 
# # Generate per-position topology matrix 
# topo_function_full <- function(fichier_data1,fichier_data2) {  
#   topologies<-c()
#   
#   for (fichier in list_gen$X1){
#     fichier_outgroup <- paste0(outgroupname, fichier, ".csv")
#     fichier_species1 <- paste0(fichier_data1, fichier, ".fasta")
#     fichier_species2 <- paste0(fichier_data2, fichier, ".fasta")
#     if (!file.exists(fichier_outgroup)) {
#       next  
#     }
#     outgroup <- read_csv(fichier_outgroup,  col_names = FALSE, show_col_types = FALSE)
#     sp1_gen<-split_function(sp1_gen,fichier_species1,"sp1")
#     sp2_gen<-split_function(sp2_gen,fichier_species2,"sp2")
#     pair <- pair_filter_function_full(sp1_gen,sp2_gen,outgroup)
#     if(nrow(pair)!=0) {
#       count<-as.data.frame(treecount_full(binerize_function_full(pair))) %>% mutate(ID = fichier) %>% relocate(ID, .before = 1)
#       topologies <- rbind(topologies,count)
#     }
#   }
#   return(topologies)
# }
# 
# # Assemble final per-position topology table with explicit gene identifiers
# topo_sp1_sp2_full <- as.data.frame(topo_function_full(dir_species_1,dir_species_2))
# names(topo_sp1_sp2_full) <- c("ID","position",apply(alltrees, 1, function(x) paste(x, collapse = ""))) 
# 
# 
# # Save detailed topology matrix for breakpoint and recombination analyses
# write_csv(as.data.frame(topo_sp1_sp2_full),paste0(outputdirectory, name , "_full.csv"))
# 
# 
# 
