#!/usr/bin/Rscript

# Author : Suzie Tallon

# Date: February 2026

#script to remove UTRs of coding sequences in XY species

library(tidyverse)
library(Biostrings)

# Retrieve the command line argument for the folder
args <- commandArgs(trailingOnly = TRUE)
dossier <- args[1]  # Le premier argument correspondra au dossier

# Import the gff3 file
fichier_gff <- "/data/work/suzie/Silene/genome_reference/S.latifolia_v4.0.gff.gff3_polished_no_Y"

# Get the list of all files in the folder with the .fasta extension
fichiers_a_traiter <- list.files(dossier, pattern = "\\.fasta$", full.names = TRUE)

donnees_gff <- read.table(fichier_gff, sep="\t", header=FALSE, comment.char="#", stringsAsFactors=FALSE)

# Add colum names
colnames(donnees_gff) <- c("Chromosome", "Source", "Type", "Start", "End", "Score", "Strand", "Phase", "Attributes")

# Split the Attributes column into separate columns
donnees_gff <- donnees_gff %>%
  separate(col = Attributes, 
           into = c("ID", "Coverage", "Sequence_ID", "Valid_ORFs", "Extra_copy_number", "Copy_num_ID"),
           sep = ";",
           convert = TRUE, 
           extra = "drop") 

# Remove the "ID=_" prefix from the "ID" column in donnees_gff
donnees_gff$ID <- gsub("^ID=", "", donnees_gff$ID)

#Build the table that allows us to know where to cut

## If there are multiple 5' and 3' UTRs, we take the min or max of their positions
# in order to select the most "internal" UTRs within the sequence

# for the 5'
gene_data_5 <- donnees_gff %>% filter(Type=="five_prime_UTR") %>%
  select(Chromosome,Type,Strand,Start,End,ID) %>%
  mutate(ID = gsub("\\.(f|t|\t|C).*", "", ID)) %>%
  mutate(Length_UTR_5=End-Start+1) %>%
  mutate(nbr_UTR_5=1) %>%
  group_by(ID) %>%
  summarize(Max_UTR_Start_5 = max(Start),
            Max_UTR_End_5 = max(End),
            Min_UTR_Start_5 = min(Start),
            Min_UTR_End_5 = min(End),
            Length_UTR_tot_5 = sum(Length_UTR_5),
            Nb_UTR_5=sum(nbr_UTR_5))

# for the 3'
gene_data_3 <- donnees_gff %>% filter(Type=="three_prime_UTR") %>%
  select(Chromosome,Type,Strand,Start,End,ID) %>% 
  mutate(ID = gsub("\\.(f|t|\t|C).*", "", ID)) %>%
  mutate(Length_UTR_3=End-Start+1) %>%
  mutate(nbr_UTR_3=1) %>%
  group_by(ID) %>%
  summarize(Max_UTR_Start_3 = max(Start),
            Max_UTR_End_3 = max(End),
            Min_UTR_Start_3 = min(Start),
            Min_UTR_End_3=min(End),
            Length_UTR_tot_3 = sum(Length_UTR_3),
            Nb_UTR_3=sum(nbr_UTR_3))

# merge the two tables
gene_data<-gene_data_5 %>% inner_join(gene_data_3)

# Get the min and max of the coding regions. It must be done differently
# for genes on the plus strand and on the minus strand
donnees_plus <- donnees_gff %>% filter(Type=="CDS") %>%
  select(Chromosome,Type,Strand,Start,End,ID) %>%
  mutate(ID = gsub("\\.(f|t|\t|C).*", "", ID)) %>%
  filter(Strand=="+") %>%
  mutate(Length_CDS=End-Start+1) %>%
  group_by(ID) %>%
  summarize(Min = min(Start),
            Max = max(End),
            Length_CDS_tot = sum(Length_CDS),
            Strand="+") %>% # this is where we determine at which position the coding region starts and ends
  inner_join(gene_data, by ="ID") %>% # join with the table containing UTR min and max values
  mutate(Start_cut=if_else(Max_UTR_Start_5!=Min_UTR_End_5,Max_UTR_End_5-Max_UTR_Start_5+2,2)) %>% 
  mutate(End_cut=if_else(Max_UTR_Start_3!=Min_UTR_End_3,Min_UTR_End_3-Min_UTR_Start_3+1,1)) %>% 
  mutate(Length_tot=Length_CDS_tot + Length_UTR_tot_5 + Length_UTR_tot_3)


donnees_moins <- donnees_gff %>% filter(Type=="CDS") %>%
  select(Chromosome,Type,Strand,Start,End,ID) %>%
  mutate(ID = gsub("\\.(f|t|\t|C).*", "", ID)) %>%
  filter(Strand=="-") %>%
  mutate(Length_CDS=End-Start+1) %>%
  group_by(ID) %>%
  summarize(Min = min(Start),
            Max = max(End),
            Length_CDS_tot = sum(Length_CDS),
            Strand = "-" ) %>%
  inner_join(gene_data, by ="ID") %>%
  mutate(Start_cut=if_else(Max_UTR_Start_5!=Min_UTR_End_5,Min_UTR_End_5-Min_UTR_Start_5+2,2)) %>%
  mutate(End_cut=if_else(Max_UTR_Start_3!=Min_UTR_End_3,Max_UTR_End_3-Max_UTR_Start_3+1,1)) %>%
  mutate(Length_tot=Length_CDS_tot + Length_UTR_tot_5 + Length_UTR_tot_3)

donnees_cut<-bind_rows(donnees_plus,donnees_moins) %>% arrange(ID)

traiter_fichier <- function(fichier) {
  tryCatch({
    donnees_cut <- donnees_cut %>% filter(ID == paste0("_", tools::file_path_sans_ext(basename(fichier))))
    sequences <- readDNAStringSet(fichier, format = "fasta")
    
    sequences_modifiees <- lapply(seq_along(sequences), function(i) {
      sequence <- sequences[[i]]
      if (donnees_cut$Nb_UTR_5>1|donnees_cut$Nb_UTR_3>1 & nchar(sequence)==donnees_cut$Length_tot){
        sequence <- subseq(sequence, start = donnees_cut$Length_UTR_tot_5+1, end = nchar(sequence) - donnees_cut$Length_UTR_tot_3)
      } else { 
        sequence <- subseq(sequence, start = donnees_cut$Start_cut, end = nchar(sequence) - donnees_cut$End_cut)
      }
      return(sequence)
    })
    
    # Convert the list of DNAString to XStringSet
    sequences_modifiees <- DNAStringSet(sequences_modifiees)
    
    names(sequences_modifiees) <- names(sequences)
    
    nouveau_dossier <- file.path(dirname(fichier), "cutted_fasta")
    
    # Create the new folder if it does not already exist
    if (!file.exists(nouveau_dossier)) {
      dir.create(nouveau_dossier)
    }
    
    # Save modified sequences in a new file
    output_file <- file.path(nouveau_dossier, paste0(tools::file_path_sans_ext(basename(fichier)), ".fasta"))
    
    # Open the file in write mode
    con <- file(output_file, "w")
    
    # Loop through modified sequences and write them to the file
    for (i in seq_along(sequences_modifiees)) {
      cat(">", names(sequences_modifiees)[i], "\n", sep = "", file = con)
      writeLines(as.character(sequences_modifiees[i]), con)
    }
    
    close(con)
    
    
    # Create the log file in the same folder as the output files
    fichier_log <- file.path(nouveau_dossier, "log.txt")
    if (!file.exists(fichier_log)) {
      file.create(fichier_log)
    }
  }, error = function(e) {
    message("Erreur avec le fichier : ", fichier)
    message("Message d'erreur : ", conditionMessage(e))
    message("Enregistrement de l'erreur dans le fichier de journal...")
    fichier_log <- file.path(dirname(fichier), "cutted_fasta", "log.txt")
    write(paste("Erreur avec le fichier:", fichier, "\n", "Message d'erreur:", conditionMessage(e), "\n\n"), file = fichier_log, append = TRUE)
  })
}

for (fichier in fichiers_a_traiter) {
  traiter_fichier(fichier)
}
