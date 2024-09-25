# Filter RNA-sequencing TE data

# Title: Filter RNA-sequencing TE data_240809

# Author: Sandra Hellberg
# Email: sandra.hellberg@liu.se

# Date: 2024-08-09
# Last modified: 2024-08-09

# Description: Filter TE RNAseq data to remove strange chromosomes

rm(list=ls()) # remove all entries in the global environment 
set.seed(524) # seed for reproducibility

#-------------------------------------------------------------------------------

              ### Set directory structure and load packages 

#-------------------------------------------------------------------------------
setwd("P:/LiU/1.Projekt/2.Nuvarande/3.CN/HMA/Analysis/RNAseq/TElocal/results/TE/rds/")
rna_dir <- "P:/LiU/1.Projekt/2.Nuvarande/3.CN/HMA/Analysis/RNAseq/TElocal/results/TE/rds/"
fig_dir <- "P:/LiU/1.Projekt/2.Nuvarande/3.CN/HMA/Analysis/Figures/pdf/"

lib_path <- .libPaths("C:/Rpackages/") # specifing the package library
#update.packages(lib.loc = lib_path, ask = FALSE, checkBuilt = TRUE, dependencies = TRUE)
#BiocManager::install(lib.loc = lib_path, update = TRUE, ask = FALSE)

pack_R <- c("stringr", "tidyr", "dplyr") # libraries to load

for (i in 1:length(pack_R)) {
  library(pack_R[i], character.only = TRUE)
}

#-------------------------------------------------------------------------------

                  ### Read in and format data 

#-------------------------------------------------------------------------------

TE_locations <- readRDS("gtf_rmsk_TElocal.rds") # locations of all TEs

TE_locations <- TE_locations[ -grep("chrUn|random|alt", TE_locations$chr), ] # remove all strange chromosomes
table(TE_locations$chr) # check so that no strange chromosomes remain


# LOUCY
loucy_input <- list.files("P:/LiU/1.Projekt/2.Nuvarande/3.CN/HMA/Analysis/RNAseq/TElocal/results/TE/rds/", pattern = "_L")

loucy_names <- loucy_input %>% # create a vector with the file names to be used in the for loop
    str_remove(".rds")

print(loucy_names)

for(i in 1:length(loucy_input)) {
  df <- readRDS(paste0("P:/LiU/1.Projekt/2.Nuvarande/3.CN/HMA/Analysis/RNAseq/TElocal/results/TE/rds/", loucy_input[i]))
  df$Name <- rownames(df)
  df <- merge(df, TE_locations, by.x = "transcript", by.y = "TE")
  x <- df[ -grep("chrY", df$chr),]
  if(nrow(x) == 0) {
    assign(loucy_names[i], df)
    saveRDS(df,paste0(rna_dir, loucy_names[i], "_clean" , ".rds"))
  } else {
    assign(loucy_names[i], x)
    saveRDS(x,paste0(rna_dir, loucy_names[i], "_clean" , ".rds"))
  }
  rm(df)
  rm(x)
}    


# SUP-T1
supt1_input <- list.files("P:/LiU/1.Projekt/2.Nuvarande/3.CN/HMA/Analysis/RNAseq/TElocal/results/TE/rds/", pattern = "_S")

supt1_names <- supt1_input %>% # create a vector with the file names to be used in the for loop
  str_remove(".rds")

print(supt1_names)

for(i in 1:length(supt1_input)) {
  df <- readRDS(paste0("P:/LiU/1.Projekt/2.Nuvarande/3.CN/HMA/Analysis/RNAseq/TElocal/results/TE/rds/", supt1_input[i]))
  df$Name <- rownames(df)
  df <- merge(df, TE_locations, by.x = "transcript", by.y = "TE")
  x <- df[ -grep("chrY", df$chr),]
  if(nrow(x) == 0) {
    assign(supt1_names[i], df)
    saveRDS(df,paste0(rna_dir, supt1_names[i], "_clean" , ".rds"))
  } else {
    assign(supt1_names[i], x)
    saveRDS(x,paste0(rna_dir, supt1_names[i], "_clean" , ".rds"))
  }
  rm(df)
  rm(x)
}    






  
 