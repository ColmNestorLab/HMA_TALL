### Generating RDS files of TEs based on RNA-seq

# Title: Generate_RDS_TE_RNAseq_240806

# Author: Sandra Hellberg
# Email: sandra.hellberg@liu.se

# Date: 2024-01-24
# Last modified: 2024-08-06


rm(list=ls()) # remove all entries in the global environment 
set.seed(121) # seed for reproducibility

#-------------------------------------------------------------------------------

               ### Set directory structure and load packages 

#-------------------------------------------------------------------------------

setwd("working directory")
rna_dir <- "input directory for rna-seq" # rna seq folder

lib_path <- .libPaths("C:/Rpackages/") # specifing the package library
#update.packages(lib.loc = lib_path, ask = FALSE, checkBuilt = TRUE, dependencies = TRUE)
#BiocManager::install(lib.loc = lib_path, update = TRUE, ask = FALSE)

pack_R <- c("dplyr", "tidyverse", "stats")

for (i in 1:length(pack_R)) {
  library(pack_R[i], character.only = TRUE)
}

#-------------------------------------------------------------------------------

                      ### Read in and subset data-----

#-------------------------------------------------------------------------------
                                  ##### LOUCY----
#-------------------------------------------------------------------------------

results_L <- list.files(rna_dir, pattern = "LOUCY") 
results_L
results_L <- results_L[1:3] # select only results from the differential analysis

treatment <- c("LD3_RNA", "LG3_RNA", "LG7_RNA")

##### ERVs----
for(i in 1:length(results_L)) {
  x <- read.delim(paste0(rna_dir, results_L[i]))
  x <- x[, !(colnames(x) %in% c("logCPM", "PValue"))]
  
  ERV <- x[grep("LTR", rownames(x)),, drop = FALSE ]
  ERV$TE <- 0
  ERV$Family <- 0
  ERV$gene <- 0
  ERV$transcript <- 0
  
  for(n in 1:nrow(ERV)) {
    y <- rownames(ERV)[n]
    z <- str_split(y, ":")
    
    w <- z[[1]][4] # add class of TE
    ERV[n,3] <- w
    
    f <- z[[1]][3] # add family
    ERV[n,4] <- f
    
    g <- z[[1]][2] # add gene
    ERV[n,5] <- g
    
    t <- z[[1]][1] # add transcript
    ERV[n,6] <- t
  }
  assign(paste0("ERVs_", treatment[i]), ERV)
  saveRDS(ERV,paste0(rna_dir,"rds/", "ERVs_", treatment[i], ".rds"))
  rm(x,y,g,n,z)
}  
  
##### LINEs----
for(i in 1:length(results_L)) {
  x <- read.delim(paste0(rna_dir, results_L[i]))
  x <- x[, !(colnames(x) %in% c("logCPM", "PValue"))]
  
  LINE <- x[grep("LINE", rownames(x)),, drop = FALSE ]
  LINE$TE <- 0
  LINE$Family <- 0
  LINE$gene <- 0
  LINE$transcript <- 0
  
  for(n in 1:nrow(LINE)) {
    y <- rownames(LINE)[n]
    z <- str_split(y, ":")
    
    w <- z[[1]][4] # add class of TE
    LINE[n,3] <- w
    
    f <- z[[1]][3] # add family
    LINE[n,4] <- f
    
    g <- z[[1]][2] # add gene
    LINE[n,5] <- g
    
    t <- z[[1]][1] # add transcript
    LINE[n,6] <- t
  }
  assign(paste0("LINEs_", treatment[i]), LINE)
  saveRDS(LINE,paste0(rna_dir,"rds/","LINEs_", treatment[i], ".rds"))
  rm(x,y,g,n,z)
}
  
##### SINEs----
for(i in 1:length(results_L)) {
  x <- read.delim(paste0(rna_dir, results_L[i]))
  x <- x[, !(colnames(x) %in% c("logCPM", "PValue"))]
  
  SINE <- x[grep("SINE", rownames(x)),, drop = FALSE ]
  SINE$TE <- 0
  SINE$Family <- 0
  SINE$gene <- 0
  SINE$transcript <- 0
  
  for(n in 1:nrow(SINE)) {
    y <- rownames(SINE)[n]
    z <- str_split(y, ":")
    
    w <- z[[1]][4] # add class of TE
    SINE[n,3] <- w
    
    f <- z[[1]][3] # add family
    SINE[n,4] <- f
    
    g <- z[[1]][2] # add gene
    SINE[n,5] <- g
    
    t <- z[[1]][1] # add transcript
    SINE[n,6] <- t
  }
  assign(paste0("SINEs_", treatment[i]), SINE)
  saveRDS(SINE,paste0(rna_dir,"rds/", "SINEs_", treatment[i], ".rds"))
  rm(x,y,g,n,z)
}

rm(LINE, SINE, ERV) # remove data created in the for loop

#-------------------------------------------------------------------------------
                              ##### SUPT1----
#-------------------------------------------------------------------------------
results_S <- list.files(rna_dir, pattern = "SUPT1")
results_S <- results_S[1:3] # select only results from the differential analysis

treatment <- c("SD3_RNA", "SG3_RNA", "SG7_RNA")

#### ERVs----
for(i in 1:length(results_S)) {
  x <- read.delim(paste0(rna_dir, results_S[i]))
  x <- x[, !(colnames(x) %in% c("logCPM", "PValue"))]
  
  ERV <- x[grep("LTR", rownames(x)),, drop = FALSE ]
  ERV$TE <- 0
  ERV$Family <- 0
  ERV$gene <- 0
  ERV$transcript <- 0
  
  for(n in 1:nrow(ERV)) {
    y <- rownames(ERV)[n]
    z <- str_split(y, ":")
    
    w <- z[[1]][4] # add class of TE
    ERV[n,3] <- w
    
    f <- z[[1]][3] # add family
    ERV[n,4] <- f
    
    g <- z[[1]][2] # add gene
    ERV[n,5] <- g
    
    t <- z[[1]][1] # add transcript
    ERV[n,6] <- t
  }
  assign(paste0("ERVs_", treatment[i]), ERV)
  saveRDS(ERV,paste0(rna_dir,"rds/", "ERVs_", treatment[i], ".rds"))
  rm(x,y,g,n,z)
}  

#### LINEs----
for(i in 1:length(results_S)) {
  x <- read.delim(paste0(rna_dir, results_S[i]))
  x <- x[, !(colnames(x) %in% c("logCPM", "PValue"))]
  
  LINE <- x[grep("LINE", rownames(x)),, drop = FALSE ]
  LINE$TE <- 0
  LINE$Family <- 0
  LINE$gene <- 0
  LINE$transcript <- 0
  
  for(n in 1:nrow(LINE)) {
    y <- rownames(LINE)[n]
    z <- str_split(y, ":")
    
    w <- z[[1]][4] # add class of TE
    LINE[n,3] <- w
    
    f <- z[[1]][3] # add family
    LINE[n,4] <- f
    
    g <- z[[1]][2] # add gene
    LINE[n,5] <- g
    
    t <- z[[1]][1] # add transcript
    LINE[n,6] <- t
  }
  assign(paste0("LINEs_", treatment[i]), LINE)
  saveRDS(LINE,paste0(rna_dir,"rds/" ,"LINEs_", treatment[i], ".rds"))
  rm(x,y,g,n,z)
  }

#### SINEs----
for(i in 1:length(results_S)) {
  x <- read.delim(paste0(rna_dir, results_S[i]))
  x <- x[, !(colnames(x) %in% c("logCPM", "PValue"))]
  
  SINE <- x[grep("SINE", rownames(x)),, drop = FALSE ]
  SINE$TE <- 0
  SINE$Family <- 0
  SINE$gene <- 0
  SINE$transcript <- 0
  
  for(n in 1:nrow(SINE)) {
    y <- rownames(SINE)[n]
    z <- str_split(y, ":")
    
    w <- z[[1]][4] # add class of TE
    SINE[n,3] <- w
    
    f <- z[[1]][3] # add family
    SINE[n,4] <- f
    
    g <- z[[1]][2] # add gene
    SINE[n,5] <- g
    
    t <- z[[1]][1] # add transcript
    SINE[n,6] <- t
  }
  assign(paste0("SINEs_", treatment[i]), SINE)
  saveRDS(SINE,paste0(rna_dir,"rds/", "SINEs_", treatment[i], ".rds"))
  rm(x,y,g,n,z)
  }

rm(LINE, SINE, ERV) # remove data created in the for loop

