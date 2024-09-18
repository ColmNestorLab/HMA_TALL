# Title: Figure 5c
# Author: Sandra Hellberg
# Email: sandra.hellberg@liu.se

# Date: 2024-02-02
# Last modified: 2024-08-14

rm(list=ls()) # remove all entries in the global environment 
set.seed(524) # seed for reproducibility

#-------------------------------------------------------------------------------

            ## Set directory structure and load packages---- 

#-------------------------------------------------------------------------------
setwd("set your working directory here")
meth_dir <- "directory for methylation files"
rna_dir <- "directory for rna-seq files"
fig_dir <- "folder for output of figures"

lib_path <- .libPaths("C:/Rpackages/") # specifing the package library
#update.packages(lib.loc = lib_path, ask = FALSE, checkBuilt = TRUE, dependencies = TRUE)
#BiocManager::install(lib.loc = lib_path, update = TRUE, ask = FALSE)

pack_R <- c("stringr", "RColorBrewer", "ggplot2", "ggfortify", "ggpubr", "grid", "reshape",
            "ggsci", "scales", "BioVenn", "gridExtra", "eulerr") # libraries to load

for (i in 1:length(pack_R)) {
  library(pack_R[i], character.only = TRUE)
}

#-------------------------------------------------------------------------------

                          ## METHYLATION----

#-------------------------------------------------------------------------------
# Read the data from csv
meth_loucy <- read.delim("Supplementary_Table_15.csv", sep = ",") # LOUCY
meth_supt1 <- read.delim("Supplementary_Table_16.csv", sep = ",") # SUP-T1

# Function to filter rows by sample pattern and class_id type, and select specific columns
filter_meth <- function(data, sample_id) {
  data %>%
    filter(grepl(sample_id, sample),
           qvalue < 0.01, 
           meth.diff <= -25) %>%
    select(transcript_id)  # Select only the `transcript_id` column
}

# Filter LOUCY dataset
TE_LD3_meth <- filter_meth(meth_loucy, "LD3")
TE_LG3_meth <- filter_meth(meth_loucy, "LG3")
TE_LG7_meth <- filter_meth(meth_loucy, "LG7")

# Filter SUPT1 dataset
TE_SD3_meth <- filter_meth(meth_supt1, "SD3")
TE_SG3_meth <- filter_meth(meth_supt1, "SG3")
TE_SG7_meth <- filter_meth(meth_supt1, "SG7")

# Define background for LOUCY and SUP-T1
background_LOUCY <- unique(meth_loucy$transcript_id)
background_SUPT1 <- unique(meth_supt1$transcript_id)

#-------------------------------------------------------------------------------
                          ### Plotting----
#-------------------------------------------------------------------------------
##### Eular diagrams of overlap----
#### LOUCY
# LD3
LD3_list <- list(
  "LD3" = TE_LD3_meth$transcript,
  "background" = background_LOUCY)

fit_LD3_meth <- euler(LD3_list)

LD3_meth_plot <- plot(fit_LD3_meth, fills = TRUE, edges = TRUE, quantities = TRUE)

# LG3
LG3_list <- list(
  "LG3" = TE_LG3_meth$transcript,
  "background" = background_LOUCY)

fit_LG3_meth <- euler(LG3_list)

LG3_meth_plot <- plot(fit_LG3_meth, fills = TRUE, edges = TRUE, quantities = TRUE)

# LG7
LG7_list <- list(
  "LG7" = TE_LG7_meth$transcript,
  "background" = background_LOUCY)

fit_LG7_meth <- euler(LG7_list)

LG7_meth_plot <- plot(fit_LG7_meth, fills = TRUE, edges = TRUE, quantities = TRUE)

#### SUPT1
# SD3
SD3_list <- list(
  "SD3" = TE_SD3_meth$transcript,
  "background" = background_SUPT1)

fit_SD3_meth <- euler(SD3_list)

SD3_meth_plot <- plot(fit_SD3_meth, fills = TRUE, edges = TRUE, quantities = TRUE)

# SG3
SG3_list <- list(
  "SG3" = TE_SG3_meth$transcript,
  "background" = background_SUPT1)

fit_SG3_meth <- euler(SG3_list)

SG3_meth_plot <- plot(fit_SG3_meth, fills = TRUE, edges = TRUE, quantities = TRUE)

# SG7
SG7_list <- list(
  "SG7" = TE_SG7_meth$transcript,
  "background" = background_SUPT1)

fit_SG7_meth <- euler(SG7_list)

SG7_meth_plot <- plot(fit_SG7_meth, fills = TRUE, edges = TRUE, quantities = TRUE)

#-------------------------------------------------------------------------------

                                  ## RNA-SEQ -----

#-------------------------------------------------------------------------------
# Read the data from csv
rna_loucy <- read.delim(paste0(rna_dir,"Supplementary_Table_17.csv"), sep = "," ) # LOUCY
rna_supt1 <- read.delim(paste0(rna_dir,"Supplementary_Table_18.csv"), sep = "," ) # SUPT1


# Function to filter rows by sample pattern and class_id type, and select specific columns
filter_rna <- function(data, sample_id) {
  data %>%
    filter(grepl(sample_id, sample), 
           logFC > 1, 
           padjust < 0.05) %>%
    select(transcript)  # Select only the `transcript_id` column
}
### Combine TEs
# LOUCY
TE_LD3_RNA <- filter_rna(rna_loucy, "LD3")
TE_LG3_RNA <- filter_rna(rna_loucy, "LG3")
TE_LG7_RNA <- filter_rna(rna_loucy, "LG7")

# SUPT1
TE_SD3_RNA <- filter_rna(rna_supt1, "SD3")
TE_SG3_RNA <- filter_rna(rna_supt1, "SG3")
TE_SG7_RNA <- filter_rna(rna_supt1, "SG7")

# Define background for LOUCY and SUP-T1
background_LOUCY_RNA <- unique(rna_loucy$transcript)
background_SUPT1_RNA <- unique(rna_supt1$transcript)
#-------------------------------------------------------------------------------
                            ### Plotting----
#-------------------------------------------------------------------------------
##### Eular diagrams of overlap----
# LOUCY
# LD3
LD3_list_RNA <- list(
  "LD3" = TE_LD3_RNA$transcript,
  "background" = background_LOUCY_RNA)

fit_LD3_RNA <- euler(LD3_list_RNA)
LD3_RNA_plot <- plot(fit_LD3_RNA, fills = TRUE, edges = TRUE, quantities = TRUE)

# LG3
LG3_list_RNA <- list(
  "LG3" = TE_LG3_RNA$transcript,
  "background" = background_LOUCY_RNA)

fit_LG3_RNA <- euler(LG3_list_RNA)
LG3_RNA_plot <- plot(fit_LG3_RNA, fills = TRUE, edges = TRUE, quantities = TRUE)

# LG7
LG7_list_RNA <- list(
  "LG7" = TE_LG7_RNA$transcript,
  "background" = background_LOUCY_RNA)

fit_LG7_RNA <- euler(LG7_list_RNA)
LG7_RNA_plot <- plot(fit_LG7_RNA, fills = TRUE, edges = TRUE, quantities = TRUE)

# SUPT1
# SD3
SD3_list_RNA <- list(
  "SD3" = TE_SD3_RNA$transcript,
  "background" = background_SUPT1_RNA)

fit_SD3_RNA <- euler(SD3_list_RNA)
SD3_RNA_plot <- plot(fit_SD3_RNA, fills = TRUE, edges = TRUE, quantities = TRUE)

# SG3
SG3_list_RNA <- list(
  "SG3" = TE_SG3_RNA$transcript,
  "background" = background_SUPT1_RNA)

fit_SG3_RNA <- euler(SG3_list_RNA)
SG3_RNA_plot <- plot(fit_SG3_RNA, fills = TRUE, edges = TRUE, quantities = TRUE)

# SG7
SG7_list_RNA <- list(
  "SG7" = TE_SG7_RNA$transcript,
  "background" = background_SUPT1_RNA)

fit_SG7_RNA <- euler(SG7_list_RNA)
SG7_RNA_plot <- plot(fit_SG7_RNA, fills = TRUE, edges = TRUE, quantities = TRUE)

#-------------------------------------------------------------------------------

                      ## Combination methylation and RNA-seq-----

#-------------------------------------------------------------------------------
# LOUCY
# LD3
LD3_list_combined <- list(
  "LD3_meth" = TE_LD3_meth$transcript,
  "LD3_RNA" = TE_LD3_RNA$transcript)

fit_LD3_combined <- euler(LD3_list_combined)
LD3_combined_plot <- plot(fit_LD3_combined, fills = TRUE, edges = TRUE, quantities = TRUE)

# LG3
LG3_list_combined <- list(
  "LG3_meth" = TE_LG3_meth$transcript,
  "LG3_RNA" = TE_LG3_RNA$transcript)

fit_LG3_combined <- euler(LG3_list_combined)
LG3_combined_plot <- plot(fit_LG3_combined, fills = TRUE, edges = TRUE, quantities = TRUE)

# LG7
LG7_list_combined <- list(
  "LG7_meth" = TE_LG7_meth$transcript,
  "LG7_RNA" = TE_LG7_RNA$transcript)

fit_LG7_combined <- euler(LG7_list_combined)
LG7_combined_plot <- plot(fit_LG7_combined, fills = TRUE, edges = TRUE, quantities = TRUE)

# SUPT1
# SD3
SD3_list_combined <- list(
  "SD3_meth" = TE_SD3_meth$transcript,
  "SD3_RNA" = TE_SD3_RNA$transcript)

fit_SD3_combined <- euler(SD3_list_combined)
SD3_combined_plot <- plot(fit_SD3_combined, fills = TRUE, edges = TRUE, quantities = TRUE)

# SG3
SG3_list_combined <- list(
  "SG3_meth" = TE_SG3_meth$transcript,
  "SG3_RNA" = TE_SG3_RNA$transcript)

fit_SG3_combined <- euler(SG3_list_combined)
SG3_combined_plot <- plot(fit_SG3_combined, fills = TRUE, edges = TRUE, quantities = TRUE)

# SG7
SG7_list_combined <- list(
  "SG7_meth" = TE_SG7_meth$transcript,
  "SG7_RNA" = TE_SG7_RNA$transcript)

fit_SG7_combined <- euler(SG7_list_combined)
SG7_combined_plot <- plot(fit_SG7_combined, fills = TRUE, edges = TRUE, quantities = TRUE)


#-------------------------------------------------------------------------------

                                ### Export plots -----

#-------------------------------------------------------------------------------

# DEC 3 days
plot <- ggarrange(LD3_meth_plot, LD3_RNA_plot, LD3_combined_plot,
                  SD3_meth_plot, SD3_RNA_plot, SD3_combined_plot, ncol = 3, nrow = 2)
pdf(paste0(fig_dir, "Eular_meth_RNA_DEC3.pdf"), width = 15, height = 10)
plot(plot)
dev.off()

# GSK300 3 days
plot <- ggarrange(LG3_meth_plot, LG3_RNA_plot, LG3_combined_plot,
                  SG3_meth_plot, SG3_RNA_plot, SG3_combined_plot, ncol = 3, nrow = 2)

pdf(paste0(fig_dir, "Eular_meth_RNA_GSK3.pdf"), width = 15, height = 10)
plot(plot)
dev.off()

# GSK300 7 days
plot <- ggarrange(LG7_meth_plot, LG7_RNA_plot, LG7_combined_plot,
                  SG7_meth_plot, SG7_RNA_plot, SG7_combined_plot, ncol = 3, nrow = 2)

pdf(paste0(fig_dir, "Eular_meth_RNA_GSK7.pdf"), width = 15, height = 10)
plot(plot)
dev.off()

# SCRIPT ENDS HERE##############################################################

