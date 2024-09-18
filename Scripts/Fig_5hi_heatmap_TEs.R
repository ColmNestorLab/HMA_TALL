# Title: Figure 5H,I
# Author: Sandra Hellberg
# Email: sandra.hellberg@liu.se

# Date: 2024-08-07
# Last modified: 2024-08-15

# Description: Heatmaps of the four largest subfamilies of ERVs, ERV1, ERVK, ERVL and ERVL-MaLR

rm(list=ls()) # remove all entries in the global environment 
set.seed(735) # seed for reproducibility

#-------------------------------------------------------------------------------

              ### Set directory structure and load packages---- 

#-------------------------------------------------------------------------------
setwd("set your working directory here")
rna_dir <- "input folder for RNAseq data"
fig_dir <- "Output folder for plots"

lib_path <- .libPaths("C:/Rpackages/") # specifing the package library
#update.packages(lib.loc = lib_path, ask = FALSE, checkBuilt = TRUE, dependencies = TRUE)
#BiocManager::install(lib.loc = lib_path, update = TRUE, ask = FALSE)

pack_R <- c("dplyr", "tidyverse", "RColorBrewer", "stats", "ggplot2", "ggfortify", 
            "ggpubr", "reshape2", "ggsci", "scales", "pheatmap")

for (i in 1:length(pack_R)) {
  library(pack_R[i], character.only = TRUE)
}

#-------------------------------------------------------------------------------

            ### Read in and format input data---- 

#-------------------------------------------------------------------------------
#### LOUCY
# Read data from csv
rna_loucy <- read.delim(paste0(rna_dir,"Supplementary_Table_17.csv"), sep = "," ) # LOUCY

# Define a function to filter by sample and TE
filter_rna_by_sample <- function(data, sample_id) {
  data %>%
    filter(sample == sample_id, grepl("LTR|LTR\\?", TE))
}

# Apply the function to each sample
ERV_LD3 <- filter_rna_by_sample(rna_loucy, "LD3")
ERV_LG3 <- filter_rna_by_sample(rna_loucy, "LG3")
ERV_LG7 <- filter_rna_by_sample(rna_loucy, "LG7")

# select only the ERVs that are upregulated with GSK day 7
LG7_sig_up <- ERV_LG7[ERV_LG7$logFC > 1 & ERV_LG7$padjust < 0.05,] # upregulated ERVs
LG7_sig_up <- LG7_sig_up$transcript

# Subset to only include significantly upregulated ERVs
ERVs_LD3_sig <- ERV_LD3[ERV_LD3$transcript %in% LG7_sig_up,]
ERVs_LG3_sig <- ERV_LG3[ERV_LG3$transcript %in% LG7_sig_up,]
ERVs_LG7_sig <- ERV_LG7[ERV_LG7$transcript %in% LG7_sig_up,]

# Subset each data frame into the different families
ERVs_LD3_ERV1 <- ERVs_LD3_sig[ERVs_LD3_sig$family=="ERV1",]
ERVs_LD3_ERVK <- ERVs_LD3_sig[ERVs_LD3_sig$family=="ERVK",]
ERVs_LD3_ERVL <- ERVs_LD3_sig[ERVs_LD3_sig$family=="ERVL",]
ERVs_LD3_ERVLMaLR <- ERVs_LD3_sig[ERVs_LD3_sig$family=="ERVL-MaLR",]

ERVs_LG3_ERV1 <- ERVs_LG3_sig[ERVs_LG3_sig$family=="ERV1",]
ERVs_LG3_ERVK <- ERVs_LG3_sig[ERVs_LG3_sig$family=="ERVK",]
ERVs_LG3_ERVL <- ERVs_LG3_sig[ERVs_LG3_sig$family=="ERVL",]
ERVs_LG3_ERVLMaLR <- ERVs_LG3_sig[ERVs_LG3_sig$family=="ERVL-MaLR",]

ERVs_LG7_ERV1 <- ERVs_LG7_sig[ERVs_LG7_sig$family=="ERV1",]
ERVs_LG7_ERVK <- ERVs_LG7_sig[ERVs_LG7_sig$family=="ERVK",]
ERVs_LG7_ERVL <- ERVs_LG7_sig[ERVs_LG7_sig$family=="ERVL",]
ERVs_LG7_ERVLMaLR <- ERVs_LG7_sig[ERVs_LG7_sig$family=="ERVL-MaLR",]

# Format data for pheatmap
# ERV1
ERV_L_ERV1 <- data.frame(cbind(D3_logFC = ERVs_LD3_ERV1$logFC, D3_adjustp = ERVs_LD3_ERV1$padjust,
                               GSK3_logFC = ERVs_LG3_ERV1$logFC, GSK3_adjustp = ERVs_LG3_ERV1$padjust,
                               GSK7_logFC = ERVs_LG7_ERV1$logFC, GSK7_adjustp = ERVs_LG7_ERV1$padjust))

rownames(ERV_L_ERV1) <- ERVs_LG7_ERV1$transcript
ERV_L_ERV1 <- ERV_L_ERV1[order(ERV_L_ERV1$GSK7_adjustp),] # sort based on adjusted pvalues for GSK7
ERV_L_ERV1 <- ERV_L_ERV1[,c(1,3,5)]

# ERVK
ERV_L_ERVK <- data.frame(cbind(D3_logFC = ERVs_LD3_ERVK$logFC, D3_adjustp = ERVs_LD3_ERVK$padjust,
                               GSK3_logFC = ERVs_LG3_ERVK$logFC, GSK3_adjustp = ERVs_LG3_ERVK$padjust,
                               GSK7_logFC = ERVs_LG7_ERVK$logFC, GSK7_adjustp = ERVs_LG7_ERVK$padjust))

rownames(ERV_L_ERVK) <- ERVs_LG7_ERVK$transcript
ERV_L_ERVK <- ERV_L_ERVK[order(ERV_L_ERVK$GSK7_adjustp),] # sort based on adjusted pvalues for GSK7
ERV_L_ERVK <- ERV_L_ERVK[,c(1,3,5)]

# ERVL
ERV_L_ERVL <- data.frame(cbind(D3_logFC = ERVs_LD3_ERVL$logFC, D3_adjustp = ERVs_LD3_ERVL$padjust,
                               GSK3_logFC = ERVs_LG3_ERVL$logFC, GSK3_adjustp = ERVs_LG3_ERVL$padjust,
                               GSK7_logFC = ERVs_LG7_ERVL$logFC, GSK7_adjustp = ERVs_LG7_ERVL$padjust))

rownames(ERV_L_ERVL) <- ERVs_LG7_ERVL$transcript
ERV_L_ERVL <- ERV_L_ERVL[order(ERV_L_ERVL$GSK7_adjustp),] # sort based on adjusted pvalues for GSK7
ERV_L_ERVL <- ERV_L_ERVL[,c(1,3,5)]


# ERVL-MaLR
ERV_L_ERVLMaLR <- data.frame(cbind(D3_logFC = ERVs_LD3_ERVLMaLR$logFC, D3_adjustp = ERVs_LD3_ERVLMaLR$padjust,
                                   GSK3_logFC = ERVs_LG3_ERVLMaLR$logFC, GSK3_adjustp = ERVs_LG3_ERVLMaLR$padjust,
                                   GSK7_logFC = ERVs_LG7_ERVLMaLR$logFC, GSK7_adjustp = ERVs_LG7_ERVLMaLR$padjust))

rownames(ERV_L_ERVLMaLR) <- ERVs_LG7_ERVLMaLR$transcript
ERV_L_ERVLMaLR <- ERV_L_ERVLMaLR[order(ERV_L_ERVLMaLR$GSK7_adjustp),] # sort based on adjusted pvalues for GSK7
ERV_L_ERVLMaLR <- ERV_L_ERVLMaLR[,c(1,3,5)]


#### SUPT1
# Read data from csv
rna_supt1 <- read.delim(paste0(rna_dir,"Supplementary_Table_18.csv"), sep = "," ) # SUPT1

# Filter each sample using the above function
ERV_SD3 <- filter_rna_by_sample(rna_supt1, "SD3")
ERV_SG3 <- filter_rna_by_sample(rna_supt1, "SG3")
ERV_SG7 <- filter_rna_by_sample(rna_supt1, "SG7")

# select only the ERVs that are upregulated with GSK day 7
SG7_sig_up <- ERV_SG7[ERV_SG7$logFC > 1 & ERV_SG7$padjust < 0.05,] # upregulated ERVs
SG7_sig_up <- SG7_sig_up$transcript

# Subset to only include significantly upregulated ERVs
ERVs_SD3_sig <- ERV_SD3[ERV_SD3$transcript %in% SG7_sig_up,]
ERVs_SG3_sig <- ERV_SG3[ERV_SG3$transcript %in% SG7_sig_up,]
ERVs_SG7_sig <- ERV_SG7[ERV_SG7$transcript %in% SG7_sig_up,]

# Subset each data frame into the different families
ERVs_SD3_ERV1 <- ERVs_SD3_sig[ERVs_SD3_sig$family=="ERV1",]
ERVs_SD3_ERVK <- ERVs_SD3_sig[ERVs_SD3_sig$family=="ERVK",]
ERVs_SD3_ERVL <- ERVs_SD3_sig[ERVs_SD3_sig$family=="ERVL",]
ERVs_SD3_ERVLMaLR <- ERVs_SD3_sig[ERVs_SD3_sig$family=="ERVL-MaLR",]

ERVs_SG3_ERV1 <- ERVs_SG3_sig[ERVs_SG3_sig$family=="ERV1",]
ERVs_SG3_ERVK <- ERVs_SG3_sig[ERVs_SG3_sig$family=="ERVK",]
ERVs_SG3_ERVL <- ERVs_SG3_sig[ERVs_SG3_sig$family=="ERVL",]
ERVs_SG3_ERVLMaLR <- ERVs_SG3_sig[ERVs_SG3_sig$family=="ERVL-MaLR",]

ERVs_SG7_ERV1 <- ERVs_SG7_sig[ERVs_SG7_sig$family=="ERV1",]
ERVs_SG7_ERVK <- ERVs_SG7_sig[ERVs_SG7_sig$family=="ERVK",]
ERVs_SG7_ERVL <- ERVs_SG7_sig[ERVs_SG7_sig$family=="ERVL",]
ERVs_SG7_ERVLMaLR <- ERVs_SG7_sig[ERVs_SG7_sig$family=="ERVL-MaLR",]

# Format data for pheatmap
# ERV1
ERV_S_ERV1 <- data.frame(cbind(D3_logFC = ERVs_SD3_ERV1$logFC, D3_adjustp = ERVs_SD3_ERV1$padjust,
                               GSK3_logFC = ERVs_SG3_ERV1$logFC, GSK3_adjustp = ERVs_SG3_ERV1$padjust,
                               GSK7_logFC = ERVs_SG7_ERV1$logFC, GSK7_adjustp = ERVs_SG7_ERV1$padjust))

rownames(ERV_S_ERV1) <- ERVs_SG7_ERV1$transcript
ERV_S_ERV1 <- ERV_S_ERV1[order(ERV_S_ERV1$GSK7_adjustp),] # sort based on adjusted pvalues for GSK7
ERV_S_ERV1 <- ERV_S_ERV1[,c(1,3,5)]

# ERVK
ERV_S_ERVK <- data.frame(cbind(D3_logFC = ERVs_SD3_ERVK$logFC, D3_adjustp = ERVs_SD3_ERVK$padjust,
                               GSK3_logFC = ERVs_SG3_ERVK$logFC, GSK3_adjustp = ERVs_SG3_ERVK$padjust,
                               GSK7_logFC = ERVs_SG7_ERVK$logFC, GSK7_adjustp = ERVs_SG7_ERVK$padjust))

rownames(ERV_S_ERVK) <- ERVs_SG7_ERVK$transcript
ERV_S_ERVK <- ERV_S_ERVK[order(ERV_S_ERVK$GSK7_adjustp),] # sort based on adjusted pvalues for GSK7
ERV_S_ERVK <- ERV_S_ERVK[,c(1,3,5)]

# ERVL
ERV_S_ERVL <- data.frame(cbind(D3_logFC = ERVs_SD3_ERVL$logFC, D3_adjustp = ERVs_SD3_ERVL$padjust,
                               GSK3_logFC = ERVs_SG3_ERVL$logFC, GSK3_adjustp = ERVs_SG3_ERVL$padjust,
                               GSK7_logFC = ERVs_SG7_ERVL$logFC, GSK7_adjustp = ERVs_SG7_ERVL$padjust))

rownames(ERV_S_ERVL) <- ERVs_SG7_ERVL$transcript
ERV_S_ERVL <- ERV_S_ERVL[order(ERV_S_ERVL$GSK7_adjustp),] # sort based on adjusted pvalues for GSK7
ERV_S_ERVL <- ERV_S_ERVL[,c(1,3,5)]


# ERVL-MaLR
ERV_S_ERVLMaLR <- data.frame(cbind(D3_logFC = ERVs_SD3_ERVLMaLR$logFC, D3_adjustp = ERVs_SD3_ERVLMaLR$padjust,
                                   GSK3_logFC = ERVs_SG3_ERVLMaLR$logFC, GSK3_adjustp = ERVs_SG3_ERVLMaLR$padjust,
                                   GSK7_logFC = ERVs_SG7_ERVLMaLR$logFC, GSK7_adjustp = ERVs_SG7_ERVLMaLR$padjust))

rownames(ERV_S_ERVLMaLR) <- ERVs_SG7_ERVLMaLR$transcript
ERV_S_ERVLMaLR <- ERV_S_ERVLMaLR[order(ERV_S_ERVLMaLR$GSK7_adjustp),] # sort based on adjusted pvalues for GSK7
ERV_S_ERVLMaLR <- ERV_S_ERVLMaLR[,c(1,3,5)]

#-------------------------------------------------------------------------------

                          ### Figure 5HI---- 

#-------------------------------------------------------------------------------
# Setting parameters for plotting
n.neg <- 4
n.pos <- 8
breakpoints <- c(-10, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 14, 16)
colors <- c(rev(RColorBrewer::brewer.pal(n.neg-1, "Blues")), "white", RColorBrewer::brewer.pal(n.pos-1, "Reds"))


ERV1_L_plot <- pheatmap(ERV_L_ERV1, cluster_cols = F, cluster_rows = F, color = colors, breaks = breakpoints, cellwidth = 30, cellheight = 1.2) # export plot
ERVK_L_plot <- pheatmap(ERV_L_ERVK, cluster_cols = F, cluster_rows = F, color = colors, breaks = breakpoints, cellwidth = 30, cellheight = 1.2) # export plot
ERVL_L_plot <- pheatmap(ERV_L_ERVL, cluster_cols = F, cluster_rows = F, color = colors, breaks = breakpoints, cellwidth = 30, cellheight = 1.2) # export plot
ERVLMaLR_L_plot <- pheatmap(ERV_L_ERVLMaLR, cluster_cols = F, cluster_rows = F, color = colors, breaks = breakpoints, cellwidth = 30, cellheight = 1.2) # export plot

ERV1_S_plot <- pheatmap(ERV_S_ERV1, cluster_cols = F, cluster_rows = F, color = colors, breaks = breakpoints, cellwidth = 30, cellheight = 1.2) # export plot
ERVK_S_plot <- pheatmap(ERV_S_ERVK, cluster_cols = F, cluster_rows = F, color = colors, breaks = breakpoints, cellwidth = 30, cellheight = 1.2) # export plot
ERVL_S_plot <- pheatmap(ERV_S_ERVL, cluster_cols = F, cluster_rows = F, color = colors, breaks = breakpoints, cellwidth = 30, cellheight = 1.2) # export plot
ERVLMaLR_S_plot <- pheatmap(ERV_S_ERVLMaLR, cluster_cols = F, cluster_rows = F, color = colors, breaks = breakpoints, cellwidth = 30, cellheight = 1.2) # export plot

