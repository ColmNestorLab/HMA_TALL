################################################################################
################### VENN: promoters keeping DNA methylation ################### 
################################################################################

# Maike Bensberg (maike.bensberg@liu.se)
# 2024-08-28


### plot a Venn diagram with overlap between promoters keeping DNA methylation
# promoters of which genes overlap with regions keeping DNA methylation
# promoter: 1000 bp up- and 500 bp downstream of the TSS
# regions keeping DBA methylation defined as having >=80% DNA methylation before and after treatment with GSK for 7 days
# comparing LOUCY and SUP-T1 cells
# see script: regions_retaining_methylation_always_80.R


### Figure: Extended Data Fig. 5a



##### initiate R environment #####
# only needs to be done once
#install.packages("renv")
#renv::init()



##### loading packages #####

library(readxl)
library(ggsci)
library(scales)
library(ggplot2)
library(eulerr)



##### set working directory and check package versions #####

# set working directory
setwd("/mnt/wwn-0x5000cca28fd30c7d-part1/HMA_TALL/EMseq/LOUCY_SUPT1_HMAtreated_12_2023/methylKit")
# print working directory
getwd()

# print version of R and attached packages
sessionInfo()



##### defining color palettes #####

npg_colors <- pal_npg()(10)
show_col (npg_colors)
nejm_colors <- pal_nejm()(8)
show_col (nejm_colors)


##### import and format data #####

### gene IDs 
# promoters with >= 80% DNA methylation before and after treatment with GSK for 7 days

promoters_retain_L_gene_ids <- read.delim("R_exports/tables/promoters_retain_L_gene_ids.txt",
                                          header = FALSE, col.names = "gene_id")
promoters_retain_S_gene_ids <- read.delim("R_exports/tables/promoters_retain_S_gene_ids.txt",
                                          header = FALSE, col.names = "gene_id")



##### plot #####

# colours from pre-defined colour paletts
# LOUCY: nejm[l] (blue-grey), SUP-T1: npg[5] (orange)

plot <-
  plot(euler(list("LOUCY" = promoters_retain_L_gene_ids$gene_id,
                  "SUP_T1"  = promoters_retain_S_gene_ids$gene_id)),
       quantities = TRUE,
       fills = list(fill = c(nejm_colors [6], npg_colors [5]), alpha = 0.5),
       main = "genes keeping DNA methylation")


### export

pdf("R_exports/plots/Fig_S5a_venn_retain_promoter_methylation_80.pdf")
plot(plot)
dev.off()


remove(plot)
remove(promoters_retain_L_gene_ids)
remove(promoters_retain_S_gene_ids)