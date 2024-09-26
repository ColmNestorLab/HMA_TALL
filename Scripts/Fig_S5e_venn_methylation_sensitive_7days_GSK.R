################################################################################
#################### VENN: methyl-sensitive LOUCY & SUP-T1 ##################### 
################################################################################

# Maike Bensberg (maike.bensberg@liu.se)
# 2024-08-28


### plot a Venn diagramm with overlap between methylation sensitive genes
# only methylation sensitive genes after traetment with GSK for 7 days
# in LOUCY and SUP_T1 cells

### methylation sensitive genes defined based on DEG and DMR analysis
# see script: combining_DNAmethylation_and_expression_data.R
# Supplementary Table 13


### Figure: Extended Data Fig. 5e



##### initiate R environment #####
# only needs to be done once
#install.packages("renv")
#renv::init()



##### loading packages #####

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

### Supplementary Table 13
# data on methylation and expression was combined and exported in script 
#   combining_DNAmethylation_and_expression_data.R

### methylation sensitive genes defined as
# silent in untreated control cell (cpm < 0.5)
# expressed and upregulated in treated cells (cpm >= 0.5, logFC > 1 and padjust < 0.05)
# promoter gets de-methylated after treatment (meth.diff <= -25 and q < 0.01)

methylation_expression_N_transcripts <- read.delim("R_exports/tables/Supplementary_Table_13.csv", sep = ",")


### extract gene lists
# methylation sensitive genes for LOUCY and SUP-T1 cells
# after treatment with GSK for 7 days

methyl_sensitive_LG7 <- unique(
  methylation_expression_N_transcripts[
    which(methylation_expression_N_transcripts$sample == "LG7" &
            methylation_expression_N_transcripts$methyl_sensitive == "yes"),"gene_id"])

methyl_sensitive_SG7 <- unique(
  methylation_expression_N_transcripts[
    which(methylation_expression_N_transcripts$sample == "SG7" &
            methylation_expression_N_transcripts$methyl_sensitive == "yes"),"gene_id"])



##### plot #####

# colours from pre-defined colour paletts
# LOUCY: nejm[l] (blue-grey), SUP-T1: npg[5] (orange)

plot <-
  plot(euler(list("LOUCY" = methyl_sensitive_LG7,
                  "SUP_T1"  = methyl_sensitive_SG7)),
       quantities = TRUE,
       fills = list(fill = c(nejm_colors [6], npg_colors [5]), alpha = 0.5),
       main = "methylation sensitive genes")


### export

pdf("R_exports/plots/Fig_S5e_venn_methylation_sensitive_7days_GSK.pdf")
plot(plot)
dev.off()


remove(plot)
remove(methyl_sensitive_LG7)
remove(methyl_sensitive_SG7)