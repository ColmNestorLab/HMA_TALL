################################################################################
################## CIRCOS PLOTS: whole genome DNA methylation ##################
################################################################################

# Maike Bensberg (maike.bensberg@liu.se)
# 2024-08-12


### make a circos plot showing whole genome DNA methylation
# based on 500 kb bins
# script for pre-processing the data: methylKit_500kb_bins.R

### Figure: Fig. 4a


### using Circlize
# following this guide: https://jokergoo.github.io/circlize_book/book/circos-heatmap.html


##### initiate R environment #####
# only needs to be done once
#install.packages("renv")
#renv::init()



##### loading packages #####

library(circlize)
library(readr)



##### set working directory and check package versions #####

# set working directory
setwd("/mnt/wwn-0x5000cca28fd30c7d-part1/HMA_TALL/EMseq/LOUCY_SUPT1_HMAtreated_12_2023/methylKit")
# print working directory
getwd()

# print version of R and attached packages
sessionInfo()



##### import data #####

### import DNA methylation data 
# average DNA methylation in 100 kb and 500 kb bins (not overlapping)
# scripts for data processing: methylKit_100kb_bins.R and methylKit_500kb_bins.R

df_EMseq_bins500kb_L_merged <- read.delim("R_exports/tables/df_EMseq_bins500kb_L_merged.txt",
                                          sep = "\t", header = TRUE)
df_EMseq_bins500kb_S_merged <- read.delim("R_exports/tables/df_EMseq_bins500kb_S_merged.txt",
                                          sep = "\t", header = TRUE)



##### plot circos plot #####

### define colours

col_fun1 <- colorRamp2(c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
                       c("#071630" ,"#1984c5",  "#a7d5ed", "#66C2A5","#ABDDA4" ,"#E6F598",
                         "#FFFFBF","#FDAE61","#F46D43","#D53E4F", "#9E0142"))


### define chromosomes to include
# exclude chr Y and any random/undefined chromosomes

chr <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
         "chr9", "chr10", "chr11", "chr12", "chr13", "chr14","chr15", "chr16",
         "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")


### plot LOUCY 500 kb bins

split <- df_EMseq_bins500kb_L_merged [which(
  df_EMseq_bins500kb_L_merged$chr %in% chr),"chr"]
split <- factor(split, levels = chr)

circos.clear()
circos.par(gap.after = c(rep(2,22),30))
# GSK 7 days
circos.heatmap(mat = df_EMseq_bins500kb_L_merged[which(
  df_EMseq_bins500kb_L_merged$chr %in% chr),"LG7"],
  split = split, col = col_fun1, 
  show.sector.labels = TRUE, cluster = FALSE,
  track.height = 0.1)
# GSK 3 days
circos.heatmap(mat = df_EMseq_bins500kb_L_merged[which(
  df_EMseq_bins500kb_L_merged$chr %in% chr),"LG3"],
  split = split, col = col_fun1, 
  show.sector.labels = TRUE, cluster = FALSE,
  track.height = 0.1)
# DEC 3 days
circos.heatmap(mat = df_EMseq_bins500kb_L_merged[which(
  df_EMseq_bins500kb_L_merged$chr %in% chr),"LD3"],
  split = split, col = col_fun1, 
  show.sector.labels = TRUE, cluster = FALSE,
  track.height = 0.1)
# control
circos.heatmap(mat = df_EMseq_bins500kb_L_merged[which(
  df_EMseq_bins500kb_L_merged$chr %in% chr),"Lctrl7"],
  split = split, col = col_fun1, 
  show.sector.labels = TRUE, cluster = FALSE,
  track.height = 0.1)

remove(split)

# export plots as pdfs manually


### plot SUP-T1 500 kb bins

split <- df_EMseq_bins500kb_S_merged [which(
  df_EMseq_bins500kb_S_merged$chr %in% chr),"chr"]
split <- factor(split, levels = chr)

circos.clear()
circos.par(gap.after = c(rep(2,22),30))
# GSK 7 days
circos.heatmap(mat = df_EMseq_bins500kb_S_merged[which(
  df_EMseq_bins500kb_S_merged$chr %in% chr),"SG7"],
  split = split, col = col_fun1, 
  show.sector.labels = TRUE, cluster = FALSE,
  track.height = 0.1)
# GSK 3 days
circos.heatmap(mat = df_EMseq_bins500kb_S_merged[which(
  df_EMseq_bins500kb_S_merged$chr %in% chr),"SG3"],
  split = split, col = col_fun1, 
  show.sector.labels = TRUE, cluster = FALSE,
  track.height = 0.1)
# DEC 3 days
circos.heatmap(mat = df_EMseq_bins500kb_S_merged[which(
  df_EMseq_bins500kb_S_merged$chr %in% chr),"SD3"],
  split = split, col = col_fun1, 
  show.sector.labels = TRUE, cluster = FALSE,
  track.height = 0.1)
# control
circos.heatmap(mat = df_EMseq_bins500kb_S_merged[which(
  df_EMseq_bins500kb_S_merged$chr %in% chr),"Sctrl7"],
  split = split, col = col_fun1, 
  show.sector.labels = TRUE, cluster = FALSE,
  track.height = 0.1)

remove(split)
remove(chr)
remove(col_fun1)

# export plots as pdfs manually