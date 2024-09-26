################################################################################
################ Combining DNA methylation and expression data ################# 
################################################################################

# Maike Bensberg (maike.bensberg@liu.se)
# 2024-08-13


### data on DNA methylation and expression
# DNA methylation and differentially methylated promoters
# gene expression, including differential gene expression, analysed by Sandra H.


### combining data
# promoter DNA methylation for individual transcripts
# gene expression summarized by gene (not transcript)
# combine based on gene ID


### exported data used in Fig. 4 and Extended Data Fig. 5-8
# Supplementary Tables 12 and 13



##### initiate R environment #####
# only needs to be done once
#install.packages("renv")
#renv::init()



##### loading packages #####

library(readr)
library(tibble)
library(dplyr)
library(reshape2)



##### set working directory and check package versions #####

# set working directory
setwd("/mnt/wwn-0x5000cca28fd30c7d-part1/HMA_TALL/EMseq/LOUCY_SUPT1_HMAtreated_12_2023/methylKit")
# print working directory
getwd()

# print version of R and attached packages
sessionInfo()



##### import and format data #####


### DNA methylation at promoters
# analysis using methylKit
# processed using script: promoter_methylation_analysis.R
# q-value: Benjamini Hochberg adjusted p-value
# annotated with gene name and transcript ID
# includes DNA methylation for treated and untreated samples
# sample column defines treatment

# Supplementary Table 8
df_DMR_promoters <- read.delim("R_exports/tables/Supplementary_Table_8.csv", sep = ",")


### Counts for RNA-seq samples
# separated by cell line
# as counts per million (cpm) and log(cpm)
# not expressed genes are removed (combined raw read count of > 15 and > 10 in at least one sample)
# counts normalized

# cpm
LOUCY_RNA_cpm <- read.delim("R_imports/count_data_Genes/LOUCY__RNAseq_cpm_genes_240806.txt")
LOUCY_RNA_cpm$gene_id <- rownames(LOUCY_RNA_cpm)
rownames(LOUCY_RNA_cpm) <- NULL
colnames(LOUCY_RNA_cpm) [1:4] <- c("Lctrl7", "LD3", "LG3", "LG7")
SUPT1_RNA_cpm <- read.delim("R_imports/count_data_Genes/SUPT1__RNAseq_cpm_genes_240806.txt")
SUPT1_RNA_cpm$gene_id <- rownames(SUPT1_RNA_cpm)
rownames(SUPT1_RNA_cpm) <- NULL
colnames(SUPT1_RNA_cpm) [1:4] <- c("Sctrl7", "SD3", "SG3", "SG7")

# log2cpm
LOUCY_RNA_log2cpm <- read.delim("R_imports/count_data_Genes/LOUCY__RNAseq_lcpm_genes_240806.txt")
LOUCY_RNA_log2cpm$gene_id <- rownames(LOUCY_RNA_log2cpm)
rownames(LOUCY_RNA_log2cpm) <- NULL
colnames(LOUCY_RNA_log2cpm) [1:4] <- c("Lctrl7", "LD3", "LG3", "LG7")
SUPT1_RNA_log2cpm <- read.delim("R_imports/count_data_Genes/SUPT1__RNAseq_lcpm_genes_240806.txt")
SUPT1_RNA_log2cpm$gene_id <- rownames(SUPT1_RNA_log2cpm)
rownames(SUPT1_RNA_log2cpm) <- NULL
colnames(SUPT1_RNA_log2cpm) [1:4] <- c("Sctrl7", "SD3", "SG3", "SG7")

# melt dataframes for cpm and log(cpm)
LOUCY_RNA_cpm <- melt(LOUCY_RNA_cpm, id.vars = "gene_id",
                      variable.name = "sample", value.name = "cpm")
SUPT1_RNA_cpm <- melt(SUPT1_RNA_cpm, id.vars = "gene_id",
                      variable.name = "sample", value.name = "cpm")
LOUCY_RNA_log2cpm <- melt(LOUCY_RNA_log2cpm, id.vars = "gene_id",
                          variable.name = "sample", value.name = "log2cpm")
SUPT1_RNA_log2cpm <- melt(SUPT1_RNA_log2cpm, id.vars = "gene_id",
                          variable.name = "sample", value.name = "log2cpm")


### differential gene expression analysis
# separate file for each treatment vs. control comparison
# normalized with expected dispersion set to 0.4
# p-value adjusted for multiple testing using Benjamini Hochberg

# import
DGE_file_names <- as.list(list.files(path = "R_imports/count_data_Genes",
                                    pattern = "_RNAseq_exactTest_genes_bcv0.4_240806.txt",
                                    full.names = TRUE))

DGE_RNAseq_files <- lapply(DGE_file_names, read.delim)

# make row names into new column for gene_id
DGE_RNAseq_files <- lapply(DGE_RNAseq_files, rownames_to_column, var = "gene_id")

names <- gsub(list.files(path = "R_imports/count_data_Genes",
                         pattern = "_RNAseq_exactTest_genes_bcv0.4_240806.txt"),
              pattern = "_RNAseq_exactTest_genes_bcv0.4_240806.txt", replacement = "")

for (i in c(1:length(names))) {
  assign(names[i],DGE_RNAseq_files[[i]])
}

remove(i)
remove(names)
remove(DGE_file_names)
remove(DGE_RNAseq_files)


### combine differential gene expression by cell line
# define sample in new column
df_DGE_L <- do.call("rbind", list(LOUCY_DECvsCtrl, LOUCY_GSK3vsCtrl, LOUCY_GSK7vsCtrl))
df_DGE_L$sample <- c(rep("LD3", nrow(LOUCY_DECvsCtrl)),
                     rep("LG3", nrow(LOUCY_GSK3vsCtrl)),
                     rep("LG7", nrow(LOUCY_GSK7vsCtrl)))

df_DGE_S <- do.call("rbind", list(SUPT1_DECvsCtrl, SUPT1_GSK3vsCtrl, SUPT1_GSK7vsCtrl))
df_DGE_S$sample <- c(rep("SD3", nrow(SUPT1_DECvsCtrl)),
                     rep("SG3", nrow(SUPT1_GSK3vsCtrl)),
                     rep("SG7", nrow(SUPT1_GSK7vsCtrl)))


### combine DGE table for all samples with cpm and log2cpm values

# add cpm and log2cpm for treated samples
df_expression_L <- inner_join(df_DGE_L,
                              inner_join(LOUCY_RNA_cpm, LOUCY_RNA_log2cpm))
colnames(df_expression_L) [7:8] <- c("cpm_treated", "log2cpm_treated")

df_expression_S <- inner_join(df_DGE_S,
                              inner_join(SUPT1_RNA_cpm, SUPT1_RNA_log2cpm))
colnames(df_expression_S) [7:8] <- c("cpm_treated", "log2cpm_treated")

# making a temporary dataframes for cpm and log2cpm of untreated samples
df_temp_L <- inner_join(LOUCY_RNA_cpm, LOUCY_RNA_log2cpm)
df_temp_L <- df_temp_L[grep(pattern = "Lctrl7", x = df_temp_L$sample),c(1,3,4)]

df_temp_S <- inner_join(SUPT1_RNA_cpm, SUPT1_RNA_log2cpm)
df_temp_S <- df_temp_S[grep(pattern = "Sctrl7", x = df_temp_S$sample),c(1,3,4)]

# adding cpm and log2cpm for untreated samples
df_expression_L <- inner_join(df_expression_L, df_temp_L)
colnames(df_expression_L) [9:10] <- c("cpm_control", "log2cpm_control")

df_expression_S <- inner_join(df_expression_S, df_temp_S)
colnames(df_expression_S) [9:10] <- c("cpm_control", "log2cpm_control")

# removing temporary dataframes
remove(df_temp_L)
remove(df_temp_S)


# combine cell lines
df_expression_combined <- rbind(df_expression_L, df_expression_S)


### add summary of expression changes to combined DEG dataframe
# upregulated: logFC > 1 and adjusted p < 0.05
# downregulated: logFC < -1 and adjusted p < 0.05
# no change: all others

df_expression_combined$expression <- "no change"
df_expression_combined[
  which(df_expression_combined$padjust < 0.05 &
          df_expression_combined$logFC > 1),]$expression <- "upregulated"
df_expression_combined[
  which(df_expression_combined$padjust < 0.05 &
          df_expression_combined$logFC < -1),]$expression <- "downregulated"




##### combine DNA methylation and expression data #####


### combine differential expression analysis and methylation analysis
# combine tables based on gene ID and sample

df_methylation_expression_N_transcripts <-
  inner_join(df_DMR_promoters, df_expression_combined[,c(1,2,4:11)])


### add columns interpreting changes in expression and DNA methylation


### DNA methylation
# gain: methylation difference >= 25 and adjusted p < 0.01
# loss: methylation difference <= -25 and adjusted p < 0.01
# no change: all others

df_methylation_expression_N_transcripts$methylation <- "no change"
df_methylation_expression_N_transcripts[
  which(df_methylation_expression_N_transcripts$qvalue < 0.01 &
          df_methylation_expression_N_transcripts$meth.diff >= 25),]$methylation <- "gain"
df_methylation_expression_N_transcripts[
  which(df_methylation_expression_N_transcripts$qvalue < 0.01 &
          df_methylation_expression_N_transcripts$meth.diff <= -25),]$methylation <- "loss"


### methylation sensitive genes
# silent in untreated control cell (cpm < 0.5)
# expressed and upregulated in treated cells (cpm >= 0.5, logFC > 1 and padjust < 0.05)
# promoter gets de-methylated after treatment (meth.diff <= -25 and q < 0.01)

df_methylation_expression_N_transcripts$methyl_sensitive <- "no"
df_methylation_expression_N_transcripts[
  which(df_methylation_expression_N_transcripts$cpm_control < 0.5 &
          df_methylation_expression_N_transcripts$cpm_treated >= 0.5 &
          df_methylation_expression_N_transcripts$logFC > 1 &
          df_methylation_expression_N_transcripts$padjust < 0.05 &
          df_methylation_expression_N_transcripts$meth.diff <= -25 &
          df_methylation_expression_N_transcripts$qvalue < 0.01),]$methyl_sensitive <- "yes"



##### export data #####
# write as comma separated files

### expression only
# will be supplementary table 12
write.table(df_expression_combined, "R_exports/tables/Supplementary_Table_12.csv",
            sep = ",", quote = FALSE, col.names = TRUE, row.names = FALSE)

### DNA methylation and expression
# Supplementary Table 13
write.table(df_methylation_expression_N_transcripts, "R_exports/tables/Supplementary_Table_13.csv",
            sep = ",", quote = FALSE, col.names = TRUE, row.names = FALSE)