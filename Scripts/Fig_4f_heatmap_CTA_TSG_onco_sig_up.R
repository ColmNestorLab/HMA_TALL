################################################################################
###### HEATMAP: expression of differentially expressed genes of interest ####### 
################################################################################

# Maike Bensberg (maike.bensberg@liu.se)
# 2024-08-22


### plot heatmap of differential expression of genes of interest
# Cancer/Testis antigens (CTAs), tumour suppressor genes (TSGs) and oncogenes
# plot all treatments for both cell lines for all gene that are upregulated in any sample
# star for significant differential expression needs to be added manually
# grey square for samples in which data is missing


### CTAs, TSGs and oncogenes summarized in Supplementary Table 14
# immune genes also included in list, but will not be plotted here
# TSGs and oncogenes based on the Network of cancer genes & healthy drivers 
# http://network-cancer-genes.org/citation.php
# CTAs based on publication by Carter et al, 2023
# https://jitc.bmj.com/content/11/12/e007935


### Figure: Fig. 4f


### enrichment analysis for same gene groups done (script enrichment_phyper.R)



##### initiate R environment #####
# only needs to be done once
#install.packages("renv")
#renv::init()



##### loading packages #####

library(pheatmap)
library(reshape2)



##### set working directory and check package versions #####

# set working directory
setwd("/mnt/wwn-0x5000cca28fd30c7d-part1/HMA_TALL/EMseq/LOUCY_SUPT1_HMAtreated_12_2023/methylKit")
# print working directory
getwd()

# print version of R and attached packages
sessionInfo()



##### import data #####


### differentially expressed genes
# Supplementary Table 12
df_expression_combined <- read.delim("R_exports/tables/Supplementary_Table_12.csv", sep = ",")

### DNA methylation sensitive genes
# Supplementary Table 13
methylation_expression_N_transcripts <- read.delim("R_exports/tables/Supplementary_Table_13.csv", sep = ",")

### lists of genes
# Supplementary Table 14
genes_of_interest <- as.data.frame(read_excel("R_imports/Supplementary_Table_14.xlsx", 
                                              skip = 1))



##### make dataframe for plotting #####
# separate dataframe for immune genes

df_for_plotting <-
  df_expression_combined[which(df_expression_combined$gene_id %in% 
                                 genes_of_interest[which(genes_of_interest$Category != "Immune"),"Gene"]),
                         c("gene_id", "logFC", "padjust", "sample")]

df_for_plotting_immune <-
  df_expression_combined[which(df_expression_combined$gene_id %in% 
                                 genes_of_interest[which(genes_of_interest$Category == "Immune"),"Gene"]),
                         c("gene_id", "logFC", "padjust", "sample")]

# find genes that are significantly upregulated in any sample
# log2FC > 1 and p-adjusted < 0.05

genes_up <- unique(df_for_plotting[which(df_for_plotting$logFC > 1 &
                                           df_for_plotting$padjust < 0.05),"gene_id"])

df_for_plotting <- df_for_plotting[which(df_for_plotting$gene_id %in% genes_up),]

genes_up_immune <- unique(df_for_plotting_immune[which(df_for_plotting_immune$logFC > 1 &
                                           df_for_plotting_immune$padjust < 0.05),"gene_id"])

df_for_plotting_immune <- df_for_plotting_immune[which(df_for_plotting_immune$gene_id %in% genes_up_immune),]

# melt dataframe to matrix
matrix_for_plotting <-
  dcast(df_for_plotting[,c("gene_id", "logFC", "sample")], gene_id ~ sample, value.var = "logFC")

matrix_for_plotting_immune <-
  dcast(df_for_plotting_immune[,c("gene_id", "logFC", "sample")], gene_id ~ sample, value.var = "logFC")



##### make separate dataframe for annotation ##### 

annotation <- matrix_for_plotting
annotation_immune <- matrix_for_plotting_immune

# define category
# if CTA overlaps with either of the other categories it is classified as TSG/Oncogene
annotation$category <- "CTA"
annotation[which(annotation$gene_id %in%
                        genes_of_interest[which(genes_of_interest$Category == "TSG"),"Gene"]),
                "category"] <- "TSG"
annotation[which(annotation$gene_id %in%
                        genes_of_interest[which(genes_of_interest$Category == "onco"),"Gene"]),
                "category"] <- "onco"

annotation_immune$category <- "immune"

# add column for methylation sensitivity after treatment with GSK for 7 days
LG7_methyl_sensitive_genes <- unique(
  methylation_expression_N_transcripts[
    which(methylation_expression_N_transcripts$methyl_sensitive == "yes" &
            methylation_expression_N_transcripts$sample == "LG7"), "gene_id"])

SG7_methyl_sensitive_genes <- unique(
  methylation_expression_N_transcripts[
    which(methylation_expression_N_transcripts$methyl_sensitive == "yes" &
            methylation_expression_N_transcripts$sample == "SG7"), "gene_id"])

annotation$methyl_sensitive <- "NA"
annotation[which(annotation$gene_id %in%
                   unique(methylation_expression_N_transcripts[
                     which(methylation_expression_N_transcripts$methyl_sensitive == "no"), "gene_id"])),
           "methyl_sensitive"] <- "no"
annotation[which(annotation$gene_id %in% LG7_methyl_sensitive_genes),
           "methyl_sensitive"] <- "LOUCY"
annotation[which(annotation$gene_id %in% SG7_methyl_sensitive_genes),
           "methyl_sensitive"] <- "SUPT1"
annotation[which(annotation$gene_id %in% LG7_methyl_sensitive_genes &
                   annotation$gene_id %in% SG7_methyl_sensitive_genes),
           "methyl_sensitive"] <- "both"

annotation_immune$methyl_sensitive <- "NA"
annotation_immune[which(annotation_immune$gene_id %in%
                          unique(methylation_expression_N_transcripts[
                            which(methylation_expression_N_transcripts$methyl_sensitive == "no"), "gene_id"])),
                  "methyl_sensitive"] <- "no"
annotation_immune[which(annotation_immune$gene_id %in% LG7_methyl_sensitive_genes),
                  "methyl_sensitive"] <- "LOUCY"
annotation_immune[which(annotation_immune$gene_id %in% SG7_methyl_sensitive_genes),
                  "methyl_sensitive"] <- "SUPT1"
annotation_immune[which(annotation_immune$gene_id %in% LG7_methyl_sensitive_genes &
                          annotation_immune$gene_id %in% SG7_methyl_sensitive_genes),
                  "methyl_sensitive"] <- "both"
  


# define gene_id as row names
rownames(matrix_for_plotting) <- matrix_for_plotting[,"gene_id"] 
rownames(annotation) <- annotation[,"gene_id"] 

rownames(matrix_for_plotting_immune) <- matrix_for_plotting_immune[,"gene_id"] 
rownames(annotation_immune) <- annotation_immune[,"gene_id"] 

# remove unneccesary columns
matrix_for_plotting <- matrix_for_plotting[,-1]
annotation <- annotation[,-c(1:7)]

matrix_for_plotting_immune <- matrix_for_plotting_immune[,-1]
annotation_immune <- annotation_immune[,-c(1:7)]

# order by category
matrix_for_plotting <- matrix_for_plotting[order(annotation$category),]




##### plot heatmap #####

# define colours
n.neg <- 4
n.pos <- 8
breakpoints <- c(seq(-5, -1, length.out=n.neg), seq(1, 15, length.out=n.pos))
colors <- c(rev(RColorBrewer::brewer.pal(n.neg-1, "Blues")), "white", RColorBrewer::brewer.pal(n.pos-1, "Reds"))


### plot

# separate heatmaps for immune genes because some immune genes and oncogenes overlap
pheatmap(matrix_for_plotting, cluster_cols = F, cluster_rows = F,
         breaks = breakpoints, color = colors, annotation_row = annotation)

pheatmap(matrix_for_plotting_immune, cluster_cols = F, cluster_rows = F,
         breaks = breakpoints, color = colors, annotation_row = annotation_immune)


remove(matrix_for_plotting)
remove(annotation)
remove(annotation_immune)
remove(LG7_methyl_sensitive_genes)
remove(SG7_methyl_sensitive_genes)
remove(n.neg)
remove(n.pos)
remove(breakpoints)
remove(colors)

### export heatmap manually



##### find significance #####

### stars for significant differential expression need to be added manually to heatmap
# significant fold change: adjusted p < 0.05

# melt dataframe to matrix
df_for_plotting$sig <- "no"
df_for_plotting[which(df_for_plotting$padjust < 0.05),"sig"] <- "yes"

df_for_plotting_immune$sig <- "no"
df_for_plotting_immune[which(df_for_plotting_immune$padjust < 0.05),"sig"] <- "yes"

# show only samples with significant differential expression
df_for_plotting[which(df_for_plotting$padjust < 0.05 &
                        df_for_plotting$sample == "SG7"),]
nrow(df_for_plotting[which(df_for_plotting$padjust < 0.05 &
                             df_for_plotting$sample == "SG7"),])

df_for_plotting_immune[which(df_for_plotting_immune$padjust < 0.05 &
                               df_for_plotting_immune$sample == "SG7"),]
nrow(df_for_plotting_immune[which(df_for_plotting_immune$padjust < 0.05 &
                                    df_for_plotting_immune$sample == "SG7"),])


remove(df_for_plotting)
remove(df_for_plotting_immune)