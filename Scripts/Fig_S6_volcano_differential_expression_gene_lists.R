################################################################################
######### VOLCANO PLOTS: differentially expressed genes from gene sets ######### 
################################################################################

# Maike Bensberg (maike.bensberg@liu.se)
# 2024-08-23


### plot volcano plots of differntially expressed CTAs, TSGs, oncogenes and immune genes
# CTA: Cancer/testis antigen, TSG: tumour suppressor genes
# list of immune genes curated by Sandra H.
# Supplementary Table 12


### Figure: Extended data Fig 6



##### initiate R environment #####
# only needs to be done once
#install.packages("renv")
#renv::init()



##### loading packages #####

library(ggsci)
library(scales)
library(ggplot2)
library(ggpubr)



##### set working directory and check package versions #####

# set working directory
setwd("/mnt/wwn-0x5000cca28fd30c7d-part1/HMA_TALL/EMseq/LOUCY_SUPT1_HMAtreated_12_2023/methylKit")
# print working directory
getwd()

# print version of R and attached packages
sessionInfo()



##### import and format data #####


### import

### differentially expressed genes
# summary dataframe generated using script combining_DNAmethylation_and_expression_data.R
# Supplementary Table 10
df_expression_combined <- read.delim("R_exports/tables/Supplementary_Table_10.csv", sep = ",")


### lists of genes
# Supplementary Table 12
genes_of_interest <- as.data.frame(read_excel("R_imports/Supplementary_Table_12.xlsx", 
                                              skip = 1))


### format

# make dataframe for plotting
# make one dataframe for genes that are CTAs, TSGs or oncogenes
# second dataframe for immune genes (because some immune genes are also TSGs/oncogenes)
# only keep information on gene, logFC, adjusted p-value and sample
df_for_plotting <-
  df_expression_combined[which(df_expression_combined$gene_id %in% 
                                 genes_of_interest[which(genes_of_interest$Category != "Immune"),"Gene"]),
                         c("gene_id", "logFC", "padjust", "sample")]
df_for_plotting_immune <-
  df_expression_combined[which(df_expression_combined$gene_id %in% 
                                 genes_of_interest[which(genes_of_interest$Category == "Immune"),"Gene"]),
                         c("gene_id", "logFC", "padjust", "sample")]

# annotate with gene category
df_for_plotting$category <- "CTA"
df_for_plotting[which(df_for_plotting$gene_id %in%
                   genes_of_interest[which(genes_of_interest$Category == "TSG"),"Gene"]),
           "category"] <- "TSG"
df_for_plotting[which(df_for_plotting$gene_id %in%
                   genes_of_interest[which(genes_of_interest$Category == "onco"),"Gene"]),
           "category"] <- "onco"

df_for_plotting_immune$category <- "immune"

# add column indicating significance
# differential expression is considered significant for adjusted p < 0.05
# logFC > 1 or < 1
df_for_plotting$sig <-"no"
df_for_plotting[which(df_for_plotting$padjust < 0.05 &
                        (df_for_plotting$logFC < -1 |
                           df_for_plotting$logFC > 1)),"sig"] <- "yes"
df_for_plotting_immune$sig <-"no"
df_for_plotting_immune[which(df_for_plotting_immune$padjust < 0.05 &
                        (df_for_plotting_immune$logFC < -1 |
                           df_for_plotting_immune$logFC > 1)),"sig"] <- "yes"



##### plot #####

# volcano plots (x: FC and y: p-value)

plot <- ggarrange(
  ggplot(data = df_for_plotting, aes(x = logFC, y = -log10(padjust), fill = sig, label = gene_id)) +
  geom_point(shape = 21, colour = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 1, color = "black") +
  geom_vline(xintercept = -1, color = "black") +
  scale_fill_manual(values = c("yes" = "#8B0000", "no" = "grey")) +
  geom_text(data = subset(df_for_plotting, sig == "yes"), hjust = -0.1, size = 3) +
  facet_grid(sample~category) +
  scale_y_continuous(limits = c(0,11)) +
  labs(title = "Differential expression of CTAs, TSGs and oncogenes") +
  theme_bw (),
  ggplot(data = df_for_plotting_immune, aes(x = logFC, y = -log10(padjust), fill = sig, label = gene_id)) +
    geom_point(shape = 21, colour = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = 1, color = "black") +
    geom_vline(xintercept = -1, color = "black") +
    scale_fill_manual(values = c("yes" = "#8B0000", "no" = "grey")) +
    geom_text(data = subset(df_for_plotting_immune, sig == "yes"), hjust = -0.1, size = 3) +
    facet_grid(sample~category) +
    scale_y_continuous(limits = c(0,11)) +
    labs(title = "Differential expression of immune genes") +
    theme_bw (),
  ncol = 2, nrow = 1, common.legend = TRUE, widths = c(2.25,1)
)

# export plot
pdf("R_exports/plots/Fig_S6_volcano_differential_expression_gene_lists.pdf")
plot(plot)
dev.off()



##### calculate numbers to add to plot #####

### information to be added manually to the plot for each sample
# number of significantly upregulated genes for each class
# percentage of all genes in list that are upregulated

# CTAs, TSGs and oncogenes
summary_diff_expression_gene_lists <-
  data.frame("list" = (c(rep("CTA",6), rep("TSG",6), rep("onco",6))))
summary_diff_expression_gene_lists$treatment <- rep(c("LD3", "LG3", "LG7", "SD3", "SG3", "SG7"),3)

count <- c()
percent <- c()

for (c in c("CTA", "TSG", "onco")) {
  for (s in c("LD3", "LG3", "LG7", "SD3", "SG3", "SG7")) {
    count <- c(count,
               nrow(df_for_plotting[which(
                 df_for_plotting$category == c &
                   df_for_plotting$sample == s &
                   df_for_plotting$sig == "yes" &
                   df_for_plotting$logFC > 1),]))
    percent <- c(percent,
                 nrow(df_for_plotting[which(
                   df_for_plotting$category == c &
                     df_for_plotting$sample == s &
                     df_for_plotting$sig == "yes" &
                     df_for_plotting$logFC > 1),])/
                   nrow(df_for_plotting[which(
                     df_for_plotting$category == c &
                       df_for_plotting$sample == s),])*100)
  }
}

summary_diff_expression_gene_lists$count <- count
summary_diff_expression_gene_lists$percentage <- percent

# immune genes
df_temp <- data.frame("list" = rep("immune",6))
df_temp$treatment <- c("LD3", "LG3", "LG7", "SD3", "SG3", "SG7")

count <- c()
percent <- c()

for (s in c("LD3", "LG3", "LG7", "SD3", "SG3", "SG7")) {
  count <- c(count,
             nrow(df_for_plotting_immune[which(
               df_for_plotting_immune$sample == s &
                 df_for_plotting_immune$sig == "yes" &
                 df_for_plotting_immune$logFC > 1),]))
  percent <- c(percent,
               nrow(df_for_plotting_immune[which(
                 df_for_plotting_immune$sample == s &
                   df_for_plotting_immune$sig == "yes" &
                   df_for_plotting_immune$logFC > 1),])/
                 nrow(df_for_plotting_immune[which(
                   df_for_plotting$sample == s),])*100)
  
}

df_temp$count <- count
df_temp$percentage <- percent

summary_diff_expression_gene_lists <- rbind(summary_diff_expression_gene_lists, df_temp)

remove(count)
remove(percent)
remove(c)
remove(s)
remove(plot)
remove(df_temp)
remove(df_for_plotting)
remove(df_for_plotting_immune)