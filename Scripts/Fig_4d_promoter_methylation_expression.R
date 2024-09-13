################################################################################
################### PROMOTER METHYLATION AND EXPRESSION PLOT ################### 
################################################################################

# Maike Bensberg (maike.bensberg@liu.se)
# 2024-08-22


### plot DNA methylation and differential expression
# promoter methylation analyzed for every covered curated transcript (NM/NR)
# scripts for promoter methylation analysis: promoter_methylation_analysis.R
# expression analyzed by Sandra H. on gene level
# scatterplot shows DNA methylation and expression for every covered promoter
# might includeseveral transcripts per gene


### Figure: Fig. 4d



##### initiate R environment #####
# only needs to be done once
#install.packages("renv")
#renv::init()



##### loading packages #####

library(ggsci)
library(scales)
library(ggplot2)



##### defining color palettes #####

npg_colors <- pal_npg()(10)
show_col (npg_colors)



##### set working directory and check package versions #####

# set working directory
setwd("/mnt/wwn-0x5000cca28fd30c7d-part1/HMA_TALL/EMseq/LOUCY_SUPT1_HMAtreated_12_2023/methylKit")
# print working directory
getwd()

# print version of R and attached packages
sessionInfo()



##### import and format data #####


### import data

# dataframe with information on DNA methylation and differential expression
# data combined in script combining_DNAmethylation_and_expression_data.R
# only curated transcripts included (NM/NR)
# Supplementary Table 11

methylation_expression_N_transcripts <- read.delim("R_exports/tables/Supplementary_Table_11.csv", sep = ",")


### make dataframe for plotting
# extract only relevant information


### add a column for plotting color
# grey: methylation difference < 25% and/or absolute expression log2(FC) <= 1
# black: methylation difference >= 25% and/or absolute expression log2(FC) > 1
# red: like black but also significant differential promoter methylation (adjusted p < 0.01) 
#       and significant differential expression (adjusted p < 0.05)

df_for_plotting <-
  methylation_expression_N_transcripts[,c("chr", "start", "end", "transcript_id",
                                          "qvalue", "meth.diff", "logFC", "padjust", "sample")]

df_for_plotting$color <- "black"
df_for_plotting[which(df_for_plotting$qvalue < 0.01 &
                        df_for_plotting$padjust < 0.05),]$color <- "red"
df_for_plotting[which((df_for_plotting$meth.diff >= -25 &
                         df_for_plotting$meth.diff <= 25) |
                        (df_for_plotting$logFC > -1 &
                           df_for_plotting$logFC < 1)),]$color <- "grey"



##### plot #####

plot <- 
  ggplot(data = df_for_plotting, aes(x = logFC, y = meth.diff, color = color)) +
  geom_point(stroke = NA, size = 1) +
  scale_color_manual(values = c("grey" = "azure3", "black" = "black", "red" = npg_colors[8])) +
  geom_hline(yintercept = c(-25,25), linetype = "dashed") +
  geom_vline(xintercept = c(-1,1), linetype = "dashed") +
  facet_wrap(~sample, ncol = 2) +
  theme_bw()

# export plot
pdf("R_exports/plots/Fig_4d_promoter_methylation_gene_expression.pdf")
plot(plot)
dev.off()

remove(plot)



##### calculate numbers to add to plot #####

### number of transcripts that
# are significantly upregulated (adjusted p < 0.05) and have a logFC of > 1
# are differentially methylated (adjusted p < 0.01) and lose > 25% DNA methylation

### count total number of transcripts and calculate percentage

df_methylation_loss_upregulated_transcripts <-
  data.frame("sample" = c("LD3", "LG3", "LG7", "SD3", "SG3", "SG7"))

count <- c()
total <- c()
for (s in c("LD3", "LG3", "LG7", "SD3", "SG3", "SG7")) {
  count <- c(count,
             nrow(df_for_plotting[which(df_for_plotting$sample == s &
                                          df_for_plotting$qvalue < 0.01 &
                                          df_for_plotting$meth.diff <= -25 &
                                          df_for_plotting$padjust < 0.05 &
                                          df_for_plotting$logFC > 1),]))
  total <- c(total,
             nrow(df_for_plotting[which(df_for_plotting$sample == s),]))
}

df_methylation_loss_upregulated_transcripts$total <- count
df_methylation_loss_upregulated_transcripts$percent <- round(count/total*100, digits = 2)


remove(df_for_plotting)