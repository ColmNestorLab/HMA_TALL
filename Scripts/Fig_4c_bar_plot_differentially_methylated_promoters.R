################################################################################
################ BAR PLOT: differentially methylated promoters ################# 
################################################################################

# Maike Bensberg (maike.bensberg@liu.se)
# 2024-08-13


### make a bar plot showing the number of promoters that lose or gain DNA methylation
# based on analysis of differentially methylated promoters
# script for pre-processing the data: promoter_methylation_analysis.R

### definition for loss/gain of DNA methylation
# methylation difference: >= 25%
# adjusted p-value (BH-adjusted) < 0.01


### Figure: Fig. 4c



##### initiate R environment #####
# only needs to be done once
#install.packages("renv")
#renv::init()



##### loading packages #####

library(ggsci)
library(scales)
library(ggplot2)



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



##### import and format data #####


### import

# Supplementary Table 8
# processed and exported in script: promoter_methylation_analysis.R
df_DMR_promoters <- read.delim("R_exports/tables/Supplementary_Table_8.csv", sep = ",")

### make new dataframe for plotting

# combine data for LOUCY and SUP-T1
df_DMR_promoters_plotting <- df_DMR_promoters

# add column for methylation change
# gain: methylation difference >= 25 and adjusted p < 0.01
# loss: methylation difference <= -25 and adjusted p < 0.01
# no change: all others

df_DMR_promoters_plotting$methylation <- "no change"
df_DMR_promoters_plotting[which(df_DMR_promoters_plotting$qvalue < 0.01 &
                                  df_DMR_promoters_plotting$meth.diff >= 25),]$methylation <- "gain"
df_DMR_promoters_plotting[which(df_DMR_promoters_plotting$qvalue < 0.01 &
                                  df_DMR_promoters_plotting$meth.diff <= -25),]$methylation <- "loss"



##### plot #####

# bar plot showing the number of differentially methylated promoters for each sample
# colours based on pre-defined colour palletts
# differential methylation divided by color (DNA methylation loss: blue, gain: red; none: grey)

plot <- 
  ggplot(data = df_DMR_promoters_plotting, aes(x = sample, fill = factor(methylation, levels = c("no change", "gain", "loss")))) +
  geom_bar(position = "fill", ) +
  scale_fill_manual(values = c("no change" = "azure3", "loss" = npg_colors[4], "gain" = npg_colors[1]),
                    name = "methylation") +
  labs(title = "promoters methylation changes") +
  geom_text(aes (label = after_stat(count)),  stat = "count", position = position_fill(vjust = 0.5)) +
  theme_bw()

pdf("R_exports/plots/Fig_4c_DMR_promoters_stacked_bar.pdf")
plot(plot)
dev.off()

remove(plot)
remove(df_DMR_promoters_plotting)
remove(npg_colors)