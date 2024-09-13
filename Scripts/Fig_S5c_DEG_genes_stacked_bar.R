################################################################################
################### BAR PLOT: differentially expressed genes ################### 
################################################################################

# Maike Bensberg (maike.bensberg@liu.se)
# 2024-08-28


### plot a bar plot showing the number of differentially expressed genes
# based on DEG analysis by Sandra H
# summarized in supplementary Table 10


### definition of differential expression
# upregulated: logFC > 1 and adjusted p < 0.05
# downregulated: logFC < -1 and adjusted p < 0.05

### Figure: Extended Data Fig. 5c



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



##### import data #####

# DEG analysis done by Sandra H.
# individual files were combined and exported in script combining_DNAmethylation_and_expression_data.R
# Supplementary Table 10

### expression change was summarized in column "expression"
# upregulated: logFC > 1 and adjusted p < 0.05
# downregulated: logFC < -1 and adjusted p < 0.05
# no change: all others

df_expression_combined <- read.delim("R_exports/tables/Supplementary_Table_10.csv", sep = ",")



##### plot ####

# bar plot showing the number of differentially expressed genes for each sample
# colours based on pre-defined colour palletts
# expression changes divided by color (downregulated: blue, upregulated: red; no change: grey)

plot <- 
  ggplot(data = df_expression_combined, aes(x = sample, fill = factor(expression, levels = c("no change", "upregulated", "downregulated")))) +
  geom_bar(position = "fill", ) +
  scale_fill_manual(values = c("no change" = "azure3", "downregulated" = npg_colors[4], "upregulated" = npg_colors[1]),
                    name = "expression") +
  labs(title = "differential gene expression") +
  geom_text(aes (label = after_stat(count)),  stat = "count", position = position_fill(vjust = 0.5)) +
  theme_bw()

pdf("R_exports/plots/Fig_S5c_DEG_genes_stacked_bar.pdf")
plot(plot)
dev.off()

remove(plot)