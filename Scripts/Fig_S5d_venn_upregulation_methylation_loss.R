################################################################################
###################### VENN: overlap between DEGs and DMR ###################### 
################################################################################

# Maike Bensberg (maike.bensberg@liu.se)
# 2024-08-28


### plot a Venn diagram with overlap between of loss of DNA methylation and expression changes
# differentially methylated promoters (losing >= 25% DNA methylation, p < 0.01)
# differentially expressed genes (absolute FC > 2, p < 0.05)

### DNA methylation and expression data combined in one file
# script: combining_DNAmethylation_and_expression_data.R
# Supplementary Table 13


### Figure: Extended Data Fig. 5d



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


##### import and format data #####

### Supplementary Table 13
# data on methylation and expression was combined and exported in script 
#   combining_DNAmethylation_and_expression_data.R

### expression change was summarized in column "expression"
# upregulated: logFC > 1 and adjusted p < 0.05
# downregulated: logFC < -1 and adjusted p < 0.05
# no change: all others

### methylation was summarized in column methylation
# gain: methylation difference > 25 and adjusted p < 0.01
# loss: methylation difference < -25 and adjusted p < 0.01
# no change: all others

methylation_expression_N_transcripts <- read.delim("R_exports/tables/Supplementary_Table_13.csv", sep = ",")


### extract gene lists
# upregulated genes and genes losing DNA methylation
# for each sample

for (t in c("LD3", "LG3", "LG7", "SD3", "SG3", "SG7")) {
  assign(paste0("expr_", t), unique(
         methylation_expression_N_transcripts[
           which(methylation_expression_N_transcripts$sample == t &
                   methylation_expression_N_transcripts$expression == "upregulated"),"gene_id"]))
  assign(paste0("methyl_", t), unique(
    methylation_expression_N_transcripts[
      which(methylation_expression_N_transcripts$sample == t &
              methylation_expression_N_transcripts$methylation == "loss"),"gene_id"]))
}

remove(t)




##### plot #####

# colours from pre-defined colour paletts
# upregulated genes: npg[1] (red), methylation loss: npg[4] (blue)

plot <- gridExtra::grid.arrange(
  plot(euler(list("upregulated" = expr_LD3,
                  "methylation loss"  = methyl_LD3)),
       quantities = TRUE,
       fills = list(fill = c(npg_colors [1], npg_colors [4]), alpha = 0.5),
       main = "LD3"),
  plot(euler(list("upregulated" = expr_LG3,
                  "methylation loss"  = methyl_LG3)),
       quantities = TRUE,
       fills = list(fill = c(npg_colors [1], npg_colors [4]), alpha = 0.5),
       main = "LG3"),
  plot(euler(list("upregulated" = expr_LG7,
                  "methylation loss"  = methyl_LG7)),
       quantities = TRUE,
       fills = list(fill = c(npg_colors [1], npg_colors [4]), alpha = 0.5),
       main = "LG7"),
  plot(euler(list("upregulated" = expr_SD3,
                  "methylation loss"  = methyl_SD3)),
       quantities = TRUE,
       fills = list(fill = c(npg_colors [1], npg_colors [4]), alpha = 0.5),
       main = "SD3"),
  plot(euler(list("upregulated" = expr_SG3,
                  "methylation loss"  = methyl_SG3)),
       quantities = TRUE,
       fills = list(fill = c(npg_colors [1], npg_colors [4]), alpha = 0.5),
       main = "SG3"),
  plot(euler(list("upregulated" = expr_SG7,
                  "methylation loss"  = methyl_SG7)),
       quantities = TRUE,
       fills = list(fill = c(npg_colors [1], npg_colors [4]), alpha = 0.5),
       main = "SG7")
)


### export

pdf("R_exports/plots/Fig_S5d_venn_upregulation_methylation_loss.pdf")
plot(plot)
dev.off()

### remove temporary objects
remove(plot)
remove(list = ls(pattern = "\\bexpr_..."))
remove(list = ls(pattern = "\\bmethyl_..."))