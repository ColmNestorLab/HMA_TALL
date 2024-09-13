################################################################################
######################## CELL CYCLE PLOT AND STATISTICS ########################
################################################################################

# Maike Bensberg (maike.bensberg@liu.se)
# 2024-07-10


# plot stacked barplots for cell cycle distributions
# statistics comparing different treatments

# Figures: Fig. 2e



##### load packages #####

library(ggsci)
library(ggplot2)
library(scales)
library(reshape2)



##### set working directory and check package versions #####

# set working directory
setwd("P:/HMA project/R_analysis/analysis_manuscript")
# print working directory
getwd()

# print version of R and attached packages
sessionInfo()



##### define color palettes #####

nejm_colors <- pal_nejm()(8)
show_col (nejm_colors)



#### import data #####


### cell cycle analysis results
# mean from biological triplicates
# column names: "Cell_line","drug","concentration","cell_cycle","percent" 
cell_cycle_treated_JS <- read.csv("data/JS_treated_cell_cycle_average_triplicates.csv")



##### plot cell cycle distribution: stacked barplot #####


### Figure 2e
# stacked bar plot
# mean of biological triplicates


### colors based on pre-defined color-palettes
# G1: nejm [2] (blue); G2/M: nejm [1] (red); S: nejm [4] (green)

plot_2e <-
  ggplot(data = cell_cycle_treated_JS, aes (x = interaction(concentration, drug), y = percent, fill = cell_cycle)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~Cell_line) +
  geom_text(aes (label = round(percent, digits = 0)), position = position_fill(vjust = 0.5), color = "white") +
  scale_x_discrete(limits = c("0.control", "300.AZA", "300.DEC", "300.GSK", "3000.GSK")) +
  scale_fill_manual(values = c("G1" = nejm_colors [2], "G2" = nejm_colors [1], "S" = nejm_colors [4])) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


### export plot

pdf("plots_for_publication/2e_cell_cycle_stacked_bar_JS_triplicates.pdf")
plot(plot_2e)
dev.off()



##### statistics #####


### comparison of cell cycle distribution for different treatments

# fisher's exact test
# compute p-values by Monte Carlo simulation (simulate.p.value)


### create data matrix
# for each cell line

# JURKAT
cell_cycle_matrix_J <- acast(data = cell_cycle_treated_JS[which(cell_cycle_treated_JS$Cell_line == "JURKAT"),],
                             cell_cycle~drug+concentration, value.var = "percent")

# SUP-T1
cell_cycle_matrix_S <- acast(data = cell_cycle_treated_JS[which(cell_cycle_treated_JS$Cell_line == "SUPT1"),],
                             cell_cycle~drug+concentration, value.var = "percent")


### compare each treatment to control

# make dataframe summarizing results for all treatments
fisher_by_treatment <- data.frame(cell_line = c(rep("JURKAT", ncol(cell_cycle_matrix_J)-1),rep("SUPT1", ncol(cell_cycle_matrix_S)-1)))
fisher_by_treatment$treatment <- c(colnames(cell_cycle_matrix_J)[!colnames(cell_cycle_matrix_J) == "control_0"],
                                   colnames(cell_cycle_matrix_S)[!colnames(cell_cycle_matrix_S) == "control_0"])

# calculate p-values
temp <- c()
for (t in colnames(cell_cycle_matrix_J)[!colnames(cell_cycle_matrix_J) == "control_0"]) {
  temp <- c(temp,
            fisher.test(cell_cycle_matrix_J[,c(t,"control_0")], simulate.p.value= TRUE)$p.value)
}
for (t in colnames(cell_cycle_matrix_S)[!colnames(cell_cycle_matrix_S) == "control_0"]) {
  temp <- c(temp,
            fisher.test(cell_cycle_matrix_S[,c(t,"control_0")], simulate.p.value= TRUE)$p.value)
}
fisher_by_treatment$p_value <- temp

remove(temp)
remove(t)

# adjust p-values for multiple testing
fisher_by_treatment$p_adjusted <- c(p.adjust(p = fisher_by_treatment[which(fisher_by_treatment$cell_line == "JURKAT"),"p_value"], method = "BH"),
                                    p.adjust(p = fisher_by_treatment[which(fisher_by_treatment$cell_line == "SUPT1"),"p_value"], method = "BH"))