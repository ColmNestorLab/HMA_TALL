################################################################################
################################## QPCR PLOTS ##################################
################################################################################

# Maike Bensberg (maike.bensberg@liu.se)
# 2024-07-11


# plot expression based on qPCR results
# plot barplots of relative expression
# 2^(delta Ct) values (normalized to GAPDH expression)

# Figures: Extended Data Fig. 1c, Extended Data Fig. 3e 



##### load packages #####

library(ggpubr)
library(ggsci)
library(scales)



##### set working directory and check package versions #####

# set working directory
setwd("P:/HMA project/R_analysis/analysis_manuscript")
# print working directory
getwd()

# print version of R and attached packages
sessionInfo()



##### define color palettes #####

npg_colors <- pal_npg()(10)
show_col (npg_colors)
nejm_colors <- pal_nejm()(8)
show_col (nejm_colors)



#### import data #####

### qPCR results
# column names: "cell_line", "drug","concentration","treatment_duration","gene","average_Ct","deltaCt","X2deltaCt","X2Ct"
# not all columns are needed (average_Ct, deltaCt and X2Ct are not used)
# empty entries when no expression was detected (no Ct value; never crossing threshold)
qPCR_treated_JLSA <- read.csv("data/qPCRs_treated_JLSA_treated.csv")



##### plot: barplots #####

# plot 2^(delta Ct) values (normalized to GAPDH expression)


### colors based on pre-defined color-palettes
# AZA: nejm [5] (purple); DEC: nejm [7] (yellow); GSK: nejm [4] (green); control: "azure3" (grey)


### Extended Data Fig. 1c

# cell lines: JURKAT, LOUCY, SUP-T1
CL <- c("JURKAT", "LOUCY", "SUPT1")

# genes: DNMT1, DNMT3A, DNMT3B
genes <- c("DNMT1", "DNMT3A", "DNMT3B")

plot_S1c <- ggarrange(
  # 3 days treatment
  ggplot(data = qPCR_treated_JLSA[which(qPCR_treated_JLSA$treatment_duration == "3days" &
                                          qPCR_treated_JLSA$gene %in% genes &
                                          qPCR_treated_JLSA$cell_line %in% CL),],
         aes(x = interaction(concentration, drug), y = X2deltaCt, fill = drug)) +
    geom_bar(stat = "identity", color = "black") +
    facet_wrap(gene~cell_line, scales = "free_y", ncol = 3) +
    scale_x_discrete(limits = c("0.control", "10.AZA", "30.AZA", "100.AZA", "300.AZA", "1000.AZA", "3000.AZA",
                                "10.DEC", "30.DEC", "100.DEC", "300.DEC", "1000.DEC", "3000.DEC",
                                "10.GSK", "30.GSK", "100.GSK", "300.GSK", "1000.GSK", "3000.GSK")) +
    scale_fill_manual(values = c("control" = "azure3", "DEC" = nejm_colors [7], "AZA" = nejm_colors [5], "GSK" = npg_colors [3])) +
    labs(title = "3 days treatment") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)),
  # 7 days treatment
  ggplot(data = qPCR_treated_JLSA[which(qPCR_treated_JLSA$treatment_duration == "7days" &
                                          qPCR_treated_JLSA$gene %in% genes &
                                          qPCR_treated_JLSA$cell_line %in% CL),],
         aes(x = interaction(concentration, drug), y = X2deltaCt, fill = drug)) +
    geom_bar(stat = "identity", color = "black") +
    facet_wrap(gene~cell_line, scales = "free_y", ncol = 3) +
    scale_x_discrete(limits = c("0.control", "10.AZA", "30.AZA","10.DEC", "30.DEC",
                                "10.GSK", "30.GSK", "100.GSK", "300.GSK", "1000.GSK", "3000.GSK")) +
    scale_fill_manual(values = c("control" = "azure3", "DEC" = nejm_colors [7], "AZA" = nejm_colors [5], "GSK" = npg_colors [3])) +
    labs(title = "7 days treatment") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)),
  common.legend = TRUE, ncol = 1, nrow = 2
)

# export plot
pdf("plots_for_publication/S1c_qPCR_DNMTs_JLS.pdf")
plot(plot_S1c)
dev.off()

remove(CL)
remove(genes)


### Extended Data Fig. 3e

# cell lines: ALL-SIL JURKAT, LOUCY, SUP-T1
CL <- c("ALLSIL", "JURKAT", "LOUCY", "SUPT1")

# genes: DAZL, GAGE12
genes <- c("DAZL", "GAGE")
  
plot_S3e <- ggarrange(
  # 3 days treatment
  ggplot(data = qPCR_treated_JLSA[which(qPCR_treated_JLSA$treatment_duration == "3days" &
                                          qPCR_treated_JLSA$gene %in% genes &
                                          qPCR_treated_JLSA$cell_line %in% CL),],
         aes(x = interaction(concentration, drug), y = X2deltaCt, fill = drug)) +
    geom_bar(stat = "identity", color = "black") +
    facet_wrap(cell_line~gene, scales = "free_y", ncol = 2) +
    scale_x_discrete(limits = c("0.control", "10.AZA", "30.AZA", "100.AZA", "300.AZA", "1000.AZA", "3000.AZA",
                                "10.DEC", "30.DEC", "100.DEC", "300.DEC", "1000.DEC", "3000.DEC",
                                "10.GSK", "30.GSK", "100.GSK", "300.GSK", "1000.GSK", "3000.GSK")) +
    scale_fill_manual(values = c("control" = "azure3", "DEC" = nejm_colors [7], "AZA" = nejm_colors [5], "GSK" = npg_colors [3])) +
    labs(title = "3 days treatment") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)),
  # 7 days treatment
  ggplot(data = qPCR_treated_JLSA[which(qPCR_treated_JLSA$treatment_duration == "7days" &
                                          qPCR_treated_JLSA$gene %in% genes &
                                          qPCR_treated_JLSA$cell_line %in% CL),],
         aes(x = interaction(concentration, drug), y = X2deltaCt, fill = drug)) +
    geom_bar(stat = "identity", color = "black") +
    facet_wrap(cell_line~gene, scales = "free_y", ncol = 2) +
    scale_x_discrete(limits = c("0.control", "10.AZA", "30.AZA","10.DEC", "30.DEC",
                                "10.GSK", "30.GSK", "100.GSK", "300.GSK", "1000.GSK", "3000.GSK")) +
    scale_fill_manual(values = c("control" = "azure3", "DEC" = nejm_colors [7], "AZA" = nejm_colors [5], "GSK" = npg_colors [3])) +
    labs(title = "7 days treatment") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)),
  common.legend = TRUE, ncol = 2, nrow = 1
)

# export plot
pdf("plots_for_publication/S3e_qPCR_DAZL_GAGE12_AJLS.pdf")
plot(plot_S3e)
dev.off()

remove(CL)
remove(genes)