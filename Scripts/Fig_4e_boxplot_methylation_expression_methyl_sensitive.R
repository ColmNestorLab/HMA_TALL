################################################################################
###### BOXPLOT: methylation and expression of methylation sensitive genes ###### 
################################################################################

# Maike Bensberg (maike.bensberg@liu.se)
# 2024-08-22


### plot DNA methylation and expression of methylation sensitive transcripts
# genes defined as methylation sensitive after treatment with GSK for 7 days


### criteria for methylation sensitive genes/transcripts
# not expressed in untreated cells (cpm < 0.5)
# expressed in treated cells (cpm >= 0.5)
# significantly upreguated (adjusted p < 0.05 and log2FC > 1)
# promoters had lost > 25% of DNA methylation
# promoters were differentially methylated (adjusted p < 0.01)


### DNA methylation data and expression data combined and methylation-sensitive genes defined
# see script combining_DNAmethylation_and_expression_data.R

### Figure: Fig. 4e



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



##### defining color palettes #####

npg_colors <- pal_npg()(10)
show_col (npg_colors)
nejm_colors <- pal_nejm()(8)
show_col (nejm_colors)



##### import and format data #####

# dataframe with information on DNA methylation and differential expression
# data combined in script combining_DNAmethylation_and_expression_data.R
# only curated transcripts included (NM/NR)
# Supplementary Table 13

methylation_expression_N_transcripts <- read.delim("R_exports/tables/Supplementary_Table_13.csv", sep = ",")


### make dataframe for plotting
# extract only relevant information

data_temp <- methylation_expression_N_transcripts

# add promoter ID as new column
data_temp$id <- paste(data_temp$chr, data_temp$start, data_temp$end, sep = ":")


### which promoters are methylation sensitive (after treatment with GSK for 7 days)?
LOUCY_methyl_sensitive_promoters <-
  data_temp[which(data_temp$methyl_sensitive == "yes" &
                    data_temp$sample == "LG7"),"id"]

SUPT1_methyl_sensitive_promoters <-
  data_temp[which(data_temp$methyl_sensitive == "yes" &
                    data_temp$sample == "SG7"),"id"]


### extract expression (log2cpm) and methylation for methylation sensitive promoters

df_for_plotting <- 
  data.frame("sample" = c(rep("Lctrl7", length(LOUCY_methyl_sensitive_promoters)),
                          rep("LD3", length(LOUCY_methyl_sensitive_promoters)),
                          rep("LG3", length(LOUCY_methyl_sensitive_promoters)),
                          rep("LG7", length(LOUCY_methyl_sensitive_promoters)),
                          rep("Sctrl7", length(SUPT1_methyl_sensitive_promoters)),
                          rep("SD3", length(SUPT1_methyl_sensitive_promoters)),
                          rep("SG3", length(SUPT1_methyl_sensitive_promoters)),
                          rep("SG7", length(SUPT1_methyl_sensitive_promoters))))

# gene_id
df_for_plotting$gene_id <-
  c(rep(data_temp[which(data_temp$id %in% LOUCY_methyl_sensitive_promoters &
                          data_temp$sample == "LG7"),"gene_id"],4),
    rep(data_temp[which(data_temp$id %in% SUPT1_methyl_sensitive_promoters &
                          data_temp$sample == "SG7"),"gene_id"],4))

# expression (log2cpm)
df_for_plotting$expression <-
  c(data_temp[which(data_temp$id %in% LOUCY_methyl_sensitive_promoters &
                      data_temp$sample == "LG7"), "log2cpm_control"],
    data_temp[which(data_temp$id %in% LOUCY_methyl_sensitive_promoters &
                      data_temp$sample == "LD3"), "log2cpm_treated"],
    data_temp[which(data_temp$id %in% LOUCY_methyl_sensitive_promoters &
                      data_temp$sample == "LG3"), "log2cpm_treated"],
    data_temp[which(data_temp$id %in% LOUCY_methyl_sensitive_promoters &
                      data_temp$sample == "LG7"), "log2cpm_treated"],
    data_temp[which(data_temp$id %in% SUPT1_methyl_sensitive_promoters &
                      data_temp$sample == "SG7"), "log2cpm_control"],
    data_temp[which(data_temp$id %in% SUPT1_methyl_sensitive_promoters &
                      data_temp$sample == "SD3"), "log2cpm_treated"],
    data_temp[which(data_temp$id %in% SUPT1_methyl_sensitive_promoters &
                      data_temp$sample == "SG3"), "log2cpm_treated"],
    data_temp[which(data_temp$id %in% SUPT1_methyl_sensitive_promoters &
                      data_temp$sample == "SG7"), "log2cpm_treated"])

# methylation across promoter (values from 0 to 1
df_for_plotting$methylation <-
  c(data_temp[which(data_temp$id %in% LOUCY_methyl_sensitive_promoters &
                      data_temp$sample == "LG7"), "methylation_control"],
    data_temp[which(data_temp$id %in% LOUCY_methyl_sensitive_promoters &
                      data_temp$sample == "LD3"), "methylation_treated"],
    data_temp[which(data_temp$id %in% LOUCY_methyl_sensitive_promoters &
                      data_temp$sample == "LG3"), "methylation_treated"],
    data_temp[which(data_temp$id %in% LOUCY_methyl_sensitive_promoters &
                      data_temp$sample == "LG7"), "methylation_treated"],
    data_temp[which(data_temp$id %in% SUPT1_methyl_sensitive_promoters &
                      data_temp$sample == "SG7"), "methylation_control"],
    data_temp[which(data_temp$id %in% SUPT1_methyl_sensitive_promoters &
                      data_temp$sample == "SD3"), "methylation_treated"],
    data_temp[which(data_temp$id %in% SUPT1_methyl_sensitive_promoters &
                      data_temp$sample == "SG3"), "methylation_treated"],
    data_temp[which(data_temp$id %in% SUPT1_methyl_sensitive_promoters &
                      data_temp$sample == "SG7"), "methylation_treated"])



##### plot #####


### fill colour defined by treatment
# control: grey; 10 nM DEC 3d (LD3): nejm_colors [7] (yellow);
# 300 nM GSK 3 days (LG3): npg_colors [10] (brown); 300 nM GSK 7 days (LG7): npg_colors [3] (green)


### statistics added to plot
# two-sample Wilcoxon test comparing treatment to control

plot <- ggarrange(
  # LOUCY methylation
  ggplot(data = df_for_plotting[which(df_for_plotting$sample %in% c("Lctrl7", "LD3", "LG3", "LG7")),],
         aes(x = sample, y = methylation, color = sample)) +
    geom_point(stroke = NA, size = 1, position = position_jitterdodge(jitter.width = 0.7)) +
    geom_boxplot(outlier.shape = NA, color = "black", fill = NA) +
    scale_color_manual(labels = c("Lctrl7" = "control", "LD3" = "10 nM DEC 3d", "LG3" = "300 nM GSK 3 days", "LG7" = "300 nM GSK 7 days"),
                       values = c("Lctrl7" = "azure3", "LD3" = nejm_colors [7], "LG3" = npg_colors [10], "LG7" = npg_colors [3])) +
    stat_compare_means(method = "wilcox.test", label = "p.signif", ref.group = "Lctrl7") +
    labs(title = "LOUCY methylation", y = "% CpG methylation") +
    theme_bw(),
  # SUP-T1 methylation
  ggplot(data = df_for_plotting[which(df_for_plotting$sample %in% c("Sctrl7", "SD3", "SG3", "SG7")),],
         aes(x = sample, y = methylation, color = sample)) +
    geom_point(stroke = NA, size = 1, position = position_jitterdodge(jitter.width = 0.7)) +
    geom_boxplot(outlier.shape = NA, color = "black", fill = NA) +
    scale_color_manual(labels = c("Sctrl7" = "control", "SD3" = "10 nM DEC 3d", "SG3" = "300 nM GSK 3 days", "SG7" = "300 nM GSK 7 days"),
                       values = c("Sctrl7" = "azure3", "SD3" = nejm_colors [7], "SG3" = npg_colors [10], "SG7" = npg_colors [3])) +
    stat_compare_means(method = "wilcox.test", label = "p.signif", ref.group = "Sctrl7") +
    labs(title = "SUP-T1 methylation", y = "% CpG methylation") +
    theme_bw(),
  # LOUCY expression
  ggplot(data = df_for_plotting[which(df_for_plotting$sample %in% c("Lctrl7", "LD3", "LG3", "LG7")),],
         aes(x = sample, y = expression, color = sample)) +
    geom_point(stroke = NA, size = 1, position = position_jitterdodge(jitter.width = 0.7)) +
    geom_boxplot(outlier.shape = NA, color = "black", fill = NA) +
    scale_color_manual(labels = c("Lctrl7" = "control", "LD3" = "10 nM DEC 3d", "LG3" = "300 nM GSK 3 days", "LG7" = "300 nM GSK 7 days"),
                       values = c("Lctrl7" = "azure3", "LD3" = nejm_colors [7], "LG3" = npg_colors [10], "LG7" = npg_colors [3])) +
    geom_hline(yintercept = log2(0.5), linetype = "dashed") +
    stat_compare_means(method = "wilcox.test", label = "p.signif", ref.group = "Lctrl7") +
    scale_y_continuous(limits = c(-7,7)) +
    labs(title = "LOUCY expression", y = "expression (cpm(log2))") +
    theme_bw(),
  # SUP-T1 expression
  ggplot(data = df_for_plotting[which(df_for_plotting$sample %in% c("Sctrl7", "SD3", "SG3", "SG7")),],
         aes(x = sample, y = expression, color = sample)) +
    geom_point(stroke = NA, size = 1, position = position_jitterdodge(jitter.width = 0.7)) +
    geom_boxplot(outlier.shape = NA, color = "black", fill = NA) +
    scale_color_manual(labels = c("Sctrl7" = "control", "SD3" = "10 nM DEC 3d", "SG3" = "300 nM GSK 3 days", "SG7" = "300 nM GSK 7 days"),
                       values = c("Sctrl7" = "azure3", "SD3" = nejm_colors [7], "SG3" = npg_colors [10], "SG7" = npg_colors [3])) +
    geom_hline(yintercept = log2(0.5), linetype = "dashed") +
    stat_compare_means(method = "wilcox.test", label = "p.signif", ref.group = "Sctrl7") +
    scale_y_continuous(limits = c(-7,7)) +
    labs(title = "SUP-T1 expression", y = "expression (cpm(log2))") +
    theme_bw(),
  nrow = 2, ncol = 2, common.legend = TRUE, legend = "bottom", align = "hv"
)

# export plot
pdf("R_exports/plots/Fig_4e_boxplot_methylation_expression_methyl_sensitive.pdf")
plot(plot)
dev.off()


remove(data_temp)
remove(LOUCY_methyl_sensitive_promoters)
remove(SUPT1_methyl_sensitive_promoters)
remove(df_for_plotting)