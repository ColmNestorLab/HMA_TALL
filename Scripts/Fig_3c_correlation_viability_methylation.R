################################################################################
###################### CORRELATION VIABILITY METHYLATION #######################
################################################################################

# Maike Bensberg (maike.bensberg@liu.se)
# 2024-07-10


# plot cell viability and global CpG methylation for treated cells
# statistics: correlation: spearman' correlation's rho

# Figure: Fig. 3c



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

### cell viability and global DNA methylation for same treatments
# column names: "cell_line","drug","concentration","treatment_duration","viability_mean","CpGmethylation" 
viability_methylation_all <- read.csv("data/viability_and_methylation_AJLS_all.csv")



##### correlation scatter plot #####


### statistics:
# fit linear and logarithmic curve to data (for all cell lines combined)
# add 95% confidence interval
# correlation analysis: spearman's rho for all cell lines combined (independent of fitted line)


### visuals:
# defined colour for each cell line based on pre-defined colour pallettes 
# ALLSIL : nejm [1] (red), JURKAT : nejm [5] (purple), LOUCY : nejm [6] (blue), SUPT1 : npg [5] (peach)


### only data for 3 days treatment

# plots for all three drugs included in Figure 3e (first 3 panels)
# linear regression fit well, so those curves were used

plot_3c1 <- ggarrange(
  # only points: coloured by cell line
  ggplot() +
    geom_point (data = viability_methylation_all[which(viability_methylation_all$treatment_duration == "3days"),], aes(x = CpGmethylation, y = viability_mean, color = cell_line, shape = treatment_duration),size = 2) +
    scale_color_manual (values = c("ALLSIL" = nejm_colors [1], "JURKAT" = nejm_colors [5], "LOUCY" = nejm_colors [6], "SUPT1" = npg_colors [5])) +
    scale_shape_manual(values = c("3days" = 16, "7days" = 15)) +
    facet_wrap(~drug) +
    scale_x_continuous(limits = c(10,90), breaks = c(25,50,75)) +
    scale_y_continuous(limits = c(0,110), breaks = c(seq(0,100,25))) +
    theme_bw (),
  # spearman's rho
  # linear regression model
  ggplot(data = viability_methylation_all[which(viability_methylation_all$treatment_duration == "3days"),], aes(x = CpGmethylation, y = viability_mean)) +
    geom_point (size = 2) +
    scale_x_continuous(limits = c(10,90), breaks = c(25,50,75)) +
    scale_y_continuous(limits = c(0,110), breaks = c(seq(0,100,25))) +
    stat_cor(method = "spearman", digits = 3) +
    geom_smooth(method = "glm", linewidth = 1.5, color = "black") +
    stat_regline_equation(label.x.npc = "middle") +
    facet_wrap(~drug) +
    theme_bw (),
  # spearman's rho
  # logarithmic regression model
  ggplot(data = viability_methylation_all[which(viability_methylation_all$treatment_duration == "3days"),], aes(x = CpGmethylation, y = viability_mean)) +
    geom_point (size = 2) +
    scale_x_continuous(limits = c(10,90), breaks = c(25,50,75)) +
    scale_y_continuous(limits = c(0,110), breaks = c(seq(0,100,25))) +
    stat_cor(method = "spearman", digits = 3) +
    geom_smooth(method = "glm", formula = "y ~ log(x)", linewidth = 1.5, color = "black") +
    stat_regline_equation(label.x.npc = "middle", formula = "y ~ log(x)") +
    facet_wrap(~drug) +
    theme_bw (),
  ncol = 1, common.legend = TRUE
)

# export plot
pdf("plots_for_publication/3c1_corr_viability_methylation_all_3daysonly.pdf")
plot(plot_3c1)
dev.off()


### data for 3 and 7 day treatment

# only plot for GSK used in Figure 3e (too few datapoints for 7 day treatment with AZA and DAC)
# logarithmic regression fit well, so this curve was used

### each drug for 3 and 7 days
plot_3c2 <- ggarrange(
  # coloured by cell line; shape defined by treatment duration: circle (3 days) and squares (7 days)
  ggplot() +
    geom_point (data = viability_methylation_all, aes(x = CpGmethylation, y = viability_mean, color = cell_line, shape = treatment_duration),size = 2) +
    scale_color_manual (values = c("ALLSIL" = nejm_colors [1], "JURKAT" = nejm_colors [5], "LOUCY" = nejm_colors [6], "SUPT1" = npg_colors [5])) +
    scale_shape_manual(values = c("3days" = 16, "7days" = 15)) +
    facet_wrap(~drug) +
    scale_x_continuous(limits = c(10,90), breaks = c(25,50,75)) +
    scale_y_continuous(limits = c(0,110), breaks = c(seq(0,100,25))) +
    theme_bw (),
  # spearman's rho
  # linear regression model
  ggplot(data = viability_methylation_all, aes(x = CpGmethylation, y = viability_mean)) +
    geom_point (size = 2) +
    scale_x_continuous(limits = c(10,90), breaks = c(25,50,75)) +
    scale_y_continuous(limits = c(0,110), breaks = c(seq(0,100,25))) +
    stat_cor(method = "spearman", digits = 3) +
    geom_smooth(method = "glm", linewidth = 1.5, color = "black") +
    stat_regline_equation(label.x.npc = "middle") +
    facet_wrap(~drug) +
    theme_bw (),
  # spearman's rho
  # logarithmic regression model
  ggplot(data = viability_methylation_all, aes(x = CpGmethylation, y = viability_mean)) +
    geom_point (size = 2) +
    scale_x_continuous(limits = c(10,90), breaks = c(25,50,75)) +
    scale_y_continuous(limits = c(0,110), breaks = c(seq(0,100,25))) +
    stat_cor(method = "spearman", digits = 3) +
    geom_smooth(method = "glm", formula = "y ~ log(x)", linewidth = 1.5, color = "black") +
    stat_regline_equation(label.x.npc = "middle", formula = "y ~ log(x)") +
    facet_wrap(~drug) +
    theme_bw (),
  ncol = 1, common.legend = TRUE
)

# export plot
pdf("plots_for_publication/3c2_corr_viability_methylation_all_3and7days.pdf")
plot(plot_3c2)
dev.off()