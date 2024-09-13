################################################################################
##################### CELL VIABILITY CURVES AND STATISTICS #####################
################################################################################

# Maike Bensberg (maike.bensberg@liu.se)
# 2024-07-10


# model cell viability curves on data from alamar blue assays
# plot dose-response curves plus mean and standard error of the mean
# statistics comparing different treatments

# Figures: Fig. 2b, Extended data Fig. 1a, Extended Data Fig. 1b


##### load packages #####

library(drc) # drc: Analysis of Dose-Response Curves (https://rdrr.io/cran/drc/)
library(tidyverse)
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

nejm_colors <- pal_nejm()(8)
show_col (nejm_colors)


#### import data #####

### cell viability for each biological replicate
# column names: "cell_line","drug", "concentration","viability","treatment_duration"
viability_all <- read.csv("data/CellViability_HMA_all.csv")

### mean cell viability, standard deviation and standard error of the mean
# column names: "cell_line","drug", "concentration","treatment_duration","viability_mean","stdev", "sem"
viability_mean <- read.csv("data/CellViability_HMA_mean.csv")


##### format data for dose response curves #####

# based on examples in Ritz et al: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0146021


### make regression curve and calculate values for plotting
# three-parameter log-logistic function where the lower limit is equal to 0
# 95% confidence interval


### one curve for all cell lines combined
# individual model for each drug (t) and separated by treatment duration (d)

for (t in c("AZA", "DEC", "GSK")) {
  for (d in c("3days", "7days")) {
    # create regression curves
    assign(paste0("Reg",t,"_",d),
           drm(viability~concentration,
               data = viability_all[which(viability_all$drug == t & viability_all$treatment_duration == d),],
               fct = LL.3()))
    # new dose levels as support for the line
    assign(paste0("nd",t,"_",d), 
           expand.grid(concentration = exp(seq(log(5),log(11000),length = 500))))
    # predictions and confidence intervals
    assign(paste0("pred",t,"_",d),
           predict(get(paste0("Reg",t,"_",d)), newdata = get(paste0("nd",t,"_",d)),
                   interval = "confidence"))
    # new data with predictions
    pred <- get(paste0("pred",t,"_",d))[,1]
    pmin <- get(paste0("pred",t,"_",d))[,2]
    pmax <- get(paste0("pred",t,"_",d))[,3]
    # add predictions to newdata dataframe
    assign(paste0("nd",t,"_",d),
           cbind(get(paste0("nd",t,"_",d)), pred, pmin, pmax))
  }
}

remove(t)
remove(d)
remove(pred)
remove(pmin)
remove(pmax)


### one curve for each cell line (CL), drug (t) and treatment duration (d)

for (CL in c("ALLSIL", "CCRFCEM", "DND41", "HPBALL", "JURKAT", "LOUCY", "MOLT3", "MOLT4", "PEER", "SUPT1", "TALL1")) {
  for (t in c("AZA", "DEC", "GSK")) {
    for (d in c("3days", "7days")) {
      # create regression curves
      assign(paste0("Reg",CL,t,d),
             drm(viability~concentration,data = viability_all[which(viability_all$drug == t & viability_all$cell_line == CL & viability_all$treatment_duration == d),],fct = LL.3()))
      # new dose levels as support for the line
      assign(paste0("nd",CL,t,d), 
             expand.grid(concentration = exp(seq(log(5),log(11000),length = 500))))
      # predictions and confidence intervals
      assign(paste0("pred",CL,t,d),
             predict(get(paste0("Reg",CL,t,d)), newdata = get(paste0("nd",CL,t,d)),
                     interval = "confidence"))
      # new data with predictions
      pred <- get(paste0("pred",CL,t,d))[,1]
      pmin <- get(paste0("pred",CL,t,d))[,2]
      pmax <- get(paste0("pred",CL,t,d))[,3]
      # add predictions to newdata dataframe
      assign(paste0("nd",CL,t,d),
             cbind(get(paste0("nd",CL,t,d)), pred, pmin, pmax))
    }
  }
}

remove(CL)
remove(t)
remove(d)
remove(pred)
remove(pmin)
remove(pmax)



##### plot dose response curves #####

# three-parameter log-logistic function where the lower limit is equal to 0
# 95% confidence interval

# mean viability and SEM plotted at each concentration


### colors based on pre-defined color-palettes
# AZA: nejm [5] (purple); DEC: nejm [7] (yellow); GSK: nejm [4] (green)



##### Figure 2b #####


### all cell lines combined; compare different treatments

plot_2b <- ggarrange(
  
  ### 3 days
  ggplot() +
    # AZA line and 95% CI
    geom_line(data = ndAZA_3days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [5]) +
    geom_ribbon(data = ndAZA_3days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [5]) +
    # DEC line and 95% CI
    geom_line(data = ndDEC_3days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [7]) +
    geom_ribbon(data = ndDEC_3days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [7]) +
    # GSK line and 95% CI
    geom_line(data = ndGSK_3days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [4]) +
    geom_ribbon(data = ndGSK_3days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [4]) +
    # mean viability and SEM at concentration that have been measured
    geom_point(data = viability_mean [which(viability_mean$treatment_duration == "3days" & viability_mean$cell_line == "all" & viability_mean$concentration != 0),],
               aes(x = concentration, y = viability_mean, shape = drug, colour = drug), size = 2.5) + 
    geom_errorbar(data = viability_mean [which(viability_mean$treatment_duration == "3days" & viability_mean$cell_line == "all" & viability_mean$concentration != 0),],
                  aes(x = concentration, y = viability_mean, ymin = viability_mean - sem, ymax = viability_mean + sem, color = drug), width = .05) +
    labs (x = "concentration [nM]", y = "Cell Viability (% of control)") + 
    # plot formatting (colors, titles, etc.)
    ggtitle("Cell Viability of HMAs in all T-ALL cell lines combined: 3 days treatment") +
    scale_color_manual (values = c("AZA" = nejm_colors [5], "DEC" = nejm_colors [7], "GSK" = nejm_colors [4])) +
    scale_shape_manual(values = c("AZA" = 17, "DEC" = 15, "GSK" = 16)) +
    scale_x_log10() +
    scale_y_continuous(limits = c(0,125), breaks = c(seq(0,100,25))) + 
    geom_hline(aes(yintercept = 50), color = "black", linetype = "dashed") + 
    theme_classic(),
  
  ### 7 days
  ggplot() +
    # AZA line and 95% CI
    geom_line(data = ndAZA_7days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [5]) +
    geom_ribbon(data = ndAZA_7days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [5]) +
    # DEC line and 95% CI
    geom_line(data = ndDEC_7days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [7]) +
    geom_ribbon(data = ndDEC_7days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [7]) +
    # GSK line and 95% CI
    geom_line(data = ndGSK_7days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [4]) +
    geom_ribbon(data = ndGSK_7days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [4]) +
    # mean viability and SEM at concentration that have been measured
    geom_point(data = viability_mean [which(viability_mean$treatment_duration == "7days" & viability_mean$cell_line == "all" & viability_mean$concentration != 0),],
               aes(x = concentration, y = viability_mean, shape = drug, colour = drug), size = 2.5) + 
    geom_errorbar(data = viability_mean [which(viability_mean$treatment_duration == "7days" & viability_mean$cell_line == "all" & viability_mean$concentration != 0),],
                  aes(x = concentration, y = viability_mean, ymin = viability_mean - sem, ymax = viability_mean + sem, color = drug), width = .05) +
    labs (x = "concentration [nM]", y = "Cell Viability (% of control)") + 
    # plot formatting (colors, titles, etc.)
    ggtitle("Cell Viability of HMAs in all T-ALL cell lines combined: 7 days treatment") +
    scale_color_manual (values = c("AZA" = nejm_colors [5], "DEC" = nejm_colors [7], "GSK" = nejm_colors [4])) +
    scale_shape_manual(values = c("AZA" = 17, "DEC" = 15, "GSK" = 16)) +
    scale_x_log10() +
    scale_y_continuous(limits = c(0,125), breaks = c(seq(0,100,25))) + 
    geom_hline(aes(yintercept = 50), color = "black", linetype = "dashed") + 
    theme_classic(),
  
  ncol = 2, legend = "bottom"
)


### export plot

pdf("plots_for_publication/2b_viability_HMAs_all_CellLines.pdf")
plot(plot_2b)
dev.off()



##### Extended Data Figure 1a and 1b #####


### individual graph for each cell line


### Extended Data Figure 1a: 3 days treatment

plot_S1a <- ggarrange (
  # ALL-SIL
  ggplot() +
    # AZA line and 95% CI
    geom_line(data = ndALLSILAZA3days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [5]) +
    geom_ribbon(data = ndALLSILAZA3days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [5]) +
    # DEC line and 95% CI
    geom_line(data = ndALLSILDEC3days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [7]) +
    geom_ribbon(data = ndALLSILDEC3days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [7]) +
    # GSK line and 95% CI
    geom_line(data = ndALLSILGSK3days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [4]) +
    geom_ribbon(data = ndALLSILGSK3days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [4]) +
    # mean viability and SEM at concentration that have been measured
    geom_point(data = viability_mean [which(viability_mean$treatment_duration == "3days" & viability_mean$cell_line == "ALLSIL" & viability_mean$concentration != 0),],
               aes(x = concentration, y = viability_mean, shape = drug, colour = drug), size = 2.5) + 
    geom_errorbar(data = viability_mean [which(viability_mean$treatment_duration == "3days" & viability_mean$cell_line == "ALLSIL" & viability_mean$concentration != 0),],
                  aes(x = concentration, y = viability_mean, ymin = viability_mean - sem, ymax = viability_mean + sem, color = drug), width = .05) +
    labs (x = "concentration [nM]", y = "Cell Viability (% of control)") + 
    # plot formatting (colors, titles, etc.)
    scale_color_manual (values = c("AZA" = nejm_colors [5], "DEC" = nejm_colors [7], "GSK" = nejm_colors [4])) +
    scale_shape_manual(values = c("AZA" = 17, "DEC" = 15, "GSK" = 16)) +
    scale_x_log10() +
    scale_y_continuous(limits = c(0,125), breaks = c(seq(0,100,25))) + 
    geom_hline(aes(yintercept = 50), color = "black", linetype = "dashed") + 
    theme_classic(),
  # CCRF-CEM
  ggplot() +
    # AZA line and 95% CI
    geom_line(data = ndCCRFCEMAZA3days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [5]) +
    geom_ribbon(data = ndCCRFCEMAZA3days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [5]) +
    # DEC line and 95% CI
    geom_line(data = ndCCRFCEMDEC3days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [7]) +
    geom_ribbon(data = ndCCRFCEMDEC3days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [7]) +
    # GSK line and 95% CI
    geom_line(data = ndCCRFCEMGSK3days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [4]) +
    geom_ribbon(data = ndCCRFCEMGSK3days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [4]) +
    # mean viability and SEM at concentration that have been measured
    geom_point(data = viability_mean [which(viability_mean$treatment_duration == "3days" & viability_mean$cell_line == "CCRFCEM" & viability_mean$concentration != 0),],
               aes(x = concentration, y = viability_mean, shape = drug, colour = drug), size = 2.5) + 
    geom_errorbar(data = viability_mean [which(viability_mean$treatment_duration == "3days" & viability_mean$cell_line == "CCRFCEM" & viability_mean$concentration != 0),],
                  aes(x = concentration, y = viability_mean, ymin = viability_mean - sem, ymax = viability_mean + sem, color = drug), width = .05) +
    labs (x = "concentration [nM]", y = "Cell Viability (% of control)") + 
    # plot formatting (colors, titles, etc.)
    scale_color_manual (values = c("AZA" = nejm_colors [5], "DEC" = nejm_colors [7], "GSK" = nejm_colors [4])) +
    scale_shape_manual(values = c("AZA" = 17, "DEC" = 15, "GSK" = 16)) +
    scale_x_log10() +
    scale_y_continuous(limits = c(0,125), breaks = c(seq(0,100,25))) + 
    geom_hline(aes(yintercept = 50), color = "black", linetype = "dashed") + 
    theme_classic(),
  # DND-41
  ggplot() +
    # AZA line and 95% CI
    geom_line(data = ndDND41AZA3days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [5]) +
    geom_ribbon(data = ndDND41AZA3days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [5]) +
    # DEC line and 95% CI
    geom_line(data = ndDND41DEC3days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [7]) +
    geom_ribbon(data = ndDND41DEC3days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [7]) +
    # GSK line and 95% CI
    geom_line(data = ndDND41GSK3days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [4]) +
    geom_ribbon(data = ndDND41GSK3days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [4]) +
    # mean viability and SEM at concentration that have been measured
    geom_point(data = viability_mean [which(viability_mean$treatment_duration == "3days" & viability_mean$cell_line == "DND41" & viability_mean$concentration != 0),],
               aes(x = concentration, y = viability_mean, shape = drug, colour = drug), size = 2.5) + 
    geom_errorbar(data = viability_mean [which(viability_mean$treatment_duration == "3days" & viability_mean$cell_line == "DND41" & viability_mean$concentration != 0),],
                  aes(x = concentration, y = viability_mean, ymin = viability_mean - sem, ymax = viability_mean + sem, color = drug), width = .05) +
    labs (x = "concentration [nM]", y = "Cell Viability (% of control)") + 
    # plot formatting (colors, titles, etc.)
    scale_color_manual (values = c("AZA" = nejm_colors [5], "DEC" = nejm_colors [7], "GSK" = nejm_colors [4])) +
    scale_shape_manual(values = c("AZA" = 17, "DEC" = 15, "GSK" = 16)) +
    scale_x_log10() +
    scale_y_continuous(limits = c(0,125), breaks = c(seq(0,100,25))) + 
    geom_hline(aes(yintercept = 50), color = "black", linetype = "dashed") + 
    theme_classic(),
  # HPB-ALL
  ggplot() +
    # AZA line and 95% CI
    geom_line(data = ndHPBALLAZA3days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [5]) +
    geom_ribbon(data = ndHPBALLAZA3days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [5]) +
    # DEC line and 95% CI
    geom_line(data = ndHPBALLDEC3days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [7]) +
    geom_ribbon(data = ndHPBALLDEC3days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [7]) +
    # GSK line and 95% CI
    geom_line(data = ndHPBALLGSK3days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [4]) +
    geom_ribbon(data = ndHPBALLGSK3days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [4]) +
    # mean viability and SEM at concentration that have been measured
    geom_point(data = viability_mean [which(viability_mean$treatment_duration == "3days" & viability_mean$cell_line == "HPBALL" & viability_mean$concentration != 0),],
               aes(x = concentration, y = viability_mean, shape = drug, colour = drug), size = 2.5) + 
    geom_errorbar(data = viability_mean [which(viability_mean$treatment_duration == "3days" & viability_mean$cell_line == "HPBALL" & viability_mean$concentration != 0),],
                  aes(x = concentration, y = viability_mean, ymin = viability_mean - sem, ymax = viability_mean + sem, color = drug), width = .05) +
    labs (x = "concentration [nM]", y = "Cell Viability (% of control)") + 
    # plot formatting (colors, titles, etc.)
    scale_color_manual (values = c("AZA" = nejm_colors [5], "DEC" = nejm_colors [7], "GSK" = nejm_colors [4])) +
    scale_shape_manual(values = c("AZA" = 17, "DEC" = 15, "GSK" = 16)) +
    scale_x_log10() +
    scale_y_continuous(limits = c(0,125), breaks = c(seq(0,100,25))) + 
    geom_hline(aes(yintercept = 50), color = "black", linetype = "dashed") + 
    theme_classic(),
  # JURKAT
  ggplot() +
    # AZA line and 95% CI
    geom_line(data = ndJURKATAZA3days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [5]) +
    geom_ribbon(data = ndJURKATAZA3days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [5]) +
    # DEC line and 95% CI
    geom_line(data = ndJURKATDEC3days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [7]) +
    geom_ribbon(data = ndJURKATDEC3days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [7]) +
    # GSK line and 95% CI
    geom_line(data = ndJURKATGSK3days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [4]) +
    geom_ribbon(data = ndJURKATGSK3days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [4]) +
    # mean viability and SEM at concentration that have been measured
    geom_point(data = viability_mean [which(viability_mean$treatment_duration == "3days" & viability_mean$cell_line == "JURKAT" & viability_mean$concentration != 0),],
               aes(x = concentration, y = viability_mean, shape = drug, colour = drug), size = 2.5) + 
    geom_errorbar(data = viability_mean [which(viability_mean$treatment_duration == "3days" & viability_mean$cell_line == "JURKAT" & viability_mean$concentration != 0),],
                  aes(x = concentration, y = viability_mean, ymin = viability_mean - sem, ymax = viability_mean + sem, color = drug), width = .05) +
    labs (x = "concentration [nM]", y = "Cell Viability (% of control)") + 
    # plot formatting (colors, titles, etc.)
    scale_color_manual (values = c("AZA" = nejm_colors [5], "DEC" = nejm_colors [7], "GSK" = nejm_colors [4])) +
    scale_shape_manual(values = c("AZA" = 17, "DEC" = 15, "GSK" = 16)) +
    scale_x_log10() +
    scale_y_continuous(limits = c(0,125), breaks = c(seq(0,100,25))) + 
    geom_hline(aes(yintercept = 50), color = "black", linetype = "dashed") + 
    theme_classic(),
  # LOUCY
  ggplot() +
    # AZA line and 95% CI
    geom_line(data = ndLOUCYAZA3days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [5]) +
    geom_ribbon(data = ndLOUCYAZA3days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [5]) +
    # DEC line and 95% CI
    geom_line(data = ndLOUCYDEC3days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [7]) +
    geom_ribbon(data = ndLOUCYDEC3days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [7]) +
    # GSK line and 95% CI
    geom_line(data = ndLOUCYGSK3days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [4]) +
    geom_ribbon(data = ndLOUCYGSK3days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [4]) +
    # mean viability and SEM at concentration that have been measured
    geom_point(data = viability_mean [which(viability_mean$treatment_duration == "3days" & viability_mean$cell_line == "LOUCY" & viability_mean$concentration != 0),],
               aes(x = concentration, y = viability_mean, shape = drug, colour = drug), size = 2.5) + 
    geom_errorbar(data = viability_mean [which(viability_mean$treatment_duration == "3days" & viability_mean$cell_line == "LOUCY" & viability_mean$concentration != 0),],
                  aes(x = concentration, y = viability_mean, ymin = viability_mean - sem, ymax = viability_mean + sem, color = drug), width = .05) +
    labs (x = "concentration [nM]", y = "Cell Viability (% of control)") + 
    # plot formatting (colors, titles, etc.)
    scale_color_manual (values = c("AZA" = nejm_colors [5], "DEC" = nejm_colors [7], "GSK" = nejm_colors [4])) +
    scale_shape_manual(values = c("AZA" = 17, "DEC" = 15, "GSK" = 16)) +
    scale_x_log10() +
    scale_y_continuous(limits = c(0,125), breaks = c(seq(0,100,25))) + 
    geom_hline(aes(yintercept = 50), color = "black", linetype = "dashed") + 
    theme_classic(),
  # MOLT-3
  ggplot() +
    # AZA line and 95% CI
    geom_line(data = ndMOLT3AZA3days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [5]) +
    geom_ribbon(data = ndMOLT3AZA3days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [5]) +
    # DEC line and 95% CI
    geom_line(data = ndMOLT3DEC3days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [7]) +
    geom_ribbon(data = ndMOLT3DEC3days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [7]) +
    # GSK line and 95% CI
    geom_line(data = ndMOLT3GSK3days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [4]) +
    geom_ribbon(data = ndMOLT3GSK3days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [4]) +
    # mean viability and SEM at concentration that have been measured
    geom_point(data = viability_mean [which(viability_mean$treatment_duration == "3days" & viability_mean$cell_line == "MOLT3" & viability_mean$concentration != 0),],
               aes(x = concentration, y = viability_mean, shape = drug, colour = drug), size = 2.5) + 
    geom_errorbar(data = viability_mean [which(viability_mean$treatment_duration == "3days" & viability_mean$cell_line == "MOLT3" & viability_mean$concentration != 0),],
                  aes(x = concentration, y = viability_mean, ymin = viability_mean - sem, ymax = viability_mean + sem, color = drug), width = .05) +
    labs (x = "concentration [nM]", y = "Cell Viability (% of control)") + 
    # plot formatting (colors, titles, etc.)
    scale_color_manual (values = c("AZA" = nejm_colors [5], "DEC" = nejm_colors [7], "GSK" = nejm_colors [4])) +
    scale_shape_manual(values = c("AZA" = 17, "DEC" = 15, "GSK" = 16)) +
    scale_x_log10() +
    scale_y_continuous(limits = c(0,125), breaks = c(seq(0,100,25))) + 
    geom_hline(aes(yintercept = 50), color = "black", linetype = "dashed") + 
    theme_classic(),
  # MOLT-4
  ggplot() +
    # AZA line and 95% CI
    geom_line(data = ndMOLT4AZA3days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [5]) +
    geom_ribbon(data = ndMOLT4AZA3days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [5]) +
    # DEC line and 95% CI
    geom_line(data = ndMOLT4DEC3days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [7]) +
    geom_ribbon(data = ndMOLT4DEC3days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [7]) +
    # GSK line and 95% CI
    geom_line(data = ndMOLT4GSK3days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [4]) +
    geom_ribbon(data = ndMOLT4GSK3days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [4]) +
    # mean viability and SEM at concentration that have been measured
    geom_point(data = viability_mean [which(viability_mean$treatment_duration == "3days" & viability_mean$cell_line == "MOLT4" & viability_mean$concentration != 0),],
               aes(x = concentration, y = viability_mean, shape = drug, colour = drug), size = 2.5) + 
    geom_errorbar(data = viability_mean [which(viability_mean$treatment_duration == "3days" & viability_mean$cell_line == "MOLT4" & viability_mean$concentration != 0),],
                  aes(x = concentration, y = viability_mean, ymin = viability_mean - sem, ymax = viability_mean + sem, color = drug), width = .05) +
    labs (x = "concentration [nM]", y = "Cell Viability (% of control)") + 
    # plot formatting (colors, titles, etc.)
    scale_color_manual (values = c("AZA" = nejm_colors [5], "DEC" = nejm_colors [7], "GSK" = nejm_colors [4])) +
    scale_shape_manual(values = c("AZA" = 17, "DEC" = 15, "GSK" = 16)) +
    scale_x_log10() +
    scale_y_continuous(limits = c(0,125), breaks = c(seq(0,100,25))) + 
    geom_hline(aes(yintercept = 50), color = "black", linetype = "dashed") + 
    theme_classic(),
  # PEER
  ggplot() +
    # AZA line and 95% CI
    geom_line(data = ndPEERAZA3days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [5]) +
    geom_ribbon(data = ndPEERAZA3days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [5]) +
    # DEC line and 95% CI
    geom_line(data = ndPEERDEC3days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [7]) +
    geom_ribbon(data = ndPEERDEC3days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [7]) +
    # GSK line and 95% CI
    geom_line(data = ndPEERGSK3days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [4]) +
    geom_ribbon(data = ndPEERGSK3days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [4]) +
    # mean viability and SEM at concentration that have been measured
    geom_point(data = viability_mean [which(viability_mean$treatment_duration == "3days" & viability_mean$cell_line == "PEER" & viability_mean$concentration != 0),],
               aes(x = concentration, y = viability_mean, shape = drug, colour = drug), size = 2.5) + 
    geom_errorbar(data = viability_mean [which(viability_mean$treatment_duration == "3days" & viability_mean$cell_line == "PEER" & viability_mean$concentration != 0),],
                  aes(x = concentration, y = viability_mean, ymin = viability_mean - sem, ymax = viability_mean + sem, color = drug), width = .05) +
    labs (x = "concentration [nM]", y = "Cell Viability (% of control)") + 
    # plot formatting (colors, titles, etc.)
    scale_color_manual (values = c("AZA" = nejm_colors [5], "DEC" = nejm_colors [7], "GSK" = nejm_colors [4])) +
    scale_shape_manual(values = c("AZA" = 17, "DEC" = 15, "GSK" = 16)) +
    scale_x_log10() +
    scale_y_continuous(limits = c(0,125), breaks = c(seq(0,100,25))) + 
    geom_hline(aes(yintercept = 50), color = "black", linetype = "dashed") + 
    theme_classic(),
  # SUP-T1
  ggplot() +
    # AZA line and 95% CI
    geom_line(data = ndSUPT1AZA3days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [5]) +
    geom_ribbon(data = ndSUPT1AZA3days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [5]) +
    # DEC line and 95% CI
    geom_line(data = ndSUPT1DEC3days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [7]) +
    geom_ribbon(data = ndSUPT1DEC3days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [7]) +
    # GSK line and 95% CI
    geom_line(data = ndSUPT1GSK3days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [4]) +
    geom_ribbon(data = ndSUPT1GSK3days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [4]) +
    # mean viability and SEM at concentration that have been measured
    geom_point(data = viability_mean [which(viability_mean$treatment_duration == "3days" & viability_mean$cell_line == "SUPT1" & viability_mean$concentration != 0),],
               aes(x = concentration, y = viability_mean, shape = drug, colour = drug), size = 2.5) + 
    geom_errorbar(data = viability_mean [which(viability_mean$treatment_duration == "3days" & viability_mean$cell_line == "SUPT1" & viability_mean$concentration != 0),],
                  aes(x = concentration, y = viability_mean, ymin = viability_mean - sem, ymax = viability_mean + sem, color = drug), width = .05) +
    labs (x = "concentration [nM]", y = "Cell Viability (% of control)") + 
    # plot formatting (colors, titles, etc.)
    scale_color_manual (values = c("AZA" = nejm_colors [5], "DEC" = nejm_colors [7], "GSK" = nejm_colors [4])) +
    scale_shape_manual(values = c("AZA" = 17, "DEC" = 15, "GSK" = 16)) +
    scale_x_log10() +
    scale_y_continuous(limits = c(0,125), breaks = c(seq(0,100,25))) + 
    geom_hline(aes(yintercept = 50), color = "black", linetype = "dashed") + 
    theme_classic(),
  # TALL-1
  ggplot() +
    # AZA line and 95% CI
    geom_line(data = ndTALL1AZA3days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [5]) +
    geom_ribbon(data = ndTALL1AZA3days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [5]) +
    # DEC line and 95% CI
    geom_line(data = ndTALL1DEC3days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [7]) +
    geom_ribbon(data = ndTALL1DEC3days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [7]) +
    # GSK line and 95% CI
    geom_line(data = ndTALL1GSK3days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [4]) +
    geom_ribbon(data = ndTALL1GSK3days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [4]) +
    # mean viability and SEM at concentration that have been measured
    geom_point(data = viability_mean [which(viability_mean$treatment_duration == "3days" & viability_mean$cell_line == "TALL1" & viability_mean$concentration != 0),],
               aes(x = concentration, y = viability_mean, shape = drug, colour = drug), size = 2.5) + 
    geom_errorbar(data = viability_mean [which(viability_mean$treatment_duration == "3days" & viability_mean$cell_line == "TALL1" & viability_mean$concentration != 0),],
                  aes(x = concentration, y = viability_mean, ymin = viability_mean - sem, ymax = viability_mean + sem, color = drug), width = .05) +
    labs (x = "concentration [nM]", y = "Cell Viability (% of control)") + 
    # plot formatting (colors, titles, etc.)
    scale_color_manual (values = c("AZA" = nejm_colors [5], "DEC" = nejm_colors [7], "GSK" = nejm_colors [4])) +
    scale_shape_manual(values = c("AZA" = 17, "DEC" = 15, "GSK" = 16)) +
    scale_x_log10() +
    scale_y_continuous(limits = c(0,125), breaks = c(seq(0,100,25))) + 
    geom_hline(aes(yintercept = 50), color = "black", linetype = "dashed") + 
    theme_classic(),
  
  labels = c("ALL-SIL", "CCRF-CEM", "DND41", "HPB-ALL", "JURKAT", "LOUCY", "MOLT-3", "MOLT-4", "PEER", "SUP-T1", "TALL-1"),
  label.y = 1.05, common.legend = TRUE, legend = "bottom",
  ncol = 6, nrow = 2 
)


### export plot

pdf("plots_for_publication/S1a_viability_HMAs_divided_by_CellLines_3days.pdf")
plot(plot_S1a)
dev.off()


### Extended Data Figure 1b: 7 days treatment
plot_S1b <- ggarrange (
  # ALL-SIL
  ggplot() +
    # AZA line and 95% CI
    geom_line(data = ndALLSILAZA7days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [5]) +
    geom_ribbon(data = ndALLSILAZA7days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [5]) +
    # DEC line and 95% CI
    geom_line(data = ndALLSILDEC7days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [7]) +
    geom_ribbon(data = ndALLSILDEC7days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [7]) +
    # GSK line and 95% CI
    geom_line(data = ndALLSILGSK7days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [4]) +
    geom_ribbon(data = ndALLSILGSK7days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [4]) +
    # mean viability and SEM at concentration that have been measured
    geom_point(data = viability_mean [which(viability_mean$treatment_duration == "7days" & viability_mean$cell_line == "ALLSIL" & viability_mean$concentration != 0),],
               aes(x = concentration, y = viability_mean, shape = drug, colour = drug), size = 2.5) + 
    geom_errorbar(data = viability_mean [which(viability_mean$treatment_duration == "7days" & viability_mean$cell_line == "ALLSIL" & viability_mean$concentration != 0),],
                  aes(x = concentration, y = viability_mean, ymin = viability_mean - sem, ymax = viability_mean + sem, color = drug), width = .05) +
    labs (x = "concentration [nM]", y = "Cell Viability (% of control)") + 
    # plot formatting (colors, titles, etc.)
    scale_color_manual (values = c("AZA" = nejm_colors [5], "DEC" = nejm_colors [7], "GSK" = nejm_colors [4])) +
    scale_shape_manual(values = c("AZA" = 17, "DEC" = 15, "GSK" = 16)) +
    scale_x_log10() +
    scale_y_continuous(limits = c(0,125), breaks = c(seq(0,100,25))) + 
    geom_hline(aes(yintercept = 50), color = "black", linetype = "dashed") + 
    theme_classic(),
  # CCRF-CEM
  ggplot() +
    # AZA line and 95% CI
    geom_line(data = ndCCRFCEMAZA7days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [5]) +
    geom_ribbon(data = ndCCRFCEMAZA7days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [5]) +
    # DEC line and 95% CI
    geom_line(data = ndCCRFCEMDEC7days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [7]) +
    geom_ribbon(data = ndCCRFCEMDEC7days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [7]) +
    # GSK line and 95% CI
    geom_line(data = ndCCRFCEMGSK7days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [4]) +
    geom_ribbon(data = ndCCRFCEMGSK7days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [4]) +
    # mean viability and SEM at concentration that have been measured
    geom_point(data = viability_mean [which(viability_mean$treatment_duration == "7days" & viability_mean$cell_line == "CCRFCEM" & viability_mean$concentration != 0),],
               aes(x = concentration, y = viability_mean, shape = drug, colour = drug), size = 2.5) + 
    geom_errorbar(data = viability_mean [which(viability_mean$treatment_duration == "7days" & viability_mean$cell_line == "CCRFCEM" & viability_mean$concentration != 0),],
                  aes(x = concentration, y = viability_mean, ymin = viability_mean - sem, ymax = viability_mean + sem, color = drug), width = .05) +
    labs (x = "concentration [nM]", y = "Cell Viability (% of control)") + 
    # plot formatting (colors, titles, etc.)
    scale_color_manual (values = c("AZA" = nejm_colors [5], "DEC" = nejm_colors [7], "GSK" = nejm_colors [4])) +
    scale_shape_manual(values = c("AZA" = 17, "DEC" = 15, "GSK" = 16)) +
    scale_x_log10() +
    scale_y_continuous(limits = c(0,125), breaks = c(seq(0,100,25))) + 
    geom_hline(aes(yintercept = 50), color = "black", linetype = "dashed") + 
    theme_classic(),
  # DND-41
  ggplot() +
    # AZA line and 95% CI
    geom_line(data = ndDND41AZA7days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [5]) +
    geom_ribbon(data = ndDND41AZA7days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [5]) +
    # DEC line and 95% CI
    geom_line(data = ndDND41DEC7days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [7]) +
    geom_ribbon(data = ndDND41DEC7days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [7]) +
    # GSK line and 95% CI
    geom_line(data = ndDND41GSK7days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [4]) +
    geom_ribbon(data = ndDND41GSK7days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [4]) +
    # mean viability and SEM at concentration that have been measured
    geom_point(data = viability_mean [which(viability_mean$treatment_duration == "7days" & viability_mean$cell_line == "DND41" & viability_mean$concentration != 0),],
               aes(x = concentration, y = viability_mean, shape = drug, colour = drug), size = 2.5) + 
    geom_errorbar(data = viability_mean [which(viability_mean$treatment_duration == "7days" & viability_mean$cell_line == "DND41" & viability_mean$concentration != 0),],
                  aes(x = concentration, y = viability_mean, ymin = viability_mean - sem, ymax = viability_mean + sem, color = drug), width = .05) +
    labs (x = "concentration [nM]", y = "Cell Viability (% of control)") + 
    # plot formatting (colors, titles, etc.)
    scale_color_manual (values = c("AZA" = nejm_colors [5], "DEC" = nejm_colors [7], "GSK" = nejm_colors [4])) +
    scale_shape_manual(values = c("AZA" = 17, "DEC" = 15, "GSK" = 16)) +
    scale_x_log10() +
    scale_y_continuous(limits = c(0,125), breaks = c(seq(0,100,25))) + 
    geom_hline(aes(yintercept = 50), color = "black", linetype = "dashed") + 
    theme_classic(),
  # HPB-ALL
  ggplot() +
    # AZA line and 95% CI
    geom_line(data = ndHPBALLAZA7days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [5]) +
    geom_ribbon(data = ndHPBALLAZA7days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [5]) +
    # DEC line and 95% CI
    geom_line(data = ndHPBALLDEC7days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [7]) +
    geom_ribbon(data = ndHPBALLDEC7days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [7]) +
    # GSK line and 95% CI
    geom_line(data = ndHPBALLGSK7days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [4]) +
    geom_ribbon(data = ndHPBALLGSK7days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [4]) +
    # mean viability and SEM at concentration that have been measured
    geom_point(data = viability_mean [which(viability_mean$treatment_duration == "7days" & viability_mean$cell_line == "HPBALL" & viability_mean$concentration != 0),],
               aes(x = concentration, y = viability_mean, shape = drug, colour = drug), size = 2.5) + 
    geom_errorbar(data = viability_mean [which(viability_mean$treatment_duration == "7days" & viability_mean$cell_line == "HPBALL" & viability_mean$concentration != 0),],
                  aes(x = concentration, y = viability_mean, ymin = viability_mean - sem, ymax = viability_mean + sem, color = drug), width = .05) +
    labs (x = "concentration [nM]", y = "Cell Viability (% of control)") + 
    # plot formatting (colors, titles, etc.)
    scale_color_manual (values = c("AZA" = nejm_colors [5], "DEC" = nejm_colors [7], "GSK" = nejm_colors [4])) +
    scale_shape_manual(values = c("AZA" = 17, "DEC" = 15, "GSK" = 16)) +
    scale_x_log10() +
    scale_y_continuous(limits = c(0,125), breaks = c(seq(0,100,25))) + 
    geom_hline(aes(yintercept = 50), color = "black", linetype = "dashed") + 
    theme_classic(),
  # JURKAT
  ggplot() +
    # AZA line and 95% CI
    geom_line(data = ndJURKATAZA7days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [5]) +
    geom_ribbon(data = ndJURKATAZA7days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [5]) +
    # DEC line and 95% CI
    geom_line(data = ndJURKATDEC7days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [7]) +
    geom_ribbon(data = ndJURKATDEC7days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [7]) +
    # GSK line and 95% CI
    geom_line(data = ndJURKATGSK7days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [4]) +
    geom_ribbon(data = ndJURKATGSK7days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [4]) +
    # mean viability and SEM at concentration that have been measured
    geom_point(data = viability_mean [which(viability_mean$treatment_duration == "7days" & viability_mean$cell_line == "JURKAT" & viability_mean$concentration != 0),],
               aes(x = concentration, y = viability_mean, shape = drug, colour = drug), size = 2.5) + 
    geom_errorbar(data = viability_mean [which(viability_mean$treatment_duration == "7days" & viability_mean$cell_line == "JURKAT" & viability_mean$concentration != 0),],
                  aes(x = concentration, y = viability_mean, ymin = viability_mean - sem, ymax = viability_mean + sem, color = drug), width = .05) +
    labs (x = "concentration [nM]", y = "Cell Viability (% of control)") + 
    # plot formatting (colors, titles, etc.)
    scale_color_manual (values = c("AZA" = nejm_colors [5], "DEC" = nejm_colors [7], "GSK" = nejm_colors [4])) +
    scale_shape_manual(values = c("AZA" = 17, "DEC" = 15, "GSK" = 16)) +
    scale_x_log10() +
    scale_y_continuous(limits = c(0,125), breaks = c(seq(0,100,25))) + 
    geom_hline(aes(yintercept = 50), color = "black", linetype = "dashed") + 
    theme_classic(),
  # LOUCY
  ggplot() +
    # AZA line and 95% CI
    geom_line(data = ndLOUCYAZA7days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [5]) +
    geom_ribbon(data = ndLOUCYAZA7days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [5]) +
    # DEC line and 95% CI
    geom_line(data = ndLOUCYDEC7days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [7]) +
    geom_ribbon(data = ndLOUCYDEC7days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [7]) +
    # GSK line and 95% CI
    geom_line(data = ndLOUCYGSK7days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [4]) +
    geom_ribbon(data = ndLOUCYGSK7days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [4]) +
    # mean viability and SEM at concentration that have been measured
    geom_point(data = viability_mean [which(viability_mean$treatment_duration == "7days" & viability_mean$cell_line == "LOUCY" & viability_mean$concentration != 0),],
               aes(x = concentration, y = viability_mean, shape = drug, colour = drug), size = 2.5) + 
    geom_errorbar(data = viability_mean [which(viability_mean$treatment_duration == "7days" & viability_mean$cell_line == "LOUCY" & viability_mean$concentration != 0),],
                  aes(x = concentration, y = viability_mean, ymin = viability_mean - sem, ymax = viability_mean + sem, color = drug), width = .05) +
    labs (x = "concentration [nM]", y = "Cell Viability (% of control)") + 
    # plot formatting (colors, titles, etc.)
    scale_color_manual (values = c("AZA" = nejm_colors [5], "DEC" = nejm_colors [7], "GSK" = nejm_colors [4])) +
    scale_shape_manual(values = c("AZA" = 17, "DEC" = 15, "GSK" = 16)) +
    scale_x_log10() +
    scale_y_continuous(limits = c(0,125), breaks = c(seq(0,100,25))) + 
    geom_hline(aes(yintercept = 50), color = "black", linetype = "dashed") + 
    theme_classic(),
  # MOLT-3
  ggplot() +
    # AZA line and 95% CI
    geom_line(data = ndMOLT3AZA7days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [5]) +
    geom_ribbon(data = ndMOLT3AZA7days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [5]) +
    # DEC line and 95% CI
    geom_line(data = ndMOLT3DEC7days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [7]) +
    geom_ribbon(data = ndMOLT3DEC7days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [7]) +
    # GSK line and 95% CI
    geom_line(data = ndMOLT3GSK7days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [4]) +
    geom_ribbon(data = ndMOLT3GSK7days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [4]) +
    # mean viability and SEM at concentration that have been measured
    geom_point(data = viability_mean [which(viability_mean$treatment_duration == "7days" & viability_mean$cell_line == "MOLT3" & viability_mean$concentration != 0),],
               aes(x = concentration, y = viability_mean, shape = drug, colour = drug), size = 2.5) + 
    geom_errorbar(data = viability_mean [which(viability_mean$treatment_duration == "7days" & viability_mean$cell_line == "MOLT3" & viability_mean$concentration != 0),],
                  aes(x = concentration, y = viability_mean, ymin = viability_mean - sem, ymax = viability_mean + sem, color = drug), width = .05) +
    labs (x = "concentration [nM]", y = "Cell Viability (% of control)") + 
    # plot formatting (colors, titles, etc.)
    scale_color_manual (values = c("AZA" = nejm_colors [5], "DEC" = nejm_colors [7], "GSK" = nejm_colors [4])) +
    scale_shape_manual(values = c("AZA" = 17, "DEC" = 15, "GSK" = 16)) +
    scale_x_log10() +
    scale_y_continuous(limits = c(0,125), breaks = c(seq(0,100,25))) + 
    geom_hline(aes(yintercept = 50), color = "black", linetype = "dashed") + 
    theme_classic(),
  # MOLT-4
  ggplot() +
    # AZA line and 95% CI
    geom_line(data = ndMOLT4AZA7days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [5]) +
    geom_ribbon(data = ndMOLT4AZA7days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [5]) +
    # DEC line and 95% CI
    geom_line(data = ndMOLT4DEC7days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [7]) +
    geom_ribbon(data = ndMOLT4DEC7days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [7]) +
    # GSK line and 95% CI
    geom_line(data = ndMOLT4GSK7days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [4]) +
    geom_ribbon(data = ndMOLT4GSK7days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [4]) +
    # mean viability and SEM at concentration that have been measured
    geom_point(data = viability_mean [which(viability_mean$treatment_duration == "7days" & viability_mean$cell_line == "MOLT4" & viability_mean$concentration != 0),],
               aes(x = concentration, y = viability_mean, shape = drug, colour = drug), size = 2.5) + 
    geom_errorbar(data = viability_mean [which(viability_mean$treatment_duration == "7days" & viability_mean$cell_line == "MOLT4" & viability_mean$concentration != 0),],
                  aes(x = concentration, y = viability_mean, ymin = viability_mean - sem, ymax = viability_mean + sem, color = drug), width = .05) +
    labs (x = "concentration [nM]", y = "Cell Viability (% of control)") + 
    # plot formatting (colors, titles, etc.)
    scale_color_manual (values = c("AZA" = nejm_colors [5], "DEC" = nejm_colors [7], "GSK" = nejm_colors [4])) +
    scale_shape_manual(values = c("AZA" = 17, "DEC" = 15, "GSK" = 16)) +
    scale_x_log10() +
    scale_y_continuous(limits = c(0,125), breaks = c(seq(0,100,25))) + 
    geom_hline(aes(yintercept = 50), color = "black", linetype = "dashed") + 
    theme_classic(),
  # PEER
  ggplot() +
    # AZA line and 95% CI
    geom_line(data = ndPEERAZA7days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [5]) +
    geom_ribbon(data = ndPEERAZA7days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [5]) +
    # DEC line and 95% CI
    geom_line(data = ndPEERDEC7days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [7]) +
    geom_ribbon(data = ndPEERDEC7days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [7]) +
    # GSK line and 95% CI
    geom_line(data = ndPEERGSK7days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [4]) +
    geom_ribbon(data = ndPEERGSK7days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [4]) +
    # mean viability and SEM at concentration that have been measured
    geom_point(data = viability_mean [which(viability_mean$treatment_duration == "7days" & viability_mean$cell_line == "PEER" & viability_mean$concentration != 0),],
               aes(x = concentration, y = viability_mean, shape = drug, colour = drug), size = 2.5) + 
    geom_errorbar(data = viability_mean [which(viability_mean$treatment_duration == "7days" & viability_mean$cell_line == "PEER" & viability_mean$concentration != 0),],
                  aes(x = concentration, y = viability_mean, ymin = viability_mean - sem, ymax = viability_mean + sem, color = drug), width = .05) +
    labs (x = "concentration [nM]", y = "Cell Viability (% of control)") + 
    # plot formatting (colors, titles, etc.)
    scale_color_manual (values = c("AZA" = nejm_colors [5], "DEC" = nejm_colors [7], "GSK" = nejm_colors [4])) +
    scale_shape_manual(values = c("AZA" = 17, "DEC" = 15, "GSK" = 16)) +
    scale_x_log10() +
    scale_y_continuous(limits = c(0,125), breaks = c(seq(0,100,25))) + 
    geom_hline(aes(yintercept = 50), color = "black", linetype = "dashed") + 
    theme_classic(),
  # SUP-T1
  ggplot() +
    # AZA line and 95% CI
    geom_line(data = ndSUPT1AZA7days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [5]) +
    geom_ribbon(data = ndSUPT1AZA7days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [5]) +
    # DEC line and 95% CI
    geom_line(data = ndSUPT1DEC7days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [7]) +
    geom_ribbon(data = ndSUPT1DEC7days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [7]) +
    # GSK line and 95% CI
    geom_line(data = ndSUPT1GSK7days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [4]) +
    geom_ribbon(data = ndSUPT1GSK7days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [4]) +
    # mean viability and SEM at concentration that have been measured
    geom_point(data = viability_mean [which(viability_mean$treatment_duration == "7days" & viability_mean$cell_line == "SUPT1" & viability_mean$concentration != 0),],
               aes(x = concentration, y = viability_mean, shape = drug, colour = drug), size = 2.5) + 
    geom_errorbar(data = viability_mean [which(viability_mean$treatment_duration == "7days" & viability_mean$cell_line == "SUPT1" & viability_mean$concentration != 0),],
                  aes(x = concentration, y = viability_mean, ymin = viability_mean - sem, ymax = viability_mean + sem, color = drug), width = .05) +
    labs (x = "concentration [nM]", y = "Cell Viability (% of control)") + 
    # plot formatting (colors, titles, etc.)
    scale_color_manual (values = c("AZA" = nejm_colors [5], "DEC" = nejm_colors [7], "GSK" = nejm_colors [4])) +
    scale_shape_manual(values = c("AZA" = 17, "DEC" = 15, "GSK" = 16)) +
    scale_x_log10() +
    scale_y_continuous(limits = c(0,125), breaks = c(seq(0,100,25))) + 
    geom_hline(aes(yintercept = 50), color = "black", linetype = "dashed") + 
    theme_classic(),
  # TALL-1
  ggplot() +
    # AZA line and 95% CI
    geom_line(data = ndTALL1AZA7days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [5]) +
    geom_ribbon(data = ndTALL1AZA7days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [5]) +
    # DEC line and 95% CI
    geom_line(data = ndTALL1DEC7days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [7]) +
    geom_ribbon(data = ndTALL1DEC7days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [7]) +
    # GSK line and 95% CI
    geom_line(data = ndTALL1GSK7days, aes(x = concentration, y = pred), linewidth = 1.2, color = nejm_colors [4]) +
    geom_ribbon(data = ndTALL1GSK7days, aes(x = concentration, y = pred, ymin = pmin, ymax = pmax), alpha = 0.2, fill = nejm_colors [4]) +
    # mean viability and SEM at concentration that have been measured
    geom_point(data = viability_mean [which(viability_mean$treatment_duration == "7days" & viability_mean$cell_line == "TALL1" & viability_mean$concentration != 0),],
               aes(x = concentration, y = viability_mean, shape = drug, colour = drug), size = 2.5) + 
    geom_errorbar(data = viability_mean [which(viability_mean$treatment_duration == "7days" & viability_mean$cell_line == "TALL1" & viability_mean$concentration != 0),],
                  aes(x = concentration, y = viability_mean, ymin = viability_mean - sem, ymax = viability_mean + sem, color = drug), width = .05) +
    labs (x = "concentration [nM]", y = "Cell Viability (% of control)") + 
    # plot formatting (colors, titles, etc.)
    scale_color_manual (values = c("AZA" = nejm_colors [5], "DEC" = nejm_colors [7], "GSK" = nejm_colors [4])) +
    scale_shape_manual(values = c("AZA" = 17, "DEC" = 15, "GSK" = 16)) +
    scale_x_log10() +
    scale_y_continuous(limits = c(0,125), breaks = c(seq(0,100,25))) + 
    geom_hline(aes(yintercept = 50), color = "black", linetype = "dashed") + 
    theme_classic(),
  labels = c("ALL-SIL", "CCRF-CEM", "DND41", "HPB-ALL", "JURKAT", "LOUCY", "MOLT-3", "MOLT-4", "PEER", "SUP-T1", "TALL-1"),
  label.y = 1.05, common.legend = TRUE, legend = "bottom",
  ncol = 6, nrow = 2 
)

pdf("plots_for_publication/S1b_viability_HMAs_divided_by_CellLines_7days.pdf")
plot(plot_S1b)
dev.off()



##### statistics #####


### comparison of response to drugs when combining all cell lines

# Mann-Whitney U/ Pairwise Wilcoxon Rank Sum test 
# with multiple testing correction (Benjamini Hochberg)
# Figure 2b

# 3 days treatment
wilcoxon_CellLines_combined_3days <- 
  pairwise.wilcox.test(x = viability_all[which(viability_all$treatment_duration == "3days"),"viability"],
                       g = viability_all[which(viability_all$treatment_duration == "3days"),"drug"], p.adjust.method = "BH")$p.value

# 7 days treatment
wilcoxon_CellLines_combined_7days <-
  pairwise.wilcox.test(x = viability_all[which(viability_all$treatment_duration == "7days"),"viability"],
                       g = viability_all[which(viability_all$treatment_duration == "7days"),"drug"], p.adjust.method = "BH")$p.value