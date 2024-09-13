################################################################################
######################### RELATIVE DNA METHYLATION PLOT ########################
################################################################################

# Maike Bensberg (maike.bensberg@liu.se)
# 2024-07-11


# based on low-coverage whole methylome sequencing
# calculate % CpG methylation relative to DNA methylation in control cells
# plot as line plot with points

# Figures: Fig. 3b (right), Fig. 3d, Extended Data Fig. 4b
# global CpG methylation and divided by genomic regions



##### load packages #####

library(ggpubr)
library(ggsci)
library(scales)
library(tidyr)



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


#### import and format data #####


### data is generated based on analysis of low-coverage whole methylome sequencing
# two scripts were used for analysis of fastq seqeuncing files
# scripts: LowPassMeth_preprocess_MB.sh & LowPassMeth_downstream.sh 


### cell lines treated with HMAs
# % CpG methylation globally (CpGmethylation) and for specified regions of interest (Supplementary Table 5)
# regions of interest based on location relative to genes: promoters, gene bodies, intergenic
# also looks at transposable elements: LINEs, SINEs, ERVs
# column names: "cell_line","drug","concentration","treatment_duration","sample_id"
# "CpGmethylation","CGI","non_CGI","promoters_all","promoters_CGI","promoters_nonCGI","genes","intergenic","LINEs","SINEs","ERVs"
treated_global_CpGme <- read.csv("data/Supplementary_Table_5.csv")


### format data for plotting

treated_global_CpGme2 <- treated_global_CpGme

# change drug "control" to "GSK" and add rows for AZA and DEC controls
temp <- treated_global_CpGme2 [which(treated_global_CpGme2$drug == "control"),]
temp$drug <- "AZA"
treated_global_CpGme2 <- rbind(treated_global_CpGme2, temp)
temp$drug <- "DEC"
treated_global_CpGme2 <- rbind(treated_global_CpGme2, temp)
remove(temp)
treated_global_CpGme2 [which(treated_global_CpGme2$drug == "control"),]$drug <- "GSK"

# add a treatment column: combination of treatment duration and concentration
treated_global_CpGme2$treatment <- paste(treated_global_CpGme2$treatment_duration, treated_global_CpGme2$concentration)

# melt dataframe based on region
treated_global_CpGme_melt <- gather(treated_global_CpGme2 [,c(1:4,6:17)], key = "region", value = "methylation", 
                                    CpGmethylation, CGI, non_CGI, promoters_CGI, promoters_nonCGI, genes, intergenic, ERVs)

# make a new dataframe showing global CpG methylation relative to initial methylation (untreated)
treated_global_CpGme2_relative <- treated_global_CpGme_melt

for (r in unique(treated_global_CpGme2_relative$region)) {
  for (CL in unique(treated_global_CpGme2_relative$cell_line)) {
    for (d in unique(treated_global_CpGme2_relative$treatment_duration)) {
      treated_global_CpGme2_relative[which(treated_global_CpGme2_relative$cell_line == CL & treated_global_CpGme2_relative$region == r
                                           & treated_global_CpGme2_relative$treatment_duration == d),]$methylation <- 
        (treated_global_CpGme2_relative[which(treated_global_CpGme2_relative$cell_line == CL & treated_global_CpGme2_relative$region == r
                                              & treated_global_CpGme2_relative$treatment_duration == d),]$methylation)/
        (treated_global_CpGme2_relative[which(treated_global_CpGme2_relative$cell_line == CL & treated_global_CpGme2_relative$region == r
                                              & treated_global_CpGme2_relative$treatment_duration == d
                                              & treated_global_CpGme2_relative$concentration == 0),]$methylation) * 100
    }
  }
}
remove(r)
remove(CL)
remove(d)



##### plot lineplots #####


### Figure 3b (right)


### plot global DNA methylation
# relative to control (untreated: 100% DNA methylation)
# only for 3 days treatment
# separate plot for each cell line


### colour, shape and linetype defined by treatment
# colours (pre-defined colour pallette): AZA: nejm [5] (purple); DEC: nejm [7] (yellow); GSK: nejm [4] (green)
# shapes: AZA: triangle; DEC: square; GSK: circle
# linetype: AZA: dotted; DEC: dashed; GSK: solid

plot_3b2 <-
  ggplot(treated_global_CpGme2_relative [which(treated_global_CpGme2_relative$treatment_duration == "3days"
                                                     & treated_global_CpGme2_relative$concentration != 0
                                                     & treated_global_CpGme2_relative$region == "CpGmethylation"),],
               aes(x = concentration, y = methylation, color = drug)) +
  geom_point(aes(shape = drug), size = 3) +
  geom_line(aes(linetype = drug), linewidth = 1) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("AZA" = nejm_colors [5], "DEC" =  nejm_colors [7], "GSK" =  nejm_colors [4])) +
  scale_shape_manual (values = c("AZA" = 17, "DEC" = 15, "GSK" = 16)) +
  scale_linetype_manual (values = c("AZA" = "dotted", "DEC" = "dashed", "GSK" = "solid")) +
  scale_x_log10() +
  scale_y_continuous(limits = c(0,110), breaks = c(0,25,50,75,100)) +
  facet_wrap(~cell_line) +
  labs (title = "methylation relative to control: 3 days treatment", y = "DNA methylation (relative to control)") +
  theme_bw()

pdf("plots_for_publication/3b2_relative_global_DNAmethylation_3days.pdf")
plot(plot_3b2)
dev.off()


### Figure 3d 


### plot DNA methylation by genomic region (CGI and non-CGI)
# relative to control (untreated: 100% DNA methylation)
# 3 days treatment for all drugs and 7 days for GSK (in a separate plot)
# separate plot for each cell line and each drug


### colour defined by region (colour pallettes)
# CGI: nejm [2] (blue); non-CGI/global: nejm [3] (orange)


### shape and linetype defined by treatment
# shapes: AZA: triangle; DEC: square; GSK: circle
# linetype: AZA: dotted; DEC: dashed; GSK: solid

plot_3d <- ggarrange(
  ggplot(treated_global_CpGme2_relative [which(treated_global_CpGme2_relative$treatment_duration == "3days"
                                               & treated_global_CpGme2_relative$concentration != 0
                                               & treated_global_CpGme2_relative$region %in% c("CGI", "non_CGI")),],
         aes(x = concentration, y = methylation, color = region)) +
    geom_point(aes(shape = drug), size = 3) +
    geom_line(aes(linetype = drug), linewidth = 1) +
    geom_hline(yintercept = 50, linetype = "dashed", color = "black") +
    scale_color_manual(breaks = c("CGI", "non_CGI"),
                       values = c(nejm_colors [2],nejm_colors [3])) +
    scale_shape_manual (values = c("AZA" = 17, "DEC" = 15, "GSK" = 16)) +
    scale_linetype_manual (values = c("AZA" = "dotted", "DEC" = "dashed", "GSK" = "solid")) +
    scale_x_log10() +
    scale_y_continuous(limits = c(0,110), breaks = c(0,25,50,75,100)) +
    facet_wrap(drug~cell_line) +
    labs (title = "methylation relative to control at different regions: 3 days treatment", y = "DNA methylation (relative to control)") +
    theme_bw(),
  ggplot(treated_global_CpGme2_relative [which(treated_global_CpGme2_relative$treatment_duration == "7days"
                                               & treated_global_CpGme2_relative$drug == "GSK"
                                               & treated_global_CpGme2_relative$concentration != 0
                                               & treated_global_CpGme2_relative$region %in% c("CGI", "non_CGI")),],
         aes(x = concentration, y = methylation, color = region)) +
    geom_point(size = 3) +
    geom_line(linewidth = 1) +
    scale_color_manual(breaks = c("CGI", "non_CGI"),
                       values = c(nejm_colors [2],nejm_colors [3])) +
    scale_x_log10() +
    scale_y_continuous(limits = c(0,110), breaks = c(0,25,50,75,100)) +
    facet_wrap(~cell_line, ncol = 4) +
    labs (title = "methylation relative to control at different regions: 7 days treatment", y = "DNA methylation (relative to control)") +
    theme_bw(),
  align = "v", common.legend = TRUE, ncol = 1, heights = c(2.63,1)
)

pdf("plots_for_publication/S4b_relative_DNAmethylation_by_CGI.pdf")
plot(plot_3d)
dev.off()


### Extended Data Figure 4b 


### plot DNA methylation by genomic region
# relative to control (untreated: 100% DNA methylation)
# 3 days treatment for all drugs and 7 days for GSK (in a separate plot)
# separate plot for each cell line


### colour defined by region (nejm colour pallette)
# CGI-promoters: nejm [6] (blue); nonCGI-promoters: nejm [5] (purple); gene bodies: nejm [4] (green); intergenic: nejm [7] (yellow); ERVs: nejm [8] (pink)


### shape and linetype defined by treatment
# shapes: AZA: triangle; DEC: square; GSK: circle
# linetype: AZA: dotted; DEC: dashed; GSK: solid

plot_S4b <- ggarrange(
  ggplot(treated_global_CpGme2_relative [which(treated_global_CpGme2_relative$treatment_duration == "3days"
                                               & treated_global_CpGme2_relative$concentration != 0
                                               & treated_global_CpGme2_relative$region %in%c("promoters_CGI", "promoters_nonCGI", "genes", "intergenic", "ERVs")),],
         aes(x = concentration, y = methylation, color = region)) +
    geom_point(aes(shape = drug), size = 3) +
    geom_line(aes(linetype = drug), linewidth = 1) +
    geom_hline(yintercept = 50, linetype = "dashed", color = "black") +
    scale_color_manual(breaks = c("promoters_CGI", "promoters_nonCGI", "genes", "intergenic", "ERVs"),
                       values = c(nejm_colors [6],nejm_colors [5],nejm_colors [4],nejm_colors [7],nejm_colors [8])) +
    scale_shape_manual (values = c("AZA" = 17, "DEC" = 15, "GSK" = 16)) +
    scale_linetype_manual (values = c("AZA" = "dotted", "DEC" = "dashed", "GSK" = "solid")) +
    scale_x_log10() +
    scale_y_continuous(limits = c(0,110), breaks = c(0,25,50,75,100)) +
    facet_wrap(drug~cell_line) +
    labs (title = "methylation relative to control at different regions: 3 days treatment", y = "DNA methylation (relative to control)") +
    theme_bw(),
  ggplot(treated_global_CpGme2_relative [which(treated_global_CpGme2_relative$treatment_duration == "7days"
                                               & treated_global_CpGme2_relative$drug == "GSK"
                                               & treated_global_CpGme2_relative$concentration != 0
                                               & treated_global_CpGme2_relative$region %in% c("promoters_CGI", "promoters_nonCGI", "genes", "intergenic", "ERVs")),],
         aes(x = concentration, y = methylation, color = region)) +
    geom_point(size = 3) +
    geom_line(linewidth = 1) +
    scale_color_manual(breaks = c("promoters_CGI", "promoters_nonCGI", "genes", "intergenic", "ERVs"),
                       values = c(nejm_colors [6],nejm_colors [5],nejm_colors [4],nejm_colors [7],nejm_colors [8])) +
    scale_x_log10() +
    scale_y_continuous(limits = c(0,110), breaks = c(0,25,50,75,100)) +
    facet_wrap(~cell_line, ncol = 4) +
    labs (title = "methylation relative to control at different regions: 7 days treatment", y = "DNA methylation (relative to control)") +
    theme_bw(),
  align = "v", common.legend = TRUE, ncol = 1, heights = c(2.63,1)
)

pdf("plots_for_publication/S4b_relative_DNAmethylation_by_region.pdf")
plot(plot_S4b)
dev.off()