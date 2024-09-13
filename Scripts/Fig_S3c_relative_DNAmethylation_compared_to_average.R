################################################################################
############# DNA METHYLATION RELATIV TO AVERAGE ACROSS CELL LINES #############
################################################################################

# Maike Bensberg (maike.bensberg@liu.se)
# 2024-07-12


# based on low-coverage whole methylome sequencing
# compare % CpG methylation for each cell line to average across all cell lines
# % CpG methylation / average of all cell lines (plot log2)

# Figure: Extended Data Fig. 3c



##### load packages #####

library(ggplot2)
library(ggsci)
library(scales)
library(dplyr)



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


### DNA methylation of untreated cell T-ALL cell lines in replicates
# % CpG methylation globally (CpGmethylation) and for specified regions of interest (Supplementary Table 2)
# regions of interest based on location relative to genes but here we will only look at CGIs
# column names: "cell_line","drug","concentration","treatment_duration","sample_id"
#   "CpGmethylation","CGI","non_CGI","promoters_all","promoters_CGI","promoters_nonCGI","genes","intergenic","LINEs","SINEs","ERVs"
untreated_global_CpGme <- read.csv("data/Supplementary_Table_2.csv")


### calculate mean from replicates
untreated_global_CpGme_mean <- 
  untreated_global_CpGme %>% group_by(cell_line) %>% 
  summarise(globalCpGmethylation = mean(CpGmethylation_global), CGI = mean(CGI))


### calculate log2 of relative difference of individual cell line to mean
untreated_global_CpGme_mean_log2ofmean <- untreated_global_CpGme_mean [,1]
untreated_global_CpGme_mean_log2ofmean$globalCpGmethylation_log2 <- log2(untreated_global_CpGme_mean$globalCpGmethylation/
                                                                           mean(untreated_global_CpGme_mean$globalCpGmethylation))
untreated_global_CpGme_mean_log2ofmean$CGI_log2 <- log2(untreated_global_CpGme_mean$CGI/
                                                          mean(untreated_global_CpGme_mean$CGI))



##### Figure S3c #####


### plot bar plot
# x-axis sorted by global methylation (ascending)


### colour defined by region (colour pallettes)
# CGI: nejm [2] (blue); non-CGI/global: nejm [3] (orange)

plot_S3c <- 
  ggplot(untreated_global_CpGme_mean_log2ofmean, aes(x = reorder(cell_line, globalCpGmethylation_log2))) +
    geom_bar(aes(y = CGI_log2), stat = "identity",
             colour = nejm_colors [2], fill = nejm_colors [2], alpha = 0.8,
             width = 0.35, position = position_nudge(x = -0.2)) +
    geom_bar(aes(y = globalCpGmethylation_log2 ), stat = "identity",
             colour = nejm_colors [3], fill = nejm_colors [3], alpha = 0.8,
             width = 0.35, position = position_nudge(x = 0.2)) +
    labs(title = "methylation relative to average of all cell lines (log2)",
         subtitle = "CGI (blue), global (oragne)",
         y = "log2 (methylation / average)", x = "cell line") +
    scale_fill_manual(values = c(nejm_colors, nejm_colors)) +
    theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# export plot
pdf("plots_for_publication/S3c_DNAmethylation_relative_average_across_cell_lines.pdf")
plot(plot_S3c)
dev.off()