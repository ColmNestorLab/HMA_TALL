################################################################################
############### DNA METHYLATION OF UNTREATED CELL LINES: AVERAGE ###############
################################################################################

# Maike Bensberg (maike.bensberg@liu.se)
# 2024-07-12


# based on low-coverage whole methylome sequencing
# plot % CpG methylation of cell lines (mean of replicates)
# divided into global and methylation at CpG islands

# Figure: Extended Data Fig. 3b



##### load packages #####

library(ggpubr)
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

npg_colors <- pal_npg()(10)
show_col (npg_colors)
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



##### Figure S3b #####


### points: global and CGI DNA methylation for individual replicates


### define colours
# mostly based on pre-defined colour pallettes
CL_colours <- c("ALLSIL" = nejm_colors[1], "CCRFCEM" = nejm_colors[2], "DND41" = "#9d6199", "HPBALL" = nejm_colors[4],
                "JURKAT" = nejm_colors[5], "LOUCY" =  nejm_colors[6], "MOLT3" = nejm_colors[7], "MOLT4" = nejm_colors[8],
                "PEER" = npg_colors[9], "SUPT1" = npg_colors[5], "TALL1" = nejm_colors[3])


plot_S3b <- ggarrange(
  # global CpG methylation
  ggplot(data = untreated_global_CpGme_mean, aes(x = reorder(cell_line, globalCpGmethylation), y = globalCpGmethylation, color = cell_line)) +
    geom_point(size = 3) +
    labs(title = "Global CpG methylation") +
    scale_color_manual(values = CL_colours) +
    scale_y_continuous(limits = c(58,92), breaks = c(50,60,70,80,90)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)),
  # CGI methylation
  ggplot(data = untreated_global_CpGme_mean, aes(x = reorder(cell_line, CGI), y = CGI, color = cell_line)) +
    geom_point(size = 3) +
    labs(title = "CGI methylation") +
    scale_color_manual(values = CL_colours) +
    scale_y_continuous(limits = c(28,62), breaks = c(30,40,50,60)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)),
  nrow = 1, ncol = 2, common.legend = TRUE
)

remove(CL_colours)

# export plot
pdf("plots_for_publication/S3b_DNAmethylation_cell_lines_average.pdf")
plot(plot_S3b)
dev.off()