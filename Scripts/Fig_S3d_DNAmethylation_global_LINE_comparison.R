################################################################################
####################### DNA METHYLATION COMPARED TO LINES ###################### 
################################################################################

# Maike Bensberg (maike.bensberg@liu.se)
# 2024-07-12


# based on low-coverage whole methylome sequencing
# compare global % CpG methylation to CpG methylation at LINE elements
# DNA methylation at LINE elements has been used to estimate global methylation

# Figure: Extended Data Fig. 3d



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


### data is generated based on analysis of low-coverage whole methylome sequencing
# two scripts were used for analysis of fastq seqeuncing files
# scripts: LowPassMeth_preprocess_MB.sh & LowPassMeth_downstream.sh 


### DNA methylation of untreated cell T-ALL cell lines in replicates
# % CpG methylation globally (CpGmethylation) and for specified regions of interest (Supplementary Table 2)
# regions of interest based on location relative to genes but here we will only look at LINEs
# column names: "cell_line","drug","concentration","treatment_duration","sample_id"
#   "CpGmethylation","CGI","non_CGI","promoters_all","promoters_CGI","promoters_nonCGI","genes","intergenic","LINEs","SINEs","ERVs"
untreated_global_CpGme <- read.csv("data/Supplementary_Table_2.csv")


### calculate mean from replicates
untreated_global_CpGme_mean <- 
  untreated_global_CpGme %>% group_by(cell_line) %>% 
  summarise(globalCpGmethylation = mean(CpGmethylation_global), LINEs = mean(LINEs))



##### Figure S3d #####


### colours.
# global DNA methylation: black; LINE elements: grey

plot_S3d <- ggarrange(
  ggplot(data = untreated_global_CpGme, aes(x = reorder(sample_id, CpGmethylation_global))) +
    geom_point(aes (y = CpGmethylation_global), color = "black", fill = "black", alpha = 0.8, size = 3) +
    geom_point(aes (y = LINEs), color = "azure3", fill = "azure3", alpha = 0.8, size = 3) +
    labs(title = "LINEs (grey) and globally (black)") +
    scale_color_manual(values = c(nejm_colors, nejm_colors)) +
    scale_y_continuous(limits = c(60,100), breaks = c(60,70,80,90,100)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)),
  ggplot(data = untreated_global_CpGme_mean, aes(x = reorder(cell_line, globalCpGmethylation))) +
    geom_point(aes (y = globalCpGmethylation), color = "black", fill = "black", alpha = 0.8, size = 3) +
    geom_point(aes (y = LINEs), color = "azure3", fill = "azure3", alpha = 0.8, size = 3) +
    labs(title = "LINEs (grey) and globally (black) mean") +
    scale_color_manual(values = c(nejm_colors, nejm_colors)) +
    scale_y_continuous(limits = c(60,100), breaks = c(60,70,80,90,100)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)),
  ncol = 2, align = "h"
)

# export plot
pdf("plots_for_publication/S3d_DNAmethylation_global_LINEs.pdf")
plot(plot_S3d)
dev.off()