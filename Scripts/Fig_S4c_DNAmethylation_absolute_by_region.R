################################################################################
###################### ABSOLUTE DNA METHYLATION BY REGION ######################
################################################################################

# Maike Bensberg (maike.bensberg@liu.se)
# 2024-07-11


# based on low-coverage whole methylome sequencing
# plot % CpG methylation (absolute) by region

# Figure: Extended Data Fig. 4c



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
                                    promoters_CGI, promoters_nonCGI, genes, ERVs, intergenic)



##### Extended Data Figure 4c #####


### plot DNA methylation by genomic region
# relative to control (untreated: 100% DNA methylation)
# 3 days treatment for all drugs and 7 days for GSK (in a separate plot)
# separate plot for each cell line


### colour defined by region (nejm colour pallette)
# CGI-promoters: nejm [6] (blue); nonCGI-promoters: nejm [5] (purple); gene bodies: nejm [4] (green); intergenic: nejm [7] (yellow); ERVs: nejm [8] (pink)

plot_S4c <-
  ggplot(treated_global_CpGme_melt[which(treated_global_CpGme_melt$drug == "GSK" &
                                           treated_global_CpGme_melt$concentration %in% c(0,3000) &
                                           treated_global_CpGme_melt$treatment_duration == "7days" &
                                           treated_global_CpGme_melt$region %in% c("promoters_CGI", "promoters_nonCGI", "genes", "intergenic", "ERVs")),],
         aes(x = concentration, y = methylation, color = region)) +
  geom_point(size = 4) +
  geom_line(linewidth = 1) +
  scale_color_manual(breaks = c("promoters_CGI", "promoters_nonCGI", "genes", "intergenic", "ERVs"),
                     values = c(nejm_colors [6],nejm_colors [5],nejm_colors [4],nejm_colors [7],nejm_colors [8])) +
  scale_y_continuous(limits = c(0,100), breaks = c(0,25,50,75,100)) +
  labs(title = "Treatment for 7 days with 3000 nM GSK", y = "% global CpG methylation") +
  facet_wrap(~cell_line, ncol = 4) +
  theme_bw()

# export plot
pdf("plots_for_publication/S4c_DNAmethylation_by_region_GSK3000.pdf")
plot(plot_S4c)
dev.off()