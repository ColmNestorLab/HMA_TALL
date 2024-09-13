################################################################################
####################### ABSOLUTE DNA METHYLATION HEATMAP #######################
################################################################################

# Maike Bensberg (maike.bensberg@liu.se)
# 2024-07-11


# based on low-coverage whole methylome sequencing
# plot DNA methylation (% CpG methylation) as heatmap
# include % CpG methylation values as text

# Figures: Fig. 3b (left), Extended Data Fig. 4a
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
npg_colors <- pal_npg()(10)
show_col (npg_colors)



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
                                    CpGmethylation, CGI, non_CGI, promoters_all, promoters_CGI, promoters_nonCGI, genes, intergenic, LINEs, SINEs, ERVs)




##### plot heatmaps #####


### tile colour based on % Cpg methylation
# colour gradient from white to blue (0 - 100% CpG methylation)
# blue colour based on pre-defined colour pallette: npg_colors [4] 


### Figure 3b (left)


### heatmap of global CpG methylation
# one heatmap for each cell line
# x-axis: drug concentration and treatment duration
# y-axis: drug

plot_3b1 <-
  ggplot(treated_global_CpGme2, aes(x = treatment, y = drug, fill = CpGmethylation)) +
  geom_tile(color = "black") +
  geom_text(aes(label = round(CpGmethylation, digits = 0)), size = 2) +
  scale_x_discrete(limits = c("3days 0", "3days 10", "3days 30", "3days 100", "3days 300", "3days 1000", "3days 3000",
                              "7days 0", "7days 10", "7days 30", "7days 100", "7days 300", "7days 1000", "7days 3000")) +
  scale_y_discrete(limits = c("GSK", "DEC", "AZA")) +
  labs(title = "global CpG methylation", x = "treatment duration and concentration", y = "drug") +
  scale_fill_gradient(low = "white", high = npg_colors[4], limits = c(0,100)) +
  facet_wrap(~cell_line, ncol = 1) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


### export plot

pdf("plots_for_publication/3b1_heatmap_treated_globalCpG.pdf")
plot(plot_3b1)
dev.off()


### Extended Data Figure 4a


### heatmap of CpG methylation by region
# one heatmap for each cell line and each drug
# x-axis: drug concentration and treatment duration
# y-axis: region


### AZA
plot_S4a_AZA <-
  ggplot(treated_global_CpGme_melt [which(treated_global_CpGme_melt$drug == "AZA"),], aes(x = treatment, y = region, fill = methylation)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "white", high = npg_colors[4], limits = c(0,100)) +
  geom_text(aes(label = round(methylation, digits = 0)), size = 2) + 
  scale_x_discrete(limits = c("3days 0", "3days 10", "3days 30", "3days 100", "3days 300", "3days 1000", "3days 3000",
                              "7days 0", "7days 10", "7days 30")) +
  scale_y_discrete(limits = c("promoters_CGI", "promoters_nonCGI", "genes", "intergenic","ERVs", "LINEs", "SINEs")) +
  facet_wrap(~cell_line) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# export
pdf("plots_for_publication/S4a_heatmap_treated_regions_AZA.pdf")
plot(plot_S4a_AZA)
dev.off()


### DEC
plot_S4a_DEC <-
  ggplot(treated_global_CpGme_melt [which(treated_global_CpGme_melt$drug == "DEC"),], aes(x = treatment, y = region, fill = methylation)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "white", high = npg_colors[4], limits = c(0,100)) +
  geom_text(aes(label = round(methylation, digits = 0)), size = 2) + 
  scale_x_discrete(limits = c("3days 0", "3days 10", "3days 30", "3days 100", "3days 300", "3days 1000", "3days 3000",
                              "7days 0", "7days 10", "7days 30")) +
  scale_y_discrete(limits = c("promoters_CGI", "promoters_nonCGI", "genes", "intergenic","ERVs", "LINEs", "SINEs")) +
  facet_wrap(~cell_line) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# export
pdf("plots_for_publication/S4a_heatmap_treated_regions_DEC.pdf")
plot(plot_S4a_DEC)
dev.off()


### GSK
plot_S4a_GSK <-
  ggplot(treated_global_CpGme_melt [which(treated_global_CpGme_melt$drug == "GSK"),], aes(x = treatment, y = region, fill = methylation)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "white", high = npg_colors[4], limits = c(0,100)) +
  geom_text(aes(label = round(methylation, digits = 0)), size = 2) + 
  scale_x_discrete(limits = c("3days 0", "3days 10", "3days 30", "3days 100", "3days 300", "3days 1000", "3days 3000",
                              "7days 0", "7days 10", "7days 30", "7days 100", "7days 300", "7days 1000", "7days 3000")) +
  scale_y_discrete(limits = c("promoters_CGI", "promoters_nonCGI", "genes", "intergenic","ERVs", "LINEs", "SINEs")) +
  facet_wrap(~cell_line) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# export
pdf("plots_for_publication/S4a_heatmap_treated_regions_GSK.pdf")
plot(plot_S4a_GSK)
dev.off()