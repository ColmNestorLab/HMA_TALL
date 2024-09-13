################################################################################
########## PIE CHART: distribution of regions keeping DNA methylation ########## 
################################################################################

# Maike Bensberg (maike.bensberg@liu.se)
# 2024-08-21


### plot a pie chart showing the distribution of regions ḱeeping DNA methylation
# across genomic locations (promoter > exon > intron > intergenic)
# only looked at regions that are 80% methylated in untreated and treated cells
# treated cells: treated with GSK for 7 days


### data analysis see script regions_retaining_methylation_always_80.R

### Figure: Fig. 4b



##### initiate R environment #####
# only needs to be done once
#install.packages("renv")
#renv::init()



##### loading packages #####

library(ggsci)
library(scales)
library(ggplot2)



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

# import
regions_80_always_annot_summary_count <- read.delim("R_exports/tables/regions_80_always_annot_summary_count.txt")

# make dataframe for plotting
df_for_plotting <- melt(regions_80_always_annot_summary_count[ which(regions_80_always_annot_summary_count$region != "TE"),],
                        id.vars = "region", variable.name = "sample", value.name = "percentage")
df_for_plotting$order <- rep(4:1,4)



##### plot #####


### pie chart
# include total number of regions overlapping a given genomic location
# comparing untreated and treated (GSK 7 days) cells
# separate chart for 


### fill colour defined by region
# promoter: nejm [5] (purple), exons: npg [3] (intense green),
# introns: npg [7] (light green), intergenic: nejm [7] (yellow)


plot <-
  ggplot(data = df_for_plotting, aes(x = "", y = percentage, fill = reorder(region, order))) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  geom_text(aes(label = round(percentage, digits = 2)), position = position_stack(vjust = 0.5)) +
  facet_wrap(~sample, scales = "free") +
  scale_fill_manual(values = c("promoter" = nejm_colors[5], "exon" = npg_colors[3],
                               "intron" = npg_colors[7], "intergenic" = nejm_colors[7])) +
  labs(title = "Number of regions that were >= 80% methylated in untreated and treated samples") +
  theme_void()

# export plot
pdf("R_exports/plots/Fig_4b_pie_regions_keeping_methylation_by_location.pdf")
plot(plot)
dev.off()

remove(plot)
remove(df_for_plotting)