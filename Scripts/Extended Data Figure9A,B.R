# Title: Extended data Fig9A-B

# Author: Sandra Hellberg
# Email: sandra.hellberg@liu.se

# Date: 2024-01-18
# Last modified: 2024-08-27

rm(list=ls()) # remove all entries in the global environment 
set.seed(2225) # seed for reproducibility

#-------------------------------------------------------------------------------

          ### Set directory structure and load packages---- 

#-------------------------------------------------------------------------------
setwd("add your working directory here")

lib_path <- .libPaths("C:/Rpackages/") # specifing the package library
# update.packages(lib.loc = lib_path, ask = FALSE, checkBuilt = TRUE, dependencies = TRUE)
# BiocManager::install(lib.loc = lib_path, update = TRUE, ask = FALSE)

pack_R <- c("dplyr", "tidyverse", "RColorBrewer", "ggplot2", "ggfortify", "ggpubr", "grid", "reshape",
            "ggsci", "scales", "vctrs") # libraries to load

for (i in 1:length(pack_R)) {
  library(pack_R[i], character.only = TRUE)
}

#-------------------------------------------------------------------------------

              ### Read in and format methylation data----

#-------------------------------------------------------------------------------

# Read the data from csv
meth_loucy <- read.delim("Supplementary_Table_15.csv", sep = ",") # LOUCY
meth_supt1 <- read.delim("Supplementary_Table_16.csv", sep = ",") # SUP-T1

# Function to filter rows by sample pattern and class_id type, and select specific columns
filter_samples <- function(data, sample_id, class_type) {
  data %>%
    filter(grepl(sample_id, sample), grepl(class_type, class_id)) %>%
    select(class_id, methylation_treated, methylation_control)
}

# Filter LOUCY dataset
LD3_meth_ERV <- filter_samples(meth_loucy, "LD3", "LTR|LTR\\?")
LD3_meth_LINE <- filter_samples(meth_loucy, "LD3", "LINE")
LD3_meth_SINE <- filter_samples(meth_loucy, "LD3", "SINE")

LG3_meth_ERV <- filter_samples(meth_loucy, "LG3", "LTR|LTR\\?")
LG3_meth_LINE <- filter_samples(meth_loucy, "LG3", "LINE")
LG3_meth_SINE <- filter_samples(meth_loucy, "LG3", "SINE")

LG7_meth_ERV <- filter_samples(meth_loucy, "LG7", "LTR|LTR\\?")
LG7_meth_LINE <- filter_samples(meth_loucy, "LG7", "LINE")
LG7_meth_SINE <- filter_samples(meth_loucy, "LG7", "SINE")

ERV_L <- cbind(LD3_meth_ERV, LG3_meth_ERV, LG7_meth_ERV) # combine into one data frame
ERV_L <- ERV_L[, -c(4,6,7,9)] # remove all but one control sample
ERV_L$class_id <- "ERV"
ERV_L <- ERV_L[, c(1,3,2,4,5)] # change order of the columns
colnames(ERV_L)[3:5] <- c("LD3", "LG3", "LG7")

LINE_L <- cbind(LD3_meth_LINE, LG3_meth_LINE, LG7_meth_LINE)
LINE_L <- LINE_L[, -c(4,6,7,9)]
LINE_L$class_id <- "LINE"
LINE_L <- LINE_L[, c(1,3,2,4,5)]
colnames(LINE_L)[3:5] <- c("LD3", "LG3", "LG7")

SINE_L <- cbind(LD3_meth_SINE, LG3_meth_SINE, LG7_meth_SINE)
SINE_L <- SINE_L[, -c(4,6,7,9)]
SINE_L$class_id <- "SINE"
SINE_L <- SINE_L[, c(1,3,2,4,5)]
colnames(SINE_L)[3:5] <- c("LD3", "LG3", "LG7")

# Filter SUP-T1 dataset
SD3_meth_ERV <- filter_samples(meth_supt1, "SD3", "LTR|LTR\\?")
SD3_meth_LINE <- filter_samples(meth_supt1, "SD3", "LINE")
SD3_meth_SINE <- filter_samples(meth_supt1, "SD3", "SINE")

SG3_meth_ERV <- filter_samples(meth_supt1, "SG3", "LTR|LTR\\?")
SG3_meth_LINE <- filter_samples(meth_supt1, "SG3", "LINE")
SG3_meth_SINE <- filter_samples(meth_supt1, "SG3", "SINE")

SG7_meth_ERV <- filter_samples(meth_supt1, "SG7","LTR|LTR\\?")
SG7_meth_LINE <- filter_samples(meth_supt1, "SG7", "LINE")
SG7_meth_SINE <- filter_samples(meth_supt1, "SG7", "SINE")

ERV_S <- cbind(SD3_meth_ERV, SG3_meth_ERV, SG7_meth_ERV)
ERV_S <- ERV_S[, -c(4,6,7,9)]
ERV_S$class_id <- "ERV"
ERV_S <- ERV_S[, c(1,3,2,4,5)]
colnames(ERV_S)[3:5] <- c("SD3", "SG3", "SG7")

LINE_S <- cbind(SD3_meth_LINE, SG3_meth_LINE, SG7_meth_LINE)
LINE_S <- LINE_S[, -c(4,6,7,9)]
LINE_S$class_id <- "LINE"
LINE_S <- LINE_S[, c(1,3,2,4,5)]
colnames(LINE_S)[3:5] <- c("SD3", "SG3", "SG7")

SINE_S <- cbind(SD3_meth_SINE, SG3_meth_SINE, SG7_meth_SINE)
SINE_S <- SINE_S[, -c(4,6,7,9)]
SINE_S$class_id <- "SINE"
SINE_S <- SINE_S[, c(1,3,2,4,5)]
colnames(SINE_S)[3:5] <- c("SD3", "SG3", "SG7")


#-------------------------------------------------------------------------------

                          ### Plotting----

#-------------------------------------------------------------------------------

#### Combing into one data frame per cell line for plotting
LOUCY <- rbind(ERV_L, LINE_L, SINE_L)
colnames(LOUCY) <- c("class_id", "untreated", "DEC", "GSK3", "GSK7")
SUPT1 <- rbind(ERV_S, LINE_S, SINE_S)
colnames(SUPT1) <- c("class_id", "untreated", "DEC", "GSK3", "GSK7")

LOUCY_melt <- melt.data.frame(LOUCY, id.vars = c("class_id")) # id.vars = variables to keep and not melt
SUPT1_melt <- melt.data.frame(SUPT1, id.vars = c("class_id")) # id.vars = variables to keep and not melt

#### Define color palette
npg_colors <- pal_npg()(10)
show_col (npg_colors)
nejm_colors <- pal_nejm()(8)
show_col (nejm_colors)


#### Defining pairwise comparisons
comparisons <- list(
  c("untreated", "DEC"),
  c("untreated", "GSK3"),
  c("untreated", "GSK7"))

#### Violin plot

LOUCY <- LOUCY_melt %>%
  ggplot(aes(x=variable, y=value)) +
  geom_violin(scale = "width", width=0.8) + 
  labs(y = "Methylation", x = "Class of TE") +
  geom_hline(yintercept = c(0.3, 0.7), col="black") +
  stat_compare_means(comparisons=comparisons, method="wilcox.test") +
  scale_y_continuous(expand = expansion(mult = .1)) +
  facet_wrap(~class_id, scales="free") +
  labs(title = "LOUCY") +
  theme_bw()


SUPT1 <- SUPT1_melt %>%
  ggplot(aes(x=variable, y=value)) +
  geom_violin(scale = "width", width=0.8) + 
  labs(y = "Methylation", x = "Class of TE") +
  geom_hline(yintercept = c(0.3, 0.7), col="black") +
  stat_compare_means(comparisons=comparisons, method="wilcox.test") +
  scale_y_continuous(expand = expansion(mult = .1)) +
  facet_wrap(~class_id, scales="free") +
  labs(title = "SUPT1") +
  theme_bw()


plot <- ggarrange(LOUCY, SUPT1, nrow=1, ncol=2)

pdf(paste0(fig_dir,"Extended data Fig9A and B.pdf"),
    width = 20,
    height = 10)
plot(plot)
dev.off()


#-------------------------------------------------------------------------------

                          ### Additional analysis----

#-------------------------------------------------------------------------------
# Calculate percent of TEs below 0.3
# Added in Extended Figure 9A and B

(sum(LD3_meth_ERV$methylation_control < 0.3)/nrow(LD3_meth_ERV)) * 100
(sum(LD3_meth_LINE$methylation_control < 0.3)/nrow(LD3_meth_LINE)) * 100
(sum(LD3_meth_SINE$methylation_control < 0.3)/nrow(LD3_meth_SINE)) * 100

# LOUCY
input <- list(LD3_meth_ERV, LG3_meth_ERV, LG7_meth_ERV,
              LD3_meth_LINE, LG3_meth_LINE, LG7_meth_LINE,
              LD3_meth_SINE, LG3_meth_SINE, LG7_meth_SINE)

for(i in 1:length(input)) {
  x <- input[[i]]
  print((sum(x[,2] < 0.3)/nrow(x)) * 100)
}

for(i in 1:length(input)) {
  x <- input[[i]]
  print((sum(x[,2] > 0.7)/nrow(x)) * 100)
}


# SUP-T1

(sum(SD3_meth_ERV$methylation_control < 0.3)/nrow(SD3_meth_ERV)) * 100
(sum(SG3_meth_LINE$methylation_control < 0.3)/nrow(SG3_meth_LINE)) * 100
(sum(SG7_meth_SINE$methylation_control < 0.3)/nrow(SG7_meth_SINE)) * 100

input <- list(SD3_meth_ERV, SG3_meth_ERV, SG7_meth_ERV,
              SD3_meth_LINE, SG3_meth_LINE, SG7_meth_LINE,
              SD3_meth_SINE, SG3_meth_SINE, SG7_meth_SINE)

for(i in 1:length(input)) {
  x <- input[[i]]
  print((sum(x[,2] < 0.3)/nrow(x)) * 100)
}

for(i in 1:length(input)) {
  x <- input[[i]]
  print((sum(x[,2] > 0.7)/nrow(x)) * 100)
}


# Calculate percent of TEs above 0.7
# Added in Extended Figure 9A and B

LD3 <- (sum(LD3_meth_ERV$methylation_control > 0.7)/nrow(LD3_meth_ERV))*100 
LG3 <- (sum(LD3_meth_LINE$methylation_control > 0.7)/nrow(LD3_meth_LINE))*100 
LG7 <- (sum(LD3_meth_SINE$methylation_control > 0.7)/nrow(LD3_meth_SINE))*100 

SD3 <- (sum(SD3_meth_ERV$methylation_control > 0.7)/nrow(SD3_meth_ERV))*100 
SG3 <- (sum(SG3_meth_LINE$methylation_control > 0.7)/nrow(SG3_meth_LINE))*100 
SG7 <- (sum(SG7_meth_SINE$methylation_control > 0.7)/nrow(SG7_meth_SINE))*100 




