# Title: Figure 5B
# Author: Sandra Hellberg
# Email: sandra.hellberg@liu.se

# Date: 2024-01-22
# Last modified: 2024-08-14

rm(list=ls()) # remove all entries in the global environment 
set.seed(568) # seed for reproducibility

#-------------------------------------------------------------------------------
                    
              ### Set directory structure and load packages 

#-------------------------------------------------------------------------------
setwd("set your working directory here")

lib_path <- .libPaths("C:/Rpackages/") # specifing the package library
#update.packages(lib.loc = lib_path, ask = FALSE, checkBuilt = TRUE, dependencies = TRUE)
#BiocManager::install(lib.loc = lib_path, update = TRUE, ask = FALSE)

pack_R <- c("dplyr", "tidyverse", "RColorBrewer", "ggplot2", "ggfortify", "ggpubr", "grid", "reshape",
            "ggsci", "scales") # libraries to load

for (i in 1:length(pack_R)) {
  library(pack_R[i], character.only = TRUE)
}

#-------------------------------------------------------------------------------

                      ## Read in and format methylation data----

#-------------------------------------------------------------------------------
# Read the data from csv
meth_loucy <- read.delim("Supplementary_Table_15.csv", sep = ",") # LOUCY
meth_supt1 <- read.delim("Supplementary_Table_16.csv", sep = ",") # SUP-T1

# Function to filter rows by sample pattern and class_id type, and select specific columns
filter_samples <- function(data, sample_id, class_type) {
  data %>%
    filter(grepl(sample_id, sample), grepl(class_type, class_id)) %>%
    select(class_id, meth.diff, qvalue)
}

# Define a list of class types to filter
class_types <- c("LTR", "LINE", "SINE")


# Filter LOUCY dataset
LD3_meth_ERV <- filter_samples(meth_loucy, "LD3", "LTR")
LD3_meth_LINE <- filter_samples(meth_loucy, "LD3", "LINE")
LD3_meth_SINE <- filter_samples(meth_loucy, "LD3", "SINE")

LG3_meth_ERV <- filter_samples(meth_loucy, "LG3", "LTR")
LG3_meth_LINE <- filter_samples(meth_loucy, "LG3", "LINE")
LG3_meth_SINE <- filter_samples(meth_loucy, "LG3", "SINE")

LG7_meth_ERV <- filter_samples(meth_loucy, "LG7", "LTR")
LG7_meth_LINE <- filter_samples(meth_loucy, "LG7", "LINE")
LG7_meth_SINE <- filter_samples(meth_loucy, "LG7", "SINE")

# Filter SUP-T1 dataset
SD3_meth_ERV <- filter_samples(meth_supt1, "SD3", "LTR")
SD3_meth_LINE <- filter_samples(meth_supt1, "SD3", "LINE")
SD3_meth_SINE <- filter_samples(meth_supt1, "SD3", "SINE")

SG3_meth_ERV <- filter_samples(meth_supt1, "SG3", "LTR")
SG3_meth_LINE <- filter_samples(meth_supt1, "SG3", "LINE")
SG3_meth_SINE <- filter_samples(meth_supt1, "SG3", "SINE")

SG7_meth_ERV <- filter_samples(meth_supt1, "SG7", "LTR")
SG7_meth_LINE <- filter_samples(meth_supt1, "SG7", "LINE")
SG7_meth_SINE <- filter_samples(meth_supt1, "SG7", "SINE")

#-------------------------------------------------------------------------------

              ## Extract data for Figure 5B----

#-------------------------------------------------------------------------------
# % of differentially methylated up and down based on total number
# of annotated TEs

#### LOUCY----
input_L <- list(LD3_meth_ERV, LG3_meth_ERV, LG7_meth_ERV, 
                LD3_meth_LINE, LG3_meth_LINE, LG7_meth_LINE,
                LD3_meth_SINE, LG3_meth_SINE, LG7_meth_SINE)

class_id <- rep(c("ERV", "LINE", "SINE"), each = 3)
treatment <- rep(c("LD3", "LG3", "LG7"), times= 3, each = 1)

summary_LOUCY <- data.frame(matrix(nrow=9, ncol=4))


for (i in 1:length(input_L)) {
  x <- input_L[[i]]
  colnames(x)[2] <- "meth.diff"
  hypo <- nrow(x[x$qvalue < 0.01 & x$meth.diff <= -25,])
  hyper <- nrow(x[x$qvalue < 0.01 & x$meth.diff >=  25,])
  
  df <- cbind(hypo, hyper, class_id[i], treatment[i])
  summary_LOUCY[i,] <- df
}

colnames(summary_LOUCY) <- c("hypo", "hyper", "class_id", "treatment")
summary_LOUCY$hypo <- as.numeric(summary_LOUCY$hypo)
summary_LOUCY$hyper <- as.numeric(summary_LOUCY$hyper)

# extract total number of TEs for each cell line
ERV_LOUCY <- nrow(LD3_meth_ERV) # number of annotated ERVs
LINE_LOUCY <- nrow(LD3_meth_LINE) # number of annotated LINEs
SINE_LOUCY <- nrow(LD3_meth_SINE) # number of annotated SINEs

# Calculate percent of hypo and hyper
summary_L_percent <- summary_LOUCY
summary_L_percent$percent_hypo <- 0
summary_L_percent$percent_hyper <- 0
summary_L_percent$percent_ns <- 0


for (i in 1:nrow(summary_L_percent)) {
  x <- summary_L_percent[i,]
    if(x$class_id=="ERV") {
    x$percent_hypo <- round(((x$hypo)/ERV_LOUCY)*100, digits = 2) # convert to positive value for hypo
    x$percent_hyper <- round(((x$hyper)/ERV_LOUCY)*100, digits = 2)
    x$percent_ns <- 100-x$percent_hypo-x$percent_hyper
  } else if(x$class_id=="LINE") {
    x$percent_hypo <- round((x$hypo/LINE_LOUCY)*100, digits = 2)
    x$percent_hyper <- round((x$hyper/LINE_LOUCY)*100, digits = 2)
    x$percent_ns <- 100-x$percent_hypo-x$percent_hyper
      } else if(x$class_id=="SINE") {
    x$percent_hypo <- round((x$hypo/SINE_LOUCY)*100, digits = 2)
    x$percent_hyper <- round((x$hyper/SINE_LOUCY)*100, digits = 2)
    x$percent_ns <- 100-x$percent_hypo-x$percent_hyper
    
  } 
  summary_L_percent[i,] <- x
  rm(x)
}


#### SUPT1----
input_S <- list(SD3_meth_ERV, SG3_meth_ERV, SG7_meth_ERV, 
                SD3_meth_LINE, SG3_meth_LINE, SG7_meth_LINE,
                SD3_meth_SINE, SG3_meth_SINE, SG7_meth_SINE)

class_id <- rep(c("ERV", "LINE", "SINE"), each = 3)
treatment <- rep(c("SD3", "SG3", "SG7"), times= 3, each = 1)

summary_SUPT1 <- data.frame(matrix(nrow=9, ncol=4))

for (i in 1:length(input_S)) {
  x <- input_S[[i]]
  colnames(x)[2] <- "meth.diff"
  hypo <- nrow(x[x$qvalue < 0.01 & x$meth.diff <= -25,])
  hyper <- nrow(x[x$qvalue < 0.01 & x$meth.diff >= 25,])
 
  df <- cbind(hypo, hyper, class_id[i], treatment[i])
  summary_SUPT1[i,] <- df
}

colnames(summary_SUPT1) <- c("hypo", "hyper", "class_id", "treatment")
summary_SUPT1$hypo <- as.numeric(summary_SUPT1$hypo)
summary_SUPT1$hyper <- as.numeric(summary_SUPT1$hyper)

ERV_SUPT1 <- nrow(SD3_meth_ERV)
LINE_SUPT1 <- nrow(SD3_meth_LINE)
SINE_SUPT1 <- nrow(SD3_meth_SINE)

summary_S_percent <- summary_SUPT1
summary_S_percent$percent_hypo <- 0
summary_S_percent$percent_hyper <- 0
summary_S_percent$percent_ns <- 0


for (i in 1:nrow(summary_S_percent)) {
  x <- summary_S_percent[i,]
  if(x$class_id=="ERV") {
    x$percent_hypo <- round((x$hypo/ERV_LOUCY)*100, digits = 2)
    x$percent_hyper <- round(((x$hyper)/ERV_LOUCY)*100, digits = 2)
    x$percent_ns <- 100-x$percent_hypo-x$percent_hyper
    
  } else if(x$class_id=="LINE") {
    x$percent_hypo <- round((x$hypo/LINE_LOUCY)*100, digits = 2)
    x$percent_hyper <- round((x$hyper/LINE_LOUCY)*100, digits = 2)
    x$percent_ns <- 100-x$percent_hypo-x$percent_hyper
    
  } else if(x$class_id=="SINE") {
    x$percent_hypo <- round((x$hypo/SINE_LOUCY)*100, digits = 2)
    x$percent_hyper <- round((x$hyper/SINE_LOUCY)*100, digits = 2)
    x$percent_ns <- 100-x$percent_hypo-x$percent_hyper
    
  } 
  summary_S_percent[i,] <- x
  rm(x)
}

# Figure 5B was done in Graphpad using summary_L_percent and summary_S_percent, columns percent_hypo, percent_hyper and percent_ns







