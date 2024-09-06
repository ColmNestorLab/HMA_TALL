# Title: Figure 5A

# Author: Sandra Hellberg
# Email: sandra.hellberg@liu.se

# Date: 2024-02-08
# Last modified: 2024-08-14

# Description: 
rm(list=ls()) # remove all entries in the global environment 
set.seed(102) # seed for reproducibility

#-------------------------------------------------------------------------------

              ## Set directory structure and load packages---- 

#-------------------------------------------------------------------------------
setwd("Set your working directory")

lib_path <- .libPaths("C:/Rpackages/") # specifing the package library
#update.packages(lib.loc = lib_path, ask = FALSE, checkBuilt = TRUE, dependencies = TRUE)
#BiocManager::install(lib.loc = lib_path, update = TRUE, ask = FALSE)

pack_R <- c("magrittr", "stringr", "circlize", "dplyr", "ComplexHeatmap") # libraries to load

for (i in 1:length(pack_R)) {
  library(pack_R[i], character.only = TRUE)
}

#-------------------------------------------------------------------------------

                  ## Read in methylation data----

#-------------------------------------------------------------------------------
# Read the data from the CSV file
meth_loucy <- read.delim("Supplementary_Table_14.csv", sep = ",") # LOUCY
meth_supt1 <- read.delim("Supplementary_Table_13.csv", sep = ",") # SUP-T1

# Filter rows that contain specific sample patterns
TE_LD3 <- meth_loucy %>% filter(grepl("LD3", sample))
TE_LG3 <- meth_loucy %>% filter(grepl("LG3", sample))
TE_LG7 <- meth_loucy %>% filter(grepl("LG7", sample))

TE_SD3 <- meth_supt1 %>% filter(grepl("SD3", sample))
TE_SG3 <- meth_supt1 %>% filter(grepl("SG3", sample))
TE_SG7 <- meth_supt1 %>% filter(grepl("SG7", sample))



#-------------------------------------------------------------------------------

          ## Combine TEs per treatment and bin into 500 kB bins----

#-------------------------------------------------------------------------------
#### LOUCY----
input_L <- list(TE_LD3, TE_LG3, TE_LG7)
chr <- unique(TE_LD3$chr) # create a vector with all chromosomes
bin_size <- 500000
result_df <- list(LD3 = data.frame(), LG3 = data.frame(), LG7 = data.frame()) 

for(i in 1:length(input_L)) {
  x <- input_L[[i]]
  
  for (n in 1:length(chr)) {
    x_subset <- x[x$chr %in% chr[n],]
    
    # Calculate bin indices manually
    x_subset$bin <- floor((x_subset$start - 1) / bin_size) + 1
    
    x_subset$bin_start <- (x_subset$bin - 1) * bin_size + 1
    x_subset$bin_end <- x_subset$bin * bin_size
    
    x_subset_summary <- x_subset %>%
      group_by(bin) %>%
      summarise(
        chr = first(chr),
        min = min(start),
        max = max(end),
        mean_treated = mean(methylation_treated),
        mean_control = mean(methylation_control),
        bin_start = first(bin_start),
        bin_end = first(bin_end)
      ) %>%
      ungroup()
    
    # Append the results to the data frame for this iteration
    result_df[[i]] <- bind_rows(result_df[[i]], x_subset_summary)
  }
}
    
LD3_bin <- result_df[["LD3"]]    
LG3_bin <- result_df[["LG3"]]   
LG7_bin <- result_df[["LG7"]]      
    
#### SUP-T1----
input_S <- list(TE_SD3, TE_SG3, TE_SG7)
chr <- unique(TE_SD3$chr) # create a vector with all chromosomes
bin_size <- 500000
result_df <- list(SD3 = data.frame(), SG3 = data.frame(), SG7 = data.frame()) 


for(i in 1:length(input_S)) {
  x <- input_S[[i]]
  
  for (n in 1:length(chr)) {
    x_subset <- x[x$chr %in% chr[n],]
    # Calculate bin indices manually
    x_subset$bin <- floor((x_subset$start - 1) / bin_size) + 1
    
    x_subset$bin_start <- (x_subset$bin - 1) * bin_size + 1
    x_subset$bin_end <- x_subset$bin * bin_size
    
    x_subset_summary <- x_subset %>%
      group_by(bin) %>%
      summarise(
        chr = first(chr),
        min = min(start),
        max = max(end),
        mean_treated = mean(methylation_treated),
        mean_control = mean(methylation_control),
        bin_start = first(bin_start),
        bin_end = first(bin_end)
      ) %>%
      ungroup()
    
    # Append the results to the data frame for this iteration
    result_df[[i]] <- bind_rows(result_df[[i]], x_subset_summary)
  }
}

SD3_bin <- result_df[["SD3"]]    
SG3_bin <- result_df[["SG3"]]   
SG7_bin <- result_df[["SG7"]]  

#-------------------------------------------------------------------------------    

                            ## Circos----

#-------------------------------------------------------------------------------    
# defining the colors for plotting
col_fun1 <- colorRamp2(c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), c("#071630" ,"#1984c5",  "#a7d5ed", "#66C2A5","#ABDDA4" ,"#E6F598","#FFFFBF","#FDAE61","#F46D43", "#D53E4F", "#9E0142"))

### LOUCY----
split <-  LG7_bin$chr 
split <- factor(split, levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
                                  "chr9", "chr10", "chr11", "chr12", "chr13", "chr14","chr15", "chr16",
                                  "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")) # define levels for the plotting

circos.clear() # clear plot
circos.par(gap.after = c(rep(2,22),30))

circ1 <- as.matrix(LG7_bin$mean_treated)
circos.heatmap(circ1 , split = split, col = col_fun1, 
               show.sector.labels = T, cluster = FALSE,
               track.height = 0.1)

circ2 <- as.matrix(LG3_bin$mean_treated)
circos.heatmap(circ2, col = col_fun1, track.height = 0.1)

circ3 <- as.matrix(LD3_bin$mean_treated)
circos.heatmap(circ3, col = col_fun1, track.height = 0.1)

circ4 <- as.matrix(LG3_bin$mean_control)
circos.heatmap(circ4, col = col_fun1, track.height = 0.1)    
    
lgd <- Legend(title = "LOUCY", col_fun = col_fun1)
grid.draw(lgd)   
  
# Export as pdf 

### SUP-T1----

circ5 <- as.matrix(SG7_bin$mean_treated)

split2 <-  SG7_bin$chr
split2 <- factor(split2, levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
                                    "chr9", "chr10", "chr11", "chr12", "chr13", "chr14","chr15", "chr16",
                                    "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"))

circos.clear()
circos.par(gap.after = c(rep(2,22),30))
circos.heatmap(circ5, split = split2, col = col_fun1, 
               show.sector.labels = TRUE, cluster = FALSE,
               track.height = 0.1)

circ6 <- as.matrix(SG3_bin$mean_treated)
circos.heatmap(circ6, col = col_fun1, track.height = 0.1)

circ7 <- as.matrix(SD3_bin$mean_treated)
circos.heatmap(circ7, col = col_fun1, track.height = 0.1)

circ8 <- as.matrix(SG3_bin$mean_control)
circos.heatmap(circ8, col = col_fun1, track.height = 0.1)

lgd <- Legend(title = "SUP-T1", col_fun = col_fun1)
grid.draw(lgd)   


# Export as pdf 


