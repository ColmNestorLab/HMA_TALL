# Supplementary Tables 15 & 16

# Author: Sandra Hellberg
# Email: sandra.hellberg@liu.se

# Date: 2024-02-02
# Last modified: 2024-08-14

# Description: 

rm(list=ls()) # remove all entries in the global environment 
set.seed(524) # seed for reproducibility

#-------------------------------------------------------------------------------

        ## Set directory structure and load packages---- 

#-------------------------------------------------------------------------------
setwd("working directory")
meth_dir <- "input directory"

lib_path <- .libPaths("C:/Rpackages/") # specifing the package library
#update.packages(lib.loc = lib_path, ask = FALSE, checkBuilt = TRUE, dependencies = TRUE)
#BiocManager::install(lib.loc = lib_path, update = TRUE, ask = FALSE)

pack_R <- c("xlsx") # libraries to load

for (i in 1:length(pack_R)) {
  library(pack_R[i], character.only = TRUE)
}


#-------------------------------------------------------------------------------

                      ## Read in and format data---- 

#-------------------------------------------------------------------------------

# LOUCY
meth_L <- list.files(meth_dir, pattern = "L.*_meth")

methylation_L <- data.frame()  # Initialize if it doesn't exist

for(i in 1:length(meth_L)) {
  x <- readRDS(paste0(meth_dir,meth_L[i]))
  methylation_L <- rbind(methylation_L, x)
  rm(x)
}

methylation_L_adjusted <- methylation_L %>% # add start +1 to make it comparable to the gft used to assign locations for the RNA-seq
  mutate(start = start + 1)


# SUPT1
meth_S <- list.files(meth_dir, pattern = "S.*_meth")

methylation_S <- data.frame()  # Initialize if it doesn't exist

for(i in 1:length(meth_S)) {
  x <- readRDS(paste0(meth_dir,meth_S[i]))
  methylation_S <- rbind(methylation_S, x)
  rm(x)
}

methylation_S_adjusted <- methylation_S %>% # add start +1 to make it comparable to the gft used to assign locations for the RNA-seq
  mutate(start = start + 1)


#-------------------------------------------------------------------------------

                    ## Export supplemental table---- 

#-------------------------------------------------------------------------------

write.csv(methylation_L_adjusted, "Supplementary_Table_15.csv", quote = F, row.names = F)
write.csv(methylation_S_adjusted, "Supplementary_Table_16.csv", quote = F, row.names = F)
