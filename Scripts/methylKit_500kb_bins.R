################################################################################
####################### MethylKit: 500 KB BINS PROCESSING ######################
################################################################################

# Maike Bensberg (maike.bensberg@liu.se)
# 2024-08-12


### process whole genome DNA methylation data
# summarize DNA methylation across 500 kb bins
# export data table to be used in making figures

# exported data used in Fig. 4a


### using MethylKit 
# https://bioconductor.org/packages/release/bioc/html/methylKit.html
# https://github.com/al2na/methylKit
# https://bioconductor.org/packages/release/bioc/manuals/methylKit/man/methylKit.pdf
# https://doi.org/10.1186/gb-2012-13-10-r87
# Akalin, A., Kormaksson, M., Li, S. et al. methylKit: a comprehensive R package for the analysis of genome-wide DNA methylation profiles. Genome Biol 13, R87 (2012). 
# Analysis based on vignette: https://bioconductor.org/packages/release/bioc/vignettes/methylKit/inst/doc/methylKit.html



##### initiate R environment #####
# only needs to be done once
#install.packages("renv")
#renv::init()



##### loading packages #####

library(methylKit)



##### set working directory and check package versions #####

# set working directory
setwd("/mnt/wwn-0x5000cca28fd30c7d-part1/HMA_TALL/EMseq/LOUCY_SUPT1_HMAtreated_12_2023/methylKit")
# print working directory
getwd()

# print version of R and attached packages
sessionInfo()



##### reading in DNA methylation call files #####
# make sure files are compatible with methylKit
# raw fastq files processed by script EMseq_CpGm_analysis.sh

EMseq_files <- as.list(list.files(path = "R_imports/methylKit_compatible",
                                  full.names = TRUE))
names <- list("SG7","LD3","LG3","Lctrl7","LG7","SD3","SG3","Sctrl7")

# minimum coverage of 5 reads for each cytosine
EMseq_object <- methRead(location = EMseq_files,
                         sample.id = names,
                         assembly = "hg38",
                         pipeline = "bismark",
                         context = "CpG",
                         resolution = "base",
                         treatment = c(1,1,1,0,1,1,1,0),
                         mincov = 5)

remove(names)



##### separate information by cell line #####

# subset
EMseq_object_L <- reorganize(EMseq_object, sample.ids = c("LD3","LG3","Lctrl7","LG7"),
                             treatment = c(1,1,0,1))
EMseq_object_S <- reorganize(EMseq_object, sample.ids = c("SG7","SD3","SG3","Sctrl7"),
                             treatment = c(1,1,1,0))



##### separate genome into bins #####


### summarize information in bins of 500 kb
# window size: 500,000 bases
# non-overlapping bins (step size = window size = 500,000)
# use 8 cores (mc.cores = 8)
EMseq_bins500kb_L <- tileMethylCounts(EMseq_object_L, mc.cores = 8,
                                      win.size = 500000, step.size = 500000)
EMseq_bins500kb_S <- tileMethylCounts(EMseq_object_S, mc.cores = 8,
                                      win.size = 500000, step.size = 500000)


### normalize coverage between samples
# needed for comparative analysis
EMseq_bins500kb_L_norm <- normalizeCoverage(EMseq_bins500kb_L, method = "median")
EMseq_bins500kb_S_norm <- normalizeCoverage(EMseq_bins500kb_S, method = "median")


### merge files
# only pick bases covered in all samples
# information per CpG (merge strands)
EMseq_bins500kb_L_merged <- methylKit::unite(object = EMseq_bins500kb_L_norm, destrand = TRUE)
EMseq_bins500kb_S_merged <- methylKit::unite(object = EMseq_bins500kb_S_norm, destrand = TRUE)



##### convert information into data frame format #####

# make dataframe
df_temp_L <- getData(EMseq_bins500kb_L_merged)
df_temp_S <- getData(EMseq_bins500kb_S_merged)


### calculate percent DNA methylation for each bin and sample

# regions
df_EMseq_bins500kb_L_merged <- df_temp_L[,1:3]
df_EMseq_bins500kb_S_merged <- df_temp_S[,1:3]

# get sample IDs
L_IDs <- getSampleID(EMseq_bins500kb_L_merged)
S_IDs <- getSampleID(EMseq_bins500kb_S_merged)

# add DNA methylation (as number between 0 and 1)
# sample names as column names
for (i in c(1:length(L_IDs))) {
  df_EMseq_bins500kb_L_merged <- cbind(df_EMseq_bins500kb_L_merged,
                                       df_temp_L[,paste0("numCs", i)]/df_temp_L[,paste0("coverage", i)])
  colnames(df_EMseq_bins500kb_L_merged)[ncol(df_EMseq_bins500kb_L_merged)] <- L_IDs[i]
}
for (i in c(1:length(S_IDs))) {
  df_EMseq_bins500kb_S_merged <- cbind(df_EMseq_bins500kb_S_merged,
                                       df_temp_S[,paste0("numCs", i)]/df_temp_S[,paste0("coverage", i)])
  colnames(df_EMseq_bins500kb_S_merged)[ncol(df_EMseq_bins500kb_S_merged)] <- S_IDs[i]
}

remove(i)
remove(L_IDs)
remove(S_IDs)
remove(df_temp_L)
remove(df_temp_S)



##### export data frame #####

# write as tab separated file
write.table(df_EMseq_bins500kb_L_merged, "R_exports/tables/df_EMseq_bins500kb_L_merged.txt",
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
write.table(df_EMseq_bins500kb_S_merged, "R_exports/tables/df_EMseq_bins500kb_S_merged.txt",
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)