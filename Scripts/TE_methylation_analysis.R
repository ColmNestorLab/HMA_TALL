################################################################################
######## MethylKit: DNA METHYLATION AT TRANSPOSABLE ELEMENTS PROCESSING ########
################################################################################

# Maike Bensberg (maike.bensberg@liu.se)
# 2024-08-12


### analysing DNA methylation at Transposable elements
# combined analysis for LINEs, SINEs and ERVs
# summarize DNA methylation at each TE
# identify differential methylation at TEs

# exported data used in Fig. 5 and Extended Data Fig. 9


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
library(genomation)
library(dplyr)



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


### separate information by cell line

EMseq_object_L <- reorganize(EMseq_object, sample.ids = c("LD3","LG3","Lctrl7","LG7"),
                             treatment = c(1,1,0,1))
EMseq_object_S <- reorganize(EMseq_object, sample.ids = c("SG7","SD3","SG3","Sctrl7"),
                             treatment = c(1,1,1,0))



##### Summarizing DNA methylation at TEs #####
# get average DNA methylation across each TE (methylKit)
# export files for DNA methylation of each sample for TEs


### load annotation file for transposable elements
# provided by Sandra H.
# same annotation used for analysis of RNA-seq

TEs_annot <- readFeatureFlank("R_imports/hg38_rmsk_TE_clean.bed")


### summarize counts
counts_TEs_L <- regionCounts(EMseq_object_L, TEs_annot$features)
counts_TEs_S <- regionCounts(EMseq_object_S, TEs_annot$features)


### normalize coverage between samples
# needed for comparative analysis

TEs_L_norm <- normalizeCoverage(counts_TEs_L, method = "median")
TEs_S_norm <- normalizeCoverage(counts_TEs_S, method = "median")


### merge files
# only keep TEs covered in all samples
# information per CpG (merge strands)

TEs_L_norm_merged <- methylKit::unite(TEs_L_norm, destrand = TRUE) 
TEs_S_norm_merged <- methylKit::unite(TEs_S_norm, destrand = TRUE)



##### annotate and export DNA methylation at TEs #####


### make dataframe
df_temp_L <- getData(TEs_L_norm_merged)
df_temp_S <- getData(TEs_S_norm_merged)

# get TE coordinates
df_TEs_DNAmethylation_L <- df_temp_L[,1:3]
df_TEs_DNAmethylation_S <- df_temp_S[,1:3]


### annotate TEs with gene_id, transcript_id, family_id and class_id

# import annotation files
gtf_TEs <- rtracklayer::import("R_imports/hg38_rmsk_TE_clean.gtf")
df_TEs <- data.frame(gtf_TEs)

# rename first column to "chr" to be named similar to df with TE coordinates from analysis
colnames(df_TEs) [1] <- "chr"

# convert chr location from 1-based, closed (gft format) to 0-based, half-open (bed format, used for annotation)
df_TEs$start <- df_TEs$start-1

# add information to dataframes for each cell line
df_TEs_DNAmethylation_L <- left_join(df_TEs_DNAmethylation_L, df_TEs[,c(1,2,3,10,11,12,13)])
df_TEs_DNAmethylation_S <- left_join(df_TEs_DNAmethylation_S, df_TEs[,c(1,2,3,10,11,12,13)])


### add DNA methylation for each sample to dataframes
# DNA methylation value between 0 and 1

# get sample IDs
L_IDs <- getSampleID(TEs_L_norm_merged)
S_IDs <- getSampleID(TEs_S_norm_merged)

# add DNA methylation for each sample
for (i in c(1:length(L_IDs))) {
  df_TEs_DNAmethylation_L <- cbind(df_TEs_DNAmethylation_L,
                                       df_temp_L[,paste0("numCs", i)]/df_temp_L[,paste0("coverage", i)])
  colnames(df_TEs_DNAmethylation_L)[ncol(df_TEs_DNAmethylation_L)] <- L_IDs[i]
}
for (i in c(1:length(S_IDs))) {
  df_TEs_DNAmethylation_S <- cbind(df_TEs_DNAmethylation_S,
                                   df_temp_S[,paste0("numCs", i)]/df_temp_S[,paste0("coverage", i)])
  colnames(df_TEs_DNAmethylation_S)[ncol(df_TEs_DNAmethylation_S)] <- S_IDs[i]
}

remove(i)
remove(L_IDs)
remove(S_IDs)


### export DNA methylation

# write as tab separated file
write.table(df_TEs_DNAmethylation_L, "R_exports/tables/df_TEs_DNAmethylation_L.txt",
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
write.table(df_TEs_DNAmethylation_S, "R_exports/tables/df_TEs_DNAmethylation_S.txt",
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)



##### Identifying differentially methylated TEs #####

# get sample IDs
L_IDs <- getSampleID(TEs_L_norm_merged)
S_IDs <- getSampleID(TEs_S_norm_merged)


### compare each treatment to the control
# Fishers exact test to calculate p-values (test = "fast.fisher")
# q-values: p-values adjusted for multiple testing with Benjamini Hochberg (adjust = "BH")
# use 8 cores (mc.cores = 8)

for (s in L_IDs[! L_IDs == "Lctrl7"]) {
  assign(paste0("DMR_TEs_",s),
         calculateDiffMeth(reorganize(TEs_L_norm_merged, sample.ids = c("Lctrl7", s), treatment = c(0,1)),
                           test = "fast.fisher", slim = FALSE, adjust = "BH", mc.cores = 8))
}
for (s in S_IDs[! S_IDs == "Sctrl7"]) {
  assign(paste0("DMR_TEs_",s),
         calculateDiffMeth(reorganize(TEs_S_norm_merged, sample.ids = c("Sctrl7", s), treatment = c(0,1)),
                           test = "fast.fisher", slim = FALSE, adjust = "BH", mc.cores = 8))
}

remove(s)
remove(L_IDs)
remove(S_IDs)



##### Annotate differentially methylated TEs #####
# make dataframe
# annotate TEs with gene_id, transcript_id, family_id and class_id
# use same annotation file as above: df_TEs

# get all files for DMR analysis of TEs (separate files for each treatment)
DMR_TE_files <- ls(pattern = "\\bDMR_TEs_...\\b")

# annotate files
for (f in DMR_TE_files) {
  assign(paste0("df_", f),
         getData(get(f)))
  assign(paste0("df_", f),
         left_join(get(paste0("df_", f)), df_TEs[,c(1,2,3,10,11,12,13)]))
}

# add methylation for treated and control samples
df_DMR_TEs_LD3$methylation_treated <- df_TEs_DNAmethylation_L$LD3
df_DMR_TEs_LG3$methylation_treated <- df_TEs_DNAmethylation_L$LG3
df_DMR_TEs_LG7$methylation_treated <- df_TEs_DNAmethylation_L$LG7
for (f in ls(pattern = "^DMR_TEs_L")) {
  assign(paste0("df_",f),
         cbind(get(paste0("df_",f)), "methylation_control" = df_TEs_DNAmethylation_L$Lctrl7))
}
df_DMR_TEs_SD3$methylation_treated <- df_TEs_DNAmethylation_S$SD3
df_DMR_TEs_SG3$methylation_treated <- df_TEs_DNAmethylation_S$SG3
df_DMR_TEs_SG7$methylation_treated <- df_TEs_DNAmethylation_S$SG7
for (f in ls(pattern = "^DMR_TEs_S")) {
  assign(paste0("df_",f),
         cbind(get(paste0("df_",f)), "methylation_control" = df_TEs_DNAmethylation_S$Sctrl7))
}

remove(f)
remove(DMR_TE_files)



##### combine data by cell line and export as dataframe #####

# combine DMR data for all treatments of each cell line
df_DMR_TEs_L <- do.call("rbind", list(df_DMR_TEs_LD3, df_DMR_TEs_LG3, df_DMR_TEs_LG7))
df_DMR_TEs_S <- do.call("rbind", list(df_DMR_TEs_SD3, df_DMR_TEs_SG3, df_DMR_TEs_SG7))

# add sample
df_DMR_TEs_L$sample <- c(rep("LD3", nrow(df_DMR_TEs_LD3)),
                         rep("LG3", nrow(df_DMR_TEs_LG3)),
                         rep("LG7", nrow(df_DMR_TEs_LG7)))
df_DMR_TEs_S$sample <- c(rep("SD3", nrow(df_DMR_TEs_SD3)),
                         rep("SG3", nrow(df_DMR_TEs_SG3)),
                         rep("SG7", nrow(df_DMR_TEs_SG7)))

# combine both cell lines
df_DMR_TEs <- rbind(df_DMR_TEs_L, df_DMR_TEs_S)


### export

# write as tab separated file
write.table(df_DMR_TEs, "R_exports/tables/df_DMR_TEs.txt",
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)