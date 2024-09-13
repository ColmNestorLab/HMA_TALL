################################################################################
############## MethylKit: DNA METHYLATION AT PROMOTERS PROCESSING ##############
################################################################################

# Maike Bensberg (maike.bensberg@liu.se)
# 2024-08-13


### analysing DNA methylation at promoters
# summarize DNA methylation at each promoter
# some genes will have several TSS and several promoters
# identify differential methylation at promoters

# exported data will become Supplementary Table 8
# exported data used in Figure 4c


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
library(rtracklayer)



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



##### Summarizing DNA methylation at promoters #####
# get average DNA methylation across each promoter (methylKit)
# done by cell line (different cell lines normalized/analyzed separately)


### load gene annotation file
# only curated transcripts (NM_* and NR_*)
# promoters defined as 1000bp up- and 500bp downstream from TSS

gene_annot <- readTranscriptFeatures("R_imports/transcripts_curated_hg38_ucsc_20240814.bed",
                                     up.flank = 1000, down.flank = 500)


### summarize counts
counts_promoters_L <- regionCounts(EMseq_object_L, gene_annot$promoters)
counts_promoters_S <- regionCounts(EMseq_object_S, gene_annot$promoters)


### normalize coverage between samples
# needed for comparative analysis

promoters_L_norm <- normalizeCoverage(counts_promoters_L, method = "median")
promoters_S_norm <- normalizeCoverage(counts_promoters_S, method = "median")


### merge files
# only keep promoters covered in all samples of a given cell line
# information per CpG (merge strands)

promoters_L_norm_merged <- methylKit::unite(promoters_L_norm, destrand = TRUE) 
promoters_S_norm_merged <- methylKit::unite(promoters_S_norm, destrand = TRUE)




##### Identifying differentially methylated promoters #####

# get sample IDs
L_IDs <- getSampleID(promoters_L_norm_merged)
S_IDs <- getSampleID(promoters_S_norm_merged)


### compare each treatment to the control
# Fishers exact test to calculate p-values (test = "fast.fisher")
# q-values: p-values adjusted for multiple testing with Benjamini Hochberg (adjust = "BH")
# use 8 cores (mc.cores = 8)

for (s in L_IDs[! L_IDs == "Lctrl7"]) {
  assign(paste0("DMR_promoters_",s),
         calculateDiffMeth(reorganize(promoters_L_norm_merged, sample.ids = c("Lctrl7", s), treatment = c(0,1)),
                           test = "fast.fisher", slim = FALSE, adjust = "BH", mc.cores = 8))
}
for (s in S_IDs[! S_IDs == "Sctrl7"]) {
  assign(paste0("DMR_promoters_",s),
         calculateDiffMeth(reorganize(promoters_S_norm_merged, sample.ids = c("Sctrl7", s), treatment = c(0,1)),
                           test = "fast.fisher", slim = FALSE, adjust = "BH", mc.cores = 8))
}

remove(s)
remove(L_IDs)
remove(S_IDs)




##### Annotate differentially methylated promoters #####

# make dataframe
# annotate promoters with gene_id and transcript_id


### import gene annotation file
# use same annotation file as used for RNA-seq analysis:
# GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf

# import annotation files
gtf_promoters <- rtracklayer::import("R_imports/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf")
df_promoters <- data.frame(gtf_promoters)

# rename first column to "chr" to be named similar to dataframe with promoter coordinates from analysis
colnames(df_promoters) [1] <- "chr"


### convert GRanges file from DMR analysis to dataframe

# get all files for DMR analysis of promoters (separate files for each teatment)
DMR_promoters_files <- ls(pattern = "\\bDMR_promoters_...\\b")

# make into dataframe
for (f in DMR_promoters_files) {
  assign(paste0("df_", f),
         getData(get(f)))
}

remove(DMR_promoters_files)
remove(list = ls(pattern = "\\bDMR_promoters_...\\b"))


### add transcript_id
# transcript_ID based on gene annotation file used above (gene_annot)

# get all dataframes for DMR promoter analysis
DMR_promoters_files <- ls(pattern = "\\bdf_DMR_promoters_...\\b")

### make new dataframe with TSS location associated with promoter
# define function
prom_to_TSS <- function(df) {
  df[which(df$strand == "+"), "start"] <-
    df[which(df$strand == "+"), "start"] + 1000
  df[which(df$strand == "+"), "end"] <-
    df[which(df$strand == "+"), "end"] - 500
  df[which(df$strand == "-"), "start"] <-
    df[which(df$strand == "-"), "start"] + 500
  df[which(df$strand == "-"), "end"] <-
    df[which(df$strand == "-"), "end"] - 1000
  return(df)
}

# make new dataframe with TSS position 
for (f in DMR_promoters_files) {
  assign(paste0(f,"_TSS"), prom_to_TSS(get(f)))
}

# assign transcript ID
for (f in DMR_promoters_files) {
  assign(f, cbind(get(f), "transcript_id" =
                    getAssociationWithTSS(annotateWithGeneParts(
                      target = as(get(paste0(f, "_TSS")), "GRanges"), feature = gene_annot))[,"feature.name"]))
}

remove(f)
remove(DMR_promoters_files)
remove(list = ls(pattern = "\\bdf_DMR_promoters_..._TSS\\b"))


### add gene ID

# annotation file imported above (gtf_promoters/df_promoters)
# same as used for RNA-seq

# get all dataframes for DMR promoter analysis
DMR_promoters_files <- ls(pattern = "\\bdf_DMR_promoters_...\\b")

# assign gene ID
for (f in DMR_promoters_files) {
  assign(f,
         left_join(get(f), df_promoters[which(df_promoters$type == "transcript"),c(10,11)]))
}

remove(DMR_promoters_files)
remove(f)



##### add methylation for treated and control samples #####

# get sample IDs
L_IDs <- getSampleID(promoters_L_norm_merged)
S_IDs <- getSampleID(promoters_S_norm_merged)

# find index for each sample in vector with samples
L_index <- which(L_IDs != "Lctrl7")
S_index <- which(S_IDs != "Sctrl7")

# add DNA methylation for the treated sample
i = 1
for (s in L_IDs[! L_IDs == "Lctrl7"]) {
  assign(paste0("df_DMR_promoters_", s),
        cbind(get(paste0("df_DMR_promoters_", s)), "methylation_treated" =
              getData(promoters_L_norm_merged)[,paste0("numCs",L_index[i])]/
                getData(promoters_L_norm_merged)[,paste0("coverage",L_index[i])]))
  assign(paste0("df_DMR_promoters_", s),
         cbind(get(paste0("df_DMR_promoters_", s)), "methylation_control" =
                 getData(promoters_L_norm_merged)[,paste0("numCs",which(L_IDs == "Lctrl7"))]/
                 getData(promoters_L_norm_merged)[,paste0("coverage",which(L_IDs == "Lctrl7"))]))
  i = i+1
}

i = 1
for (s in S_IDs[! S_IDs == "Sctrl7"]) {
  assign(paste0("df_DMR_promoters_", s),
         cbind(get(paste0("df_DMR_promoters_", s)), "methylation_treated" =
                 getData(promoters_S_norm_merged)[,paste0("numCs",S_index[i])]/
                 getData(promoters_S_norm_merged)[,paste0("coverage",S_index[i])]))
  assign(paste0("df_DMR_promoters_", s),
         cbind(get(paste0("df_DMR_promoters_", s)), "methylation_control" =
                 getData(promoters_S_norm_merged)[,paste0("numCs",which(S_IDs == "Sctrl7"))]/
                 getData(promoters_S_norm_merged)[,paste0("coverage",which(S_IDs == "Sctrl7"))]))
  i = i+1
}

remove(i)
remove(s)
remove(L_IDs)
remove(S_IDs)
remove(L_index)
remove(S_index)



##### combine data and export as dataframe #####

# combine DMR data for all treatments of each cell line
df_DMR_promoters_L_Ntranscripts <- do.call("rbind", list(df_DMR_promoters_LD3, df_DMR_promoters_LG3, df_DMR_promoters_LG7))
df_DMR_promoters_S_Ntranscripts <- do.call("rbind", list(df_DMR_promoters_SD3, df_DMR_promoters_SG3, df_DMR_promoters_SG7))

# add sample
df_DMR_promoters_L_Ntranscripts$sample <- c(rep("LD3", nrow(df_DMR_promoters_LD3)),
                               rep("LG3", nrow(df_DMR_promoters_LG3)),
                               rep("LG7", nrow(df_DMR_promoters_LG7)))
df_DMR_promoters_S_Ntranscripts$sample <- c(rep("SD3", nrow(df_DMR_promoters_SD3)),
                               rep("SG3", nrow(df_DMR_promoters_SG3)),
                               rep("SG7", nrow(df_DMR_promoters_SG7)))

# combine cell lines
df_DMR_promoters <- rbind(df_DMR_promoters_L_Ntranscripts, df_DMR_promoters_S_Ntranscripts)

remove(df_DMR_promoters_L_Ntranscripts)
remove(df_DMR_promoters_S_Ntranscripts)
remove(list = ls(pattern = "\\bdf_DMR_promoters_...\\b"))


### export

# write as tab separated file
write.table(df_DMR_promoters, "R_exports/tables/Supplementary_Table_8.csv",
            sep = ",", quote = FALSE, col.names = TRUE, row.names = FALSE)