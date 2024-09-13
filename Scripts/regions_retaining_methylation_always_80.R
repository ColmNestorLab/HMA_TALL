################################################################################
################## IDENTIFY REGIONS RETAINING DNA METHYLATION ################## 
################################################################################

# Maike Bensberg (maike.bensberg@liu.se)
# 2024-08-14


### identify regions retaining DNA methylation
# In what regions is DNA methylation essential?
# retention of DNA methylation defined as losing less DNA methylation than other regions
# only interested in regions that are > 80% methylated in untreated and treated cells
# only look at samples retaining DNA methylation when treated with GSK for 7 days


### using MethylKit for dividing genome into bins
# https://bioconductor.org/packages/release/bioc/html/methylKit.html
# https://github.com/al2na/methylKit
# https://bioconductor.org/packages/release/bioc/manuals/methylKit/man/methylKit.pdf
# https://doi.org/10.1186/gb-2012-13-10-r87
# Akalin, A., Kormaksson, M., Li, S. et al. methylKit: a comprehensive R package for the analysis of genome-wide DNA methylation profiles. Genome Biol 13, R87 (2012). 
# Analysis based on vignette: https://bioconductor.org/packages/release/bioc/vignettes/methylKit/inst/doc/methylKit.html


### exported data used in Figure 4b and Extended Data Fig. 5a-b
# exported gene IDs will be combined to become supplementary table 9



##### initiate R environment #####
# only needs to be done once
#install.packages("renv")
#renv::init()



##### loading packages #####

library(methylKit)
library(genomation)
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
# only untreated samples or treated for 7 days with GSK

EMseq_files <- as.list(list.files(path = "R_imports/methylKit_compatible",
                                  full.names = TRUE, pattern = "_7d_"))
names <- list("SG7","Lctrl7","LG7","Sctrl7")

# minimum coverage of 5 reads for each cytosine
EMseq_object <- methRead(location = EMseq_files,
                         sample.id = names,
                         assembly = "hg38",
                         pipeline = "bismark",
                         context = "CpG",
                         resolution = "base",
                         treatment = c(1,0,1,0),
                         mincov = 5)

remove(names)



##### separate information by cell line #####

# subset
EMseq_object_L <- reorganize(EMseq_object, sample.ids = c("Lctrl7","LG7"),
                             treatment = c(0,1))
EMseq_object_S <- reorganize(EMseq_object, sample.ids = c("SG7","Sctrl7"),
                             treatment = c(1,0))



##### separate genome into bins #####


### summarize information in bins of 1 kb
# window size: 1,000 bases
# bins overlapping by 500 bases(step size = window size = 500)
# use 8 cores (mc.cores = 8)
EMseq_bins1kb_L <- tileMethylCounts(EMseq_object_L, win.size = 1000,
                                    step.size = 500, mc.cores = 8)
EMseq_bins1kb_S <- tileMethylCounts(EMseq_object_S, win.size = 1000,
                                    step.size = 500, mc.cores = 8)


### normalize coverage between samples
# needed for comparative analysis
EMseq_bins1kb_L_norm <- normalizeCoverage(EMseq_bins1kb_L, method = "median")
EMseq_bins1kb_S_norm <- normalizeCoverage(EMseq_bins1kb_S, method = "median")


### merge files
# only pick bases covered in all samples
# information per CpG (merge strands)
EMseq_bins1kb_L_merged <- methylKit::unite(object = EMseq_bins1kb_L_norm, destrand = TRUE)
EMseq_bins1kb_S_merged <- methylKit::unite(object = EMseq_bins1kb_S_norm, destrand = TRUE)



##### convert information into data frame format #####

# make dataframe
df_temp_L <- getData(EMseq_bins1kb_L_merged)
df_temp_S <- getData(EMseq_bins1kb_S_merged)


### calculate percent DNA methylation for each bin and sample

# regions
df_EMseq_bins1kb_L_merged <- df_temp_L[,1:3]
df_EMseq_bins1kb_S_merged <- df_temp_S[,1:3]

# get sample IDs
L_IDs <- getSampleID(EMseq_bins1kb_L_merged)
S_IDs <- getSampleID(EMseq_bins1kb_S_merged)

# add DNA methylation (as number between 0 and 1)
# sample names as column names
for (i in c(1:length(L_IDs))) {
  df_EMseq_bins1kb_L_merged <- cbind(df_EMseq_bins1kb_L_merged,
                                     df_temp_L[,paste0("numCs", i)]/df_temp_L[,paste0("coverage", i)])
  colnames(df_EMseq_bins1kb_L_merged)[ncol(df_EMseq_bins1kb_L_merged)] <- L_IDs[i]
}
for (i in c(1:length(S_IDs))) {
  df_EMseq_bins1kb_S_merged <- cbind(df_EMseq_bins1kb_S_merged,
                                     df_temp_S[,paste0("numCs", i)]/df_temp_S[,paste0("coverage", i)])
  colnames(df_EMseq_bins1kb_S_merged)[ncol(df_EMseq_bins1kb_S_merged)] <- S_IDs[i]
}

remove(i)
remove(L_IDs)
remove(S_IDs)
remove(df_temp_L)
remove(df_temp_S)



##### identify regions retaining DNA methylation #####
# only for regions with >= 80% DNA methylation in control and treated cells


### extract regions with >=80% DNA methylation

# in control cells only
regions_80_LG7_ctrl <- df_EMseq_bins1kb_L_merged[which(df_EMseq_bins1kb_L_merged$Lctrl7 >= 0.8),]
regions_80_SG7_ctrl <- df_EMseq_bins1kb_S_merged[which(df_EMseq_bins1kb_S_merged$Sctrl7 >= 0.8),]

# in both control and treated cells
regions_80_LG7_ctrl_treated <- df_EMseq_bins1kb_L_merged[which(df_EMseq_bins1kb_L_merged$Lctrl7 >= 0.8 &
                                                                 df_EMseq_bins1kb_L_merged$LG7 >= 0.8),]
regions_80_SG7_ctrl_treated<- df_EMseq_bins1kb_S_merged[which(df_EMseq_bins1kb_S_merged$Sctrl7 >= 0.8 &
                                                                df_EMseq_bins1kb_S_merged$SG7 >= 0.8),]


### exclude random/undefined chromosomes and chr Y

# chromosomes to include
c <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
       "chr9", "chr10", "chr11", "chr12", "chr13", "chr14","chr15", "chr16",
       "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

# subset dataframes
regions_80_LG7_ctrl<- regions_80_LG7_ctrl[which(regions_80_LG7_ctrl$chr %in% c),]
regions_80_SG7_ctrl <- regions_80_SG7_ctrl[which(regions_80_SG7_ctrl$chr %in% c),]

regions_80_LG7_ctrl_treated <- regions_80_LG7_ctrl_treated[which(regions_80_LG7_ctrl_treated$chr %in% c),]
regions_80_SG7_ctrl_treated <- regions_80_SG7_ctrl_treated[which(regions_80_SG7_ctrl_treated$chr %in% c),]

remove(c)


### combine overlapping and adjacent regions
# also make dataframe into GRanges object

# >= 80% DNA methylation in control cells only
regions_80_LG7_ctrl_GRanges <-
  reduce(makeGRangesFromDataFrame(regions_80_LG7_ctrl))

regions_80_SG7_ctrl_GRanges <-
  reduce(makeGRangesFromDataFrame(regions_80_SG7_ctrl))

# >= 80% DNA methylation in control and treated cells
regions_80_LG7_ctrl_treated_GRanges <-
  reduce(makeGRangesFromDataFrame(regions_80_LG7_ctrl_treated))

regions_80_SG7_ctrl_treated_GRanges <-
  reduce(makeGRangesFromDataFrame(regions_80_SG7_ctrl_treated))



##### find genomic location of regions retaining DNA methylation #####

### Where can regions be found (promoter, exon, intron, intergenic)?
# promoters: 1000 bp up- and 500 bp downstream from TSS
# when a region covers several parts, locations will be prioritized as follows:
# promoter > exon > intron > intergenic


### import location of transcripts
# only curated transcripts (NM_/NR_)
gene_annot <- readTranscriptFeatures("R_imports/transcripts_curated_hg38_ucsc_20240814.bed",
                                     up.flank = 1000, down.flank = 500)


### assign genomic location
# using genomation

# all regions with >80% DNA methylation in controls
regions_80_LG7_ctrl_GRanges_annot_gene <- 
  annotateWithGeneParts(target = regions_80_LG7_ctrl_GRanges, feature = gene_annot)
regions_80_SG7_ctrl_GRanges_annot_gene <- 
  annotateWithGeneParts(target = regions_80_SG7_ctrl_GRanges, feature = gene_annot)

# regions retaining 80% DNA methylation
regions_80_LG7_ctrl_treated_GRanges_annot_gene <- 
  annotateWithGeneParts(target = regions_80_LG7_ctrl_treated_GRanges, feature = gene_annot)
regions_80_SG7_ctrl_treated_GRanges_annot_gene <- 
  annotateWithGeneParts(target = regions_80_SG7_ctrl_treated_GRanges, feature = gene_annot)


### Are regions found at transposable elements (TEs)?


### import bedfile with location of the features

# transposable elements
# provided by Sandra H.
# same annotation used for analysis of RNA-seq
TEs_annot <- readFeatureFlank("R_imports/hg38_rmsk_TE_clean.bed")


### assign feature: transposable elements

# all regions with >80% DNA methylation in controls
regions_80_LG7_ctrl_GRanges_annot_TE <- 
  annotateWithFeature(target = regions_80_LG7_ctrl_GRanges,
                      feature = TEs_annot$features, feature.name = "TE")
regions_80_SG7_ctrl_GRanges_annot_TE <- 
  annotateWithFeature(target = regions_80_SG7_ctrl_GRanges,
                      feature = TEs_annot$features, feature.name = "TE")

# regions retaining 80% DNA methylation
regions_80_LG7_ctrl_treated_GRanges_annot_TE <- 
  annotateWithFeature(target = regions_80_LG7_ctrl_treated_GRanges,
                      feature = TEs_annot$features, feature.name = "TE")
regions_80_SG7_ctrl_treated_GRanges_annot_TE <- 
  annotateWithFeature(target = regions_80_SG7_ctrl_treated_GRanges,
                      feature = TEs_annot$features, feature.name = "TE")



##### annotate regions with associated genes #####

### make dataframe
# of the reduced GRanges object (adjacent/overlapping regions are merged)

df_regions_80_LG7_ctrl <- as(regions_80_LG7_ctrl_GRanges, "data.frame")[1:4]
df_regions_80_SG7_ctrl <- as(regions_80_SG7_ctrl_GRanges, "data.frame")[1:4]

df_regions_80_LG7_ctrl_treated <- as(regions_80_LG7_ctrl_treated_GRanges, "data.frame")[1:4]
df_regions_80_SG7_ctrl_treated <- as(regions_80_SG7_ctrl_treated_GRanges, "data.frame")[1:4]

# change name of the first column to chr
colnames(df_regions_80_LG7_ctrl) [1] <- "chr"
colnames(df_regions_80_SG7_ctrl) [1] <- "chr"

colnames(df_regions_80_LG7_ctrl_treated) [1] <- "chr"
colnames(df_regions_80_SG7_ctrl_treated) [1] <- "chr"


### add genomic location
# promoter > exon > intron > intergenic

# LOUCY
df_regions_80_LG7_ctrl$location <- "intergenic"
i = 3
for (l in c("intron","exon","promoter")) {
  df_regions_80_LG7_ctrl[which(
    getMembers(regions_80_LG7_ctrl_GRanges_annot_gene)[,i] == 1),"location"] <- l
  i = i-1
}

df_regions_80_LG7_ctrl_treated$location <- "intergenic"
i = 3
for (l in c("intron","exon","promoter")) {
  df_regions_80_LG7_ctrl_treated[which(
    getMembers(regions_80_LG7_ctrl_treated_GRanges_annot_gene)[,i] == 1),"location"] <- l
  i = i-1
}

#SUP-T1
df_regions_80_SG7_ctrl$location <- "intergenic"
i = 3
for (l in c("intron","exon","promoter")) {
  df_regions_80_SG7_ctrl[which(
    getMembers(regions_80_SG7_ctrl_GRanges_annot_gene)[,i] == 1),"location"] <- l
  i = i-1
}

df_regions_80_SG7_ctrl_treated$location <- "intergenic"
i = 3
for (l in c("intron","exon","promoter")) {
  df_regions_80_SG7_ctrl_treated[which(
    getMembers(regions_80_SG7_ctrl_treated_GRanges_annot_gene)[,i] == 1),"location"] <- l
  i = i-1
}

remove(i)
remove(l)


### add transcript_id
# transcript_ID based on gene annotation file used above (gene_annot)

# get all dataframes
files <- ls(pattern = "\\bdf_regions_80_..._ctrl")

# assign transcript ID
for (f in files) {
  assign(f, cbind(get(f), "transcript_id" =
                    getAssociationWithTSS(annotateWithGeneParts(
                      target = as(get(f), "GRanges"), feature = gene_annot))
                  [,c("dist.to.feature","feature.name")]))
}

remove(files)
remove(f)

colnames(df_regions_80_LG7_ctrl)[6:7] <- c("distance","transcript_id")
colnames(df_regions_80_LG7_ctrl_treated)[6:7] <- c("distance","transcript_id")

colnames(df_regions_80_SG7_ctrl)[6:7] <- c("distance","transcript_id")
colnames(df_regions_80_SG7_ctrl_treated)[6:7] <- c("distance","transcript_id")


### add gene ID

# import gene annotation file
# use same annotation file as used for RNA-seq analysis:
# GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf

gtf <- rtracklayer::import("R_imports/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf")
df_gtf <- data.frame(gtf)

# rename first column to "chr" to be named similar to df with TE coordinates from analysis
colnames(df_gtf) [1] <- "chr"

# get all dataframes for DMR promoter analysis
files <- ls(pattern = "\\bdf_regions_80_..._ctrl")

# assign gene ID
for (f in files) {
  assign(f,
         left_join(get(f), df_gtf[which(df_gtf$type == "transcript"),c(10,11)]))
}

remove(files)
remove(f)



##### make summary dataframes #####


### percentage of regions overlapping with
# promoters > exons > introns > intergenic
# TEs
# for regions retaining DNA methylation and all regions with > 80% methylation

# make dataframe
regions_80_always_annot_summary_percent <- data.frame("region" = c("promoter","exon","intron","intergenic","TE"))

# find file names for analyzed samples
samples <- c("regions_80_LG7", "regions_80_SG7")

# add percentages for defined regions
i = 2
for (s in samples) {
  for (n in c("ctrl","ctrl_treated")) {
    regions_80_always_annot_summary_percent <- cbind(regions_80_always_annot_summary_percent,
                                                     c(get(paste0(s,"_",n,"_GRanges_annot_gene"))@precedence[["promoter"]],
                                                       get(paste0(s,"_",n,"_GRanges_annot_gene"))@precedence[["exon"]],
                                                       get(paste0(s,"_",n,"_GRanges_annot_gene"))@precedence[["intron"]],
                                                       get(paste0(s,"_",n,"_GRanges_annot_gene"))@precedence[["intergenic"]],
                                                       get(paste0(s,"_",n,"_GRanges_annot_TE"))@precedence[["TE"]]))
    colnames(regions_80_always_annot_summary_percent) [i] <- paste0(substring(s,nchar(s)-2),"_",n) 
    i = i+1
  }
}

remove(samples)
remove(i)
remove(s)
remove(n)


### total number of regions overlapping with
# promoters > exons > introns > intergenic
# TEs
# for regions retaining DNA methylation and all regions with > 80% methylation

# make dataframe
regions_80_always_annot_summary_count <- data.frame("region" = c("promoter","exon","intron","intergenic"))

# find file names for analyzed samples
samples <- ls(pattern = "^df_regions_80_")

# add percentages for defined regions
for (s in samples) {
  temp <- c()
  for (r in c("promoter","exon","intron","intergenic")) {
    temp <- c(temp,
              table(get(s)$location)[[r]])
  }
  regions_80_always_annot_summary_count <- cbind(regions_80_always_annot_summary_count, temp)
}

colnames(regions_80_always_annot_summary_count) [2:5] <- samples

remove (samples)
remove(temp)
remove(s)
remove(r)

remove(list = ls(pattern = "GRanges"))


##### export files #####

### summary files
write.table(regions_80_always_annot_summary_percent, "R_exports/tables/regions_80_always_annot_summary_percent.txt",
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
write.table(regions_80_always_annot_summary_count, "R_exports/tables/regions_80_always_annot_summary_count.txt",
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

### annotated regions retaining DNA methyltion

# as text files
write.table(df_regions_80_LG7_ctrl, "R_exports/tables/df_regions_80_LG7_ctrl.txt",
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
write.table(df_regions_80_LG7_ctrl_treated, "R_exports/tables/df_regions_80_LG7_ctrl_treated.txt",
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
write.table(df_regions_80_SG7_ctrl, "R_exports/tables/df_regions_80_SG7_ctrl.txt",
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
write.table(df_regions_80_SG7_ctrl_treated, "R_exports/tables/df_regions_80_SG7_ctrl_treated.txt",
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)




##### have a quick look at promoters retaining DNA methylation #####

### What genes retain DNA methylation?


### still >= 80% DNA methylation after treatment with GSK for 7 days
# at promoters
promoters_retain_L_gene_ids <- unique(
  df_regions_80_LG7_ctrl_treated[grep("promoter", df_regions_80_LG7_ctrl_treated$location),
                           "gene_id"])
promoters_retain_S_gene_ids <- unique(
  df_regions_80_SG7_ctrl_treated[grep("promoter", df_regions_80_SG7_ctrl_treated$location),
                           "gene_id"])
promoters_retain_both_gene_ids <- intersect(promoters_retain_L_gene_ids, promoters_retain_S_gene_ids)

# how many gene do not lose DNA methylation at their promoters?
length(promoters_retain_L_gene_ids) # 854
length(promoters_retain_S_gene_ids) # 734
length(promoters_retain_both_gene_ids) # 60


### what regions have >= 80% DNA methylation in control cells
# at promoters
promoters_80ctrl_L_gene_ids <- unique(
  df_regions_80_LG7_ctrl[grep("promoter", df_regions_80_LG7_ctrl$location),
                                 "gene_id"])
promoters_80ctrl_S_gene_ids <- unique(
  df_regions_80_SG7_ctrl[grep("promoter", df_regions_80_SG7_ctrl$location),
                                 "gene_id"])
promoters_80ctrl_both_gene_ids <- intersect(promoters_80ctrl_L_gene_ids, promoters_80ctrl_S_gene_ids)

# how many gene do not lose DNA methylation at their promoters?
length(promoters_80ctrl_L_gene_ids) # 18298
length(promoters_80ctrl_S_gene_ids) # 19000
length(promoters_80ctrl_both_gene_ids) # 14531


##### export list of genes #####

### will together become supplementary table 9

# whose promoters are >= 80% methylated in control cells
# retain methylation 
# LOUCY cells

write.table(promoters_80ctrl_L_gene_ids, "R_exports/tables/promoters_80ctrl_L_gene_ids.txt",
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(promoters_retain_L_gene_ids, "R_exports/tables/promoters_retain_L_gene_ids.txt",
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# whose promoters are >= 80% methylated in control cells
# retain methylation 
# LOUCY cells

write.table(promoters_80ctrl_S_gene_ids, "R_exports/tables/promoters_80ctrl_S_gene_ids.txt",
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(promoters_retain_S_gene_ids, "R_exports/tables/promoters_retain_S_gene_ids.txt",
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)