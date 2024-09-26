# Supplementary Table 17 & 18

# Author: Sandra Hellberg
# Email: sandra.hellberg@liu.se

# Date: 2024-02-02
# Last modified: 2024-08-14


rm(list=ls()) # remove all entries in the global environment 
set.seed(524) # seed for reproducibility

#-------------------------------------------------------------------------------

            ### Set directory structure and load packages---- 

#-------------------------------------------------------------------------------
setwd("working directory")
rna_dir <- "input directory"

lib_path <- .libPaths("C:/Rpackages/") # specifing the package library
#update.packages(lib.loc = lib_path, ask = FALSE, checkBuilt = TRUE, dependencies = TRUE)
#BiocManager::install(lib.loc = lib_path, update = TRUE, ask = FALSE)

pack_R <- c("stringr", "dplyr", "tidyr") # libraries to load

for (i in 1:length(pack_R)) {
  library(pack_R[i], character.only = TRUE)
}

#-------------------------------------------------------------------------------

                ### Read in and format input data----

#-------------------------------------------------------------------------------
#### Read in results from differential analysis
rds_input <- list.files("RDS files from RNA_seq", pattern = "_clean")

rds_names <- rds_input %>% # create a vector with the file names to be used in the for loop
  str_remove("_clean") %>%
  str_remove(".rds")

print(rds_names)

samples <- rep(c("LD3", "LG3", "LG7", "SD3", "SG3", "SG7"), 3)

for(i in 1:length(rds_input)) {
  df <- readRDS(paste0("Directory for RNA_seq RDS files", rds_input[i]))
  df$sig <- "No"
  df$sig[df$logFC > 1 & df$padjust < 0.05] <- "Upregulated"
  df$sig[df$logFC < -1 & df$padjust < 0.05] <- "Downregulated"
  df <- df[, c("Name", "TE", "chr", "start", "end" , "strand","Family", "gene", "transcript", "logFC", "padjust",  "sig")]
  df$sample <- samples[i]
  assign(rds_names[i], df)
  rm(df)
}

#### LOUCY

# Combine TEs
TE_LD3_RNA <- rbind(ERVs_LD3_RNA, LINEs_LD3_RNA, SINEs_LD3_RNA)
TE_LG3_RNA <- rbind(ERVs_LG3_RNA, LINEs_LG3_RNA, SINEs_LG3_RNA)
TE_LG7_RNA <- rbind(ERVs_LG7_RNA, LINEs_LG7_RNA, SINEs_LG7_RNA)

# Read in cpm
LOUCY_cpm <- read.delim("LOUCY_RNAseq_cpm_TE_240806.txt")

LD3_cpm <- LOUCY_cpm[,c(1,2)]
colnames(LD3_cpm) <- c("cpm_untreated", "cpm_treated")

LG3_cpm <- LOUCY_cpm[,c(1,3)]
colnames(LG3_cpm) <- c("cpm_untreated", "cpm_treated")

LG7_cpm <- LOUCY_cpm[,c(1,4)]
colnames(LG7_cpm) <- c("cpm_untreated", "cpm_treated")

# Read in lcpm
LOUCY_lcpm <- read.delim("LOUCY_RNAseq_lcpm_TE_240806.txt")

LD3_lcpm <- LOUCY_lcpm[,c(1,2)]
colnames(LD3_lcpm) <- c("lcpm_untreated", "lcpm_treated")

LG3_lcpm <- LOUCY_lcpm[,c(1,3)]
colnames(LG3_lcpm) <- c("lcpm_untreated", "lcpm_treated")

LG7_lcpm <- LOUCY_lcpm[,c(1,4)]
colnames(LG7_lcpm) <- c("lcpm_untreated", "lcpm_treated")

# combined cpm and lcpm for each treatment
LD3_combined <- cbind(LD3_cpm, LD3_lcpm)
LD3_combined$Name <- rownames(LD3_combined)

LG3_combined <- cbind(LG3_cpm, LG3_lcpm)
LG3_combined$Name <- rownames(LG3_combined)

LG7_combined <- cbind(LG7_cpm, LG7_lcpm)
LG7_combined$Name <- rownames(LG7_combined)

# merge combined with results from differential analysis
LD3 <- merge(LD3_combined, TE_LD3_RNA, by = "Name")
LD3 <- LD3[,c("Name", "TE", "chr", "start", "end", "strand", "Family", "gene", "transcript", "cpm_untreated", "cpm_treated", "lcpm_untreated","lcpm_treated", 
              "logFC", "padjust", "sig", "sample")]

LG3 <- merge(LG3_combined, TE_LG3_RNA, by = "Name")
LG3 <- LG3[,c("Name", "TE", "chr", "start", "end", "strand", "Family", "gene", "transcript", "cpm_untreated", "cpm_treated", "lcpm_untreated","lcpm_treated", 
              "logFC", "padjust", "sig", "sample")]


LG7 <- merge(LG7_combined, TE_LG7_RNA, by = "Name")
LG7 <- LG7[,c("Name", "TE", "chr", "start", "end", "strand", "Family", "gene", "transcript", "cpm_untreated", "cpm_treated", "lcpm_untreated","lcpm_treated", 
              "logFC", "padjust", "sig", "sample")]

# Read in methylation LOUCY and add relevant information regarding loss of methylation in the table
LOUCY_meth <- read.delim("Supplementary_Table_15.csv", sep = "," )

# Function to filter rows by sample pattern and significance
filter_meth <- function(data, sample_id) {
  data %>%
    filter(grepl(sample_id, sample),
           qvalue < 0.01, 
           meth.diff <= -25)
}

# Filter
# LD3
TE_LD3_meth <- filter_meth(LOUCY_meth, "LD3")
TE_LD3_meth <- TE_LD3_meth$transcript_id

LD3$methylation_loss <- ifelse(LD3$transcript %in% TE_LD3_meth, "Yes", "No")

LD3_sig <- LD3[LD3$sig == "Upregulated",] # Significantly up-regulated genes
LD3_sig_subset <- LD3_sig[LD3_sig$cpm_untreated < 0.5 & LD3_sig$cpm_treated >= 0.5,]

LD3_common <- intersect(LD3_sig_subset$transcript, TE_LD3_meth) # no common entries, none of the significantly upregulated TEs lose methylation

LD3$reexpressed <- "No"

# LG3
TE_LG3_meth <- filter_meth(LOUCY_meth, "LG3")
TE_LG3_meth <- TE_LG3_meth$transcript_id

LG3$methylation_loss <- ifelse(LG3$transcript %in% TE_LG3_meth, "Yes", "No")

LG3_sig <- LG3[LG3$sig == "Upregulated",] # Significantly up-regulated genes
LG3_sig_subset <- LG3_sig[LG3_sig$cpm_untreated < 0.5 & LG3_sig$cpm_treated >= 0.5,]

LG3_common <- intersect(LG3_sig_subset$transcript, TE_LG3_meth) 
LG3$reexpressed <- "No" # add a column called re-expressed
LG3$reexpressed[LG3$transcript %in% LG3_common] <- "Yes"


# LG7
TE_LG7_meth <- filter_meth(LOUCY_meth, "LG7")
TE_LG7_meth <- TE_LG7_meth$transcript_id

LG7$methylation_loss <- ifelse(LG7$transcript %in% TE_LG7_meth, "Yes", "No")


LG7_sig <- LG7[LG7$sig == "Upregulated",] # Significantly up-regulated genes
LG7_sig_subset <- LG7_sig[LG7_sig$cpm_untreated < 0.5 & LG7_sig$cpm_treated >= 0.5,]

LG7_common <- intersect(LG7_sig_subset$transcript, TE_LG7_meth) 
LG7$reexpressed <- "No" # add a column called re-expressed
LG7$reexpressed[LG7$transcript %in% LG7_common] <- "Yes"

LOUCY_RNA <- rbind(LD3, LG3, LG7)
colnames(LOUCY_RNA)[1] <- "full_name"
colnames(LOUCY_RNA)[8] <- "name"
colnames(LOUCY_RNA)[7] <- "family"

write.csv(LOUCY_RNA, "Supplementary_Table_17.csv", quote = F, row.names = F)

#### SUPT1

# Combine TEs
TE_SD3_RNA <- rbind(ERVs_SD3_RNA, LINEs_SD3_RNA, SINEs_SD3_RNA)
TE_SG3_RNA <- rbind(ERVs_SG3_RNA, LINEs_SG3_RNA, SINEs_SG3_RNA)
TE_SG7_RNA <- rbind(ERVs_SG7_RNA, LINEs_SG7_RNA, SINEs_SG7_RNA)

# Read in cpm
SUPT1_cpm <- read.delim("SUPT1_RNAseq_cpm_TE_240806.txt")

SD3_cpm <- SUPT1_cpm[,c(1,2)]
colnames(SD3_cpm) <- c("cpm_untreated", "cpm_treated")

SG3_cpm <- SUPT1_cpm[,c(1,3)]
colnames(SG3_cpm) <- c("cpm_untreated", "cpm_treated")

SG7_cpm <- SUPT1_cpm[,c(1,4)]
colnames(SG7_cpm) <- c("cpm_untreated", "cpm_treated")

# Read in lcpm
SUPT1_lcpm <- read.delim("SUPT1_RNAseq_lcpm_TE_240806.txt")

SD3_lcpm <- SUPT1_lcpm[,c(1,2)]
colnames(SD3_lcpm) <- c("lcpm_untreated", "lcpm_treated")

SG3_lcpm <- SUPT1_lcpm[,c(1,3)]
colnames(SG3_lcpm) <- c("lcpm_untreated", "lcpm_treated")

SG7_lcpm <- SUPT1_lcpm[,c(1,4)]
colnames(SG7_lcpm) <- c("lcpm_untreated", "lcpm_treated")

# combined cpm and lcpm for each treatment
SD3_combined <- cbind(SD3_cpm, SD3_lcpm)
SD3_combined$Name <- rownames(SD3_combined)

SG3_combined <- cbind(SG3_cpm, SG3_lcpm)
SG3_combined$Name <- rownames(SG3_combined)

SG7_combined <- cbind(SG7_cpm, SG7_lcpm)
SG7_combined$Name <- rownames(SG7_combined)

# merge combined with results from differential analysis
SD3 <- merge(SD3_combined, TE_SD3_RNA, by = "Name")
SD3 <- SD3[,c("Name", "TE", "chr", "start", "end", "strand", "Family", "gene", "transcript", "cpm_untreated", "cpm_treated", "lcpm_untreated","lcpm_treated", 
              "logFC", "padjust", "sig", "sample")]

SG3 <- merge(SG3_combined, TE_SG3_RNA, by = "Name")
SG3 <- SG3[,c("Name", "TE", "chr", "start", "end", "strand", "Family", "gene", "transcript", "cpm_untreated", "cpm_treated", "lcpm_untreated","lcpm_treated", 
              "logFC", "padjust", "sig", "sample")]

SG7 <- merge(SG7_combined, TE_SG7_RNA, by = "Name")
SG7 <- SG7[,c("Name", "TE", "chr", "start", "end", "strand", "Family", "gene", "transcript", "cpm_untreated", "cpm_treated", "lcpm_untreated","lcpm_treated", 
              "logFC", "padjust", "sig", "sample")]

# Read in methylation LOUCY and add relevant information regarding loss of methylation in the table
SUPT1_meth <- read.delim("Supplementary_Table_16.csv", sep = "," )

# Filter
# SD3
TE_SD3_meth <- filter_meth(SUPT1_meth, "SD3")
TE_SD3_meth <- TE_SD3_meth$transcript_id

SD3$methylation_loss <- ifelse(SD3$transcript %in% TE_SD3_meth, "Yes", "No")

SD3_sig <- SD3[SD3$sig == "Upregulated",] # Significantly up-regulated genes
SD3_sig_subset <- SD3_sig[SD3_sig$cpm_untreated < 0.5 & SD3_sig$cpm_treated >= 0.5,]

SD3_common <- intersect(SD3_sig_subset$transcript, TE_SD3_meth) # no common entries, none of the significantly upregulated TEs lose methylation

SD3$reexpressed <- "No"

# SG3
TE_SG3_meth <- filter_meth(SUPT1_meth, "SG3")
TE_SG3_meth <- TE_SG3_meth$transcript_id

SG3$methylation_loss <- ifelse(SG3$transcript %in% TE_SG3_meth, "Yes", "No")

SG3_sig <- SG3[SG3$sig == "Upregulated",] # Significantly up-regulated genes
SG3_sig_subset <- SG3_sig[SG3_sig$cpm_untreated < 0.5 & SG3_sig$cpm_treated >= 0.5,]

SG3_common <- intersect(SG3_sig_subset$transcript, TE_SG3_meth) 
SG3$reexpressed <- "No" # add a column called re-expressed
SG3$reexpressed[SG3$transcript %in% SG3_common] <- "Yes"

# SG7
TE_SG7_meth <- filter_meth(SUPT1_meth, "SG7")
TE_SG7_meth <- TE_SG7_meth$transcript_id

SG7$methylation_loss <- ifelse(SG7$transcript %in% TE_SG7_meth, "Yes", "No")

SG7_sig <- SG7[SG7$sig == "Upregulated",] # Significantly up-regulated genes
SG7_sig_subset <- SG7_sig[SG7_sig$cpm_untreated < 0.5 & SG7_sig$cpm_treated >= 0.5,]

SG7_common <- intersect(SG7_sig_subset$transcript, TE_SG7_meth) 
SG7$reexpressed <- "No" # add a column called re-expressed
SG7$reexpressed[SG7$transcript %in% SG7_common] <- "Yes"

SUPT1_RNA <- rbind(SD3, SG3, SG7)
colnames(SUPT1_RNA)[1] <- "full_name"
colnames(SUPT1_RNA)[8] <- "name"
colnames(SUPT1_RNA)[7] <- "family"

write.csv(SUPT1_RNA, "Supplementary_Table_18.csv", quote = F, row.names = F)

# END OF SCRIPT
