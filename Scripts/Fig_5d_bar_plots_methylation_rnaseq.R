# Title: Figure 5D
# Author: Sandra Hellberg
# Email: sandra.hellberg@liu.se

# Date: 2024-08-07
# Last modified: 2024-08-14

# Description: 

rm(list=ls()) # remove all entries in the global environment 
set.seed(1002) # seed for reproducibility

#-------------------------------------------------------------------------------

              ### Set directory structure and load packages 

#-------------------------------------------------------------------------------
setwd("set your working directory")
meth_dir <- "input folder for methylation data"
rna_dir <- "input folder for RNAseq data"
fig_dir <- "Output folder for figures"

lib_path <- .libPaths("C:/Rpackages/") # specifing the package library
#update.packages(lib.loc = lib_path, ask = FALSE, checkBuilt = TRUE, dependencies = TRUE)
#BiocManager::install(lib.loc = lib_path, update = TRUE, ask = FALSE)

pack_R <- c("dplyr", "tidyverse", "RColorBrewer", "stats", "ggplot2", "ggfortify", "ggpubr", "reshape2", "ggsci", "scales")

for (i in 1:length(pack_R)) {
  library(pack_R[i], character.only = TRUE)
}

#-------------------------------------------------------------------------------

                              ## LOUCY----

#-------------------------------------------------------------------------------
#### Read in and format input data
rna_loucy <- read.delim(paste0(rna_dir,"Supplementary_Table_17.csv"), sep = "," ) # LOUCY
rna_LG7 <- rna_loucy[rna_loucy$sample == "LG7",]
LG7_sig <- rna_LG7[rna_LG7$logFC > 1 & rna_LG7$padjust < 0.05,]

meth_loucy <- read.delim(paste0(meth_dir,"Supplementary_Table_15.csv"), sep = ",") # LOUCY
meth_LG7 <- meth_loucy[meth_loucy$sample == "LG7",]
meth_LG7_sig <- meth_LG7[meth_LG7$qvalue < 0.01 & meth_LG7$meth.diff <= -25,]

common_L <- merge(LG7_sig, meth_LG7_sig, by.x = "transcript", by.y = "transcript_id") # overlapping TEs between methylation and RNA-seq

#### RNAseq
# Function to filter rows by sample pattern and class_id type, and select specific columns
filter_rna <- function(data, sample_id) {
  data %>%
    filter(grepl(sample_id, sample)) %>%
    select(family, transcript, lcpm_treated, lcpm_untreated)  # Select only the `transcript_id` column
}
# Filter
TE_LD3_RNA <- filter_rna(rna_loucy, "LD3")
TE_LD3_RNA <- TE_LD3_RNA[TE_LD3_RNA$transcript %in% common_L$transcript,]

TE_LG3_RNA <- filter_rna(rna_loucy, "LG3")
TE_LG3_RNA <- TE_LG3_RNA[TE_LG3_RNA$transcript %in% common_L$transcript,]

TE_LG7_RNA <- filter_rna(rna_loucy, "LG7")
TE_LG7_RNA <- TE_LG7_RNA[TE_LG7_RNA$transcript %in% common_L$transcript,]

# combined lcpm data frame 
combined_rna <- data.frame(cbind(TE_LD3_RNA$lcpm_untreated, TE_LD3_RNA$lcpm_treated, TE_LG3_RNA$lcpm_treated, TE_LG7_RNA$lcpm_treated))
colnames(combined_rna) <- c("untreated", "D3", "GSK3", "GSK7")
rownames(combined_rna) <- TE_LG7_RNA$transcript


#### Methylation
# Function to filter rows by sample pattern and class_id type, and select specific columns
filter_meth <- function(data, sample_id) {
  data %>%
    filter(grepl(sample_id, sample)) %>%
    select(class_id, transcript_id, qvalue, meth.diff, methylation_treated, methylation_control)  # Select columns to keep
}

# Filter
TE_LD3_meth <- filter_meth(meth_loucy, "LD3")
TE_LD3_meth <- TE_LD3_meth[TE_LD3_meth$transcript_id %in% common_L$transcript,]

TE_LG3_meth <- filter_meth(meth_loucy, "LG3")
TE_LG3_meth <- TE_LG3_meth[TE_LG3_meth$transcript_id %in% common_L$transcript,]

TE_LG7_meth <- filter_meth(meth_loucy, "LG7")
TE_LG7_meth <- TE_LG7_meth[TE_LG7_meth$transcript_id %in% common_L$transcript,]


# combined methylation data frame 
combined_meth <- data.frame(cbind(TE_LD3_meth$methylation_control, TE_LD3_meth$methylation_treated, TE_LG3_meth$methylation_treated, TE_LG7_meth$methylation_treated))
colnames(combined_meth) <- c("untreated", "D3", "GSK3", "GSK7")
rownames(combined_meth) <- TE_LG7_meth$transcript_id

#-------------------------------------------------------------------------------
                         ## Plotting
#-------------------------------------------------------------------------------
npg_colors <- pal_npg()(10)
show_col (npg_colors)
nejm_colors <- pal_nejm()(8)
show_col (nejm_colors)

## Box plot of methylation-------------------------------------------------------
combined_meth$family <- rownames(combined_meth) 
L_meth_melt <- melt(combined_meth,  id.vars = "family")


p1 <- ggplot(data = L_meth_melt, aes(x = variable, y = value)) +
  geom_jitter(aes(color = variable), size = 0.5) +  # Map 'variable' to color
  geom_boxplot(outlier.shape = NA, color = "black", fill = NA) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", ref.group = "untreated") +
  labs(title = "LOUCY methylation", y = "% CpG methylation") +
  scale_color_manual(
    values = c(
      "untreated" = "azure3",
      "D3" = nejm_colors[7],
      "GSK3" = npg_colors[10],
      "GSK7" = npg_colors[3]
    )
  ) +
  theme_bw()

print(p1)

## Box plot of RNA-seq-----------------------------------------------------------
combined_rna$family <- rownames(combined_rna)
L_rna_melt <- melt(combined_rna,  id.vars = "family")

p2 <- ggplot(data = L_rna_melt, aes(x = variable, y = value)) +
  geom_jitter(aes(color = variable), size = 0.5) +  # Map 'variable' to color
  geom_boxplot(outlier.shape = NA, color = "black", fill = NA) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", ref.group = "untreated") +
  labs(title = "LOUCY expression", y = "expression (log2(cpm))") +
  ylim(-5,15) +
  scale_color_manual(
    values = c(
      "untreated" = "azure3",
      "D3" = nejm_colors[7],
      "GSK3" = npg_colors[10],
      "GSK7" = npg_colors[3]
    )
  ) +
  theme_bw()


print(p2)


#-------------------------------------------------------------------------------

                                ## SUP-T1----

#-------------------------------------------------------------------------------
#### Read in and format input data
rna_supt1 <- read.delim(paste0(rna_dir,"Supplementary_Table_18.csv"), sep = "," ) # SUPT1
rna_SG7 <- rna_supt1[rna_supt1$sample == "SG7",]
SG7_sig <- rna_SG7[rna_SG7$logFC > 1 & rna_SG7$padjust < 0.05,]

meth_supt1 <- read.delim(paste0(meth_dir,"Supplementary_Table_16.csv"), sep = ",") # SUPT1
meth_SG7 <- meth_supt1[meth_supt1$sample == "SG7",]
meth_SG7_sig <- meth_SG7[meth_SG7$qvalue < 0.01 & meth_SG7$meth.diff <= -25,]

common_S <- merge(SG7_sig, meth_SG7_sig, by.x = "transcript", by.y = "transcript_id") # overlapping TEs between methylation and RNA-seq

#### RNAseq
# Filter
TE_SD3_RNA <- filter_rna(rna_supt1, "SD3")
TE_SD3_RNA <- TE_SD3_RNA[TE_SD3_RNA$transcript %in% common_S$transcript,]

TE_SG3_RNA <- filter_rna(rna_supt1, "SG3")
TE_SG3_RNA <- TE_SG3_RNA[TE_SG3_RNA$transcript %in% common_S$transcript,]

TE_SG7_RNA <- filter_rna(rna_supt1, "SG7")
TE_SG7_RNA <- TE_SG7_RNA[TE_SG7_RNA$transcript %in% common_S$transcript,]

# combined lcpm data frame 
combined_rna <- data.frame(cbind(TE_SD3_RNA$lcpm_untreated, TE_SD3_RNA$lcpm_treated, TE_SG3_RNA$lcpm_treated, TE_SG7_RNA$lcpm_treated))
colnames(combined_rna) <- c("untreated", "D3", "GSK3", "GSK7")
rownames(combined_rna) <- TE_SG7_RNA$transcript


#### Methylation
# Filter
TE_SD3_meth <- filter_meth(meth_supt1, "SD3")
TE_SD3_meth <- TE_SD3_meth[TE_SD3_meth$transcript_id %in% common_S$transcript,]

TE_SG3_meth <- filter_meth(meth_supt1, "SG3")
TE_SG3_meth <- TE_SG3_meth[TE_SG3_meth$transcript_id %in% common_S$transcript,]

TE_SG7_meth <- filter_meth(meth_supt1, "SG7")
TE_SG7_meth <- TE_SG7_meth[TE_SG7_meth$transcript_id %in% common_S$transcript,]


# combined methylation data frame 
combined_meth <- data.frame(cbind(TE_SD3_meth$methylation_control, TE_SD3_meth$methylation_treated, TE_SG3_meth$methylation_treated, TE_SG7_meth$methylation_treated))
colnames(combined_meth) <- c("untreated", "D3", "GSK3", "GSK7")
rownames(combined_meth) <- TE_SG7_meth$transcript_id

#-------------------------------------------------------------------------------
                                      ## Plotting
#-------------------------------------------------------------------------------

## Box plot of methylation-------------------------------------------------------
combined_meth$family <- rownames(combined_meth) 
S_meth_melt <- melt(combined_meth,  id.vars = "family")


p3 <- ggplot(data = S_meth_melt, aes(x = variable, y = value)) +
  geom_jitter(aes(color = variable), size = 0.5) +  # Map 'variable' to color
  geom_boxplot(outlier.shape = NA, color = "black", fill = NA) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", ref.group = "untreated") +
  labs(title = "SUPT1 methylation", y = "% CpG methylation") +
  scale_color_manual(
    values = c(
      "untreated" = "azure3",
      "D3" = nejm_colors[7],
      "GSK3" = npg_colors[10],
      "GSK7" = npg_colors[3]
    )
  ) +
  theme_bw()

print(p3)

## Box plot of RNA-seq-----------------------------------------------------------
combined_rna$family <- rownames(combined_rna)
S_rna_melt <- melt(combined_rna,  id.vars = "family")

p4 <- ggplot(data = S_rna_melt, aes(x = variable, y = value)) +
  geom_jitter(aes(color = variable), size = 0.5) +  # Map 'variable' to color
  geom_boxplot(outlier.shape = NA, color = "black", fill = NA) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", ref.group = "untreated") +
  labs(title = "SUPT1 expression", y = "expression (log2(cpm))") +
  ylim(-5,15) +
  scale_color_manual(
    values = c(
      "untreated" = "azure3",
      "D3" = nejm_colors[7],
      "GSK3" = npg_colors[10],
      "GSK7" = npg_colors[3]
    )
  ) +
  theme_bw()


print(p4)

#-------------------------------------------------------------------------------
                        ## Export plots----
#-------------------------------------------------------------------------------

plot <- ggarrange(p1,p3,p2,p4, ncol= 2, nrow=2)
pdf(paste0(fig_dir, "Fig5D_barplots.pdf"), width = 8, height = 7)
plot(plot)
dev.off()


##### SCRIPT ENDS HERE


