# Title: Figure 5e,f,g
# Author: Sandra Hellberg
# Email: sandra.hellberg@liu.se

# Author: Sandra Hellberg
# Email: sandra.hellberg@liu.se

# Date: 2024-08-07
# Last modified: 2024-08-15

# Description: 

rm(list=ls()) # remove all entries in the global environment 
set.seed(1002) # seed for reproducibility

#-------------------------------------------------------------------------------

          ### Set directory structure and load packages---- 

#-------------------------------------------------------------------------------
setwd("Set your working directory")
rna_dir <- "input folder for RNA-seq data"
fig_dir <- "output folder for figures"

lib_path <- .libPaths("C:/Rpackages/") # specifing the package library
#update.packages(lib.loc = lib_path, ask = FALSE, checkBuilt = TRUE, dependencies = TRUE)
#BiocManager::install(lib.loc = lib_path, update = TRUE, ask = FALSE)

pack_R <- c("dplyr", "tidyverse", "RColorBrewer", "stats", "ggplot2", "ggfortify", 
            "ggpubr", "reshape2", "ggsci", "scales", "eulerr")

for (i in 1:length(pack_R)) {
  library(pack_R[i], character.only = TRUE)
}

#-------------------------------------------------------------------------------

            ### Read in and format input data -----

#-------------------------------------------------------------------------------
#### Read data from csv
# LOUCY
rna_loucy <- read.delim(paste0(rna_dir,"Supplementary_Table_17.csv"), sep = "," ) # LOUCY
rna_LG7 <- rna_loucy[rna_loucy$sample == "LG7",] # select GSK day 7

LG7_sig_up <- rna_LG7[rna_LG7$logFC > 1 & rna_LG7$padjust < 0.05,] # upregulated TEs
LG7_sig_down <- rna_LG7[rna_LG7$logFC < -1 & rna_LG7$padjust < 0.05,] # downregulated TEs
LG7_sig <- rbind(LG7_sig_up, LG7_sig_down) # combine up and downregulated into one dataframe

ERV_LG7 <- rna_LG7[grepl("LTR|LTR\\?", rna_LG7$TE), ]
LINE_LG7 <- rna_LG7[grepl("LINE", rna_LG7$TE), ]
SINE_LG7 <- rna_LG7[grepl("SINE", rna_LG7$TE), ]


ERVs_LG7_sig <- LG7_sig[grepl("LTR|LTR\\?", LG7_sig$TE), ]
LINEs_LG7_sig <- LG7_sig[grepl("LINE", LG7_sig$TE), ]
SINEs_LG7_sig <- LG7_sig[grepl("SINE", LG7_sig$TE), ]

# SUPT1
rna_supt1 <- read.delim(paste0(rna_dir,"Supplementary_Table_18.csv"), sep = "," ) # SUPT1
rna_SG7 <- rna_supt1[rna_supt1$sample == "SG7",] # select GSK day 7

SG7_sig_up <- rna_SG7[rna_SG7$logFC > 1 & rna_SG7$padjust < 0.05,] # upregulated TEs
SG7_sig_down <- rna_SG7[rna_SG7$logFC < -1 & rna_SG7$padjust < 0.05,] # downregulated TEs
SG7_sig <- rbind(SG7_sig_up, SG7_sig_down) # combine up and downregulated into one dataframe

ERV_SG7 <- rna_SG7[grepl("LTR|LTR\\?", rna_SG7$TE), ]
LINE_SG7 <- rna_SG7[grepl("LINE", rna_SG7$TE), ]
SINE_SG7 <- rna_SG7[grepl("SINE", rna_SG7$TE), ]

ERVs_SG7_sig <- SG7_sig[grepl("LTR|LTR\\?", SG7_sig$TE), ]
LINEs_SG7_sig <- SG7_sig[grepl("LINE", SG7_sig$TE), ]
SINEs_SG7_sig <- SG7_sig[grepl("SINE", SG7_sig$TE), ]


#-------------------------------------------------------------------------------

                  ### # Test difference between distributions-----

#-------------------------------------------------------------------------------
# LOUCY
# ERV vs LINE
upregulated_L<- c(sum(ERVs_LG7_sig$sig == "Upregulated"), sum(LINEs_LG7_sig$sig == "Upregulated"))  # Number of successes in each group
sample_sizes_L <- c(nrow(ERV_LG7), nrow(LINE_LG7))  # Sample sizes for each group

# Perform the test
result_L <- prop.test(upregulated_L, sample_sizes_L, alternative = "greater")

# Print the result
print(result_L)

# ERV vs SINE
upregulated_L<- c(sum(ERVs_LG7_sig$sig == "Upregulated"), sum(SINEs_LG7_sig$sig == "Upregulated"))  # Number of successes in each group
sample_sizes_L <- c(nrow(ERV_LG7), nrow(SINE_LG7))  # Sample sizes for each group

# Perform the test
result_L <- prop.test(upregulated_L, sample_sizes_L, alternative = "greater")

# Print the result
print(result_L)

# SUP-T1
# ERV vs LINE
upregulated_S<- c(sum(ERVs_SG7_sig$sig == "Upregulated"), sum(LINEs_SG7_sig$sig == "Upregulated"))  # Number of successes in each group
sample_sizes_S <- c(nrow(ERV_SG7), nrow(LINE_SG7))  # Sample sizes for each group

# Perform the test
result_S <- prop.test(upregulated_S, sample_sizes_S, alternative = "greater")

# Print the result
print(result_S)

# ERV vs SINE
upregulated_S<- c(sum(ERVs_SG7_sig$sig == "Upregulated"), sum(SINEs_SG7_sig$sig == "Upregulated"))  # Number of successes in each group
sample_sizes_S <- c(nrow(ERV_SG7), nrow(SINE_SG7))  # Sample sizes for each group

# Perform the test
result_S <- prop.test(upregulated_S, sample_sizes_S, alternative = "greater")

# Print the result
print(result_S)

#-------------------------------------------------------------------------------

                      ### Plotting -----

#-------------------------------------------------------------------------------
#### Fig 5E----
LG7_TE <- rbind(ERVs_LG7_sig, LINEs_LG7_sig, SINEs_LG7_sig)
LG7_TE <- LG7_TE[LG7_TE$sig == "Upregulated",]

SG7_TE <- rbind(ERVs_SG7_sig, LINEs_SG7_sig, SINEs_SG7_sig)
SG7_TE <- SG7_TE[SG7_TE$sig == "Upregulated",]

combined_list <- list(
  "LG7" = LG7_TE$transcript,
  "SG7" = SG7_TE$transcript)

fit_combined <- euler(combined_list)
combined_plot <- plot(fit_combined, fills = TRUE, edges = TRUE, quantities = TRUE)

print(combined_plot)

# Export as pdf

# check classes of intersecting TEs

class_TE <- merge(LG7_TE, SG7_TE, by.x = "transcript", by.y = "transcript")
table(class_TE$TE.x)

#### Fig 5F----
cols <- c("Upregulated" = "#8B0000", "Downregulated" = "blue") 

p1 <- ggplot(data=ERVs_LG7_sig, aes(x=logFC, y=-log10(padjust), col=sig)) + 
  geom_point(shape = 21, # Specify shape and colour as fixed local parameters
             fill = "black",
             size = 2) +  # change colour to fill
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  geom_vline(xintercept = 1) + 
  theme_bw() +
  scale_color_manual(values = cols) +  # Use scale_color_manual instead of scale_fill_manual
  ylim(-log10(0.05), 7) +
  xlim(-13, 15) +
  labs(title = "LOUCY_ERV_G7")

print(p1)

p2 <- ggplot(data=LINEs_LG7_sig, aes(x=logFC, y=-log10(padjust), col=sig)) + 
  geom_point(shape = 21, # Specify shape and colour as fixed local parameters
             fill = "black",
             size = 2) +  # change colour to fill
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  geom_vline(xintercept = 1) + 
  theme_bw() +
  scale_color_manual(values = cols) +  # Use scale_color_manual instead of scale_fill_manual
  ylim(-log10(0.05), 7) +
  xlim(-13, 15) +
  labs(title = "LOUCY_LINE_G7")

print(p2)

p3 <- ggplot(data=SINEs_LG7_sig, aes(x=logFC, y=-log10(padjust), col=sig)) + 
  geom_point(shape = 21, # Specify shape and colour as fixed local parameters
             fill = "black",
             size = 2) +  # change colour to fill
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  geom_vline(xintercept = 1) + 
  theme_bw() +
  scale_color_manual(values = cols) +  # Use scale_color_manual instead of scale_fill_manual
  ylim(-log10(0.05), 7) +
  xlim(-13, 15) +
  labs(title = "LOUCY_SINE_G7")

print(p3)


p4 <- ggplot(data=ERVs_SG7_sig, aes(x=logFC, y=-log10(padjust), col=sig)) + 
  geom_point(shape = 21, # Specify shape and colour as fixed local parameters
             fill = "black",
             size = 2) +  # change colour to fill
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  geom_vline(xintercept = 1) + 
  theme_bw() +
  scale_color_manual(values = cols) +  # Use scale_color_manual instead of scale_fill_manual
  ylim(-log10(0.05), 7) +
  xlim(-13, 15) +
  labs(title = "SUPT1_ERV_G7")

print(p4)


p5 <- ggplot(data=LINEs_SG7_sig, aes(x=logFC, y=-log10(padjust), col=sig)) + 
  geom_point(shape = 21, # Specify shape and colour as fixed local parameters
             fill = "black",
             size = 2) +  # change colour to fill
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  geom_vline(xintercept = 1) + 
  theme_bw() +
  scale_color_manual(values = cols) +  # Use scale_color_manual instead of scale_fill_manual
  ylim(-log10(0.05), 7) +
  xlim(-13, 15) +
  labs(title = "SUPT1_LINE_G7")

print(p5)


p6 <- ggplot(data=SINEs_SG7_sig, aes(x=logFC, y=-log10(padjust), col=sig)) + 
  geom_point(shape = 21, # Specify shape and colour as fixed local parameters
             fill = "black",
             size = 2) +  # change colour to fill
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  geom_vline(xintercept = 1) + 
  theme_bw() +
  scale_color_manual(values = cols) +  # Use scale_color_manual instead of scale_fill_manual
  ylim(-log10(0.05), 7) +
  xlim(-13, 15) +
  labs(title = "SUPT1_SINE_G7")

print(p6)


# Export as pdf
plot <- ggarrange(p1, p2, p3, p4, p5, p6, nrow=2, ncol=3)

pdf(paste0(fig_dir,"Fig5F_240815.pdf"),
    width = 13,
    height = 5)
plot(plot)
dev.off()


