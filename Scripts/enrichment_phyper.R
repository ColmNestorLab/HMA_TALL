################################################################################
############ ENRICHMENT FOR IMMUNE GENES, CTAs, TSGs AND ONCOGENES ############# 
################################################################################

# Maike Bensberg (maike.bensberg@liu.se)
# 2024-08-20


### lists of genes to test
# differentially expressed genes (defined in script: combining_DNAmethylation_and_expression_data.R)
# genes sensitive to DNA de-methylation (defined in script: combining_DNAmethylation_and_expression_data.R)
# test all treatments for differential expression and only methylation sensitive genes after 7 days treatment


### check enrichment of gene lists for
# immune genes
# Cancer/testis anigenes (CTAs)
# tumour suppressor genes (TSGs)
# Oncogenes


### Immune genes, CTAs, TSGs and oncogenes summarized in Supplementary Table 12
# list of immune genes curated by Sandra H. (citations ????)
# TSGs and oncogenes based on the Network of cancer genes & healthy drivers 
# http://network-cancer-genes.org/citation.php
# CTAs based on publication by Carter et al, 2023
# https://jitc.bmj.com/content/11/12/e007935


### hypergeometric test used (phyper in base R)

### also generate a plot the adjusted p-values
# Extended Data Fig. 7a



##### initiate R environment #####
# only needs to be done once
#install.packages("renv")
#renv::init()



##### loading packages #####

library(readxl)
library(ggpubr)
library(ggsci)
library(scales)


##### set working directory and check package versions #####

# set working directory
setwd("/mnt/wwn-0x5000cca28fd30c7d-part1/HMA_TALL/EMseq/LOUCY_SUPT1_HMAtreated_12_2023/methylKit")
# print working directory
getwd()

# print version of R and attached packages
sessionInfo()



##### reading in data #####

### differentially expressed genes
# Supplementary Table 10
df_expression_combined <- read.delim("R_exports/tables/Supplementary_Table_10.csv", sep = ",")

### DNA methylation sensitive genes
# Supplementary Table 11
methylation_expression_N_transcripts <- read.delim("R_exports/tables/Supplementary_Table_11.csv", sep = ",")

### lists of genes
# Supplementary Table 12
genes_of_interest <- as.data.frame(read_excel("R_imports/Supplementary_Table_12.xlsx", 
                                              skip = 1))



##### hypergeometric test: differentially expressed genes #####
# hypergeometric test (phyper in base R)


### create dataframe
# list: list of genes of interest (immune/CTA/TSG/onco)
# treatment

df_enrichment_phyper_DEG <-
  data.frame("list" = (c(rep("Immune",6), rep("CTA",6), rep("TSG",6), rep("onco",6))))
df_enrichment_phyper_DEG$treatment <- rep(c("LD3", "LG3", "LG7", "SD3", "SG3", "SG7"),4)


### define values for
# q: number of genes in the gene list that are upregulated
# m: number of genes in the gene list that are covered by RNA-seq
# n: number of genes that are covered by RNA-seq but are not part of the gene list
# k: total number of genes that are upregulated

q <- c()
m <- c()
n <- c()
k <- c()

for (c in c("Immune", "CTA", "TSG", "onco")) {
  for (t in c("LD3", "LG3", "LG7", "SD3", "SG3", "SG7")) {
    q <- c(q,
           length(unique(df_expression_combined[which(
             df_expression_combined$sample == t &
               df_expression_combined$expression == "upregulated" &
               df_expression_combined$gene_id %in%
               genes_of_interest[which(genes_of_interest$Category == c),
                                 "Gene"]),"gene_id"])))
    m <- c(m,
           length(unique(df_expression_combined[which(
             df_expression_combined$sample == t &
               df_expression_combined$gene_id %in%
               genes_of_interest[which(genes_of_interest$Category == c),
                                 "Gene"]),"gene_id"])))
    n <- c(n,
           length(unique(df_expression_combined[which(
             df_expression_combined$sample == t &
               !(df_expression_combined$gene_id %in%
               genes_of_interest[which(genes_of_interest$Category == c),
                                 "Gene"])),"gene_id"])))
    k <- c(k,
           length(unique(df_expression_combined[which(
             df_expression_combined$sample == t &
               df_expression_combined$expression == "upregulated"),"gene_id"])))
  }
}

df_enrichment_phyper_DEG$q <- q
df_enrichment_phyper_DEG$m <- m
df_enrichment_phyper_DEG$n <- n
df_enrichment_phyper_DEG$k <- k

remove(c)
remove(q)
remove(m)
remove(n)
remove(k)


### run test (calculate p-value)

p <- c()

for (i in 1:nrow(df_enrichment_phyper_DEG)) {
  p <- c(p,
         phyper(q = df_enrichment_phyper_DEG[i, "q"]-1,
                m = df_enrichment_phyper_DEG[i, "m"],
                n = df_enrichment_phyper_DEG[i, "n"],
                k = df_enrichment_phyper_DEG[i, "k"],
                lower.tail = FALSE))
}

df_enrichment_phyper_DEG$p <- p

remove(p)


### correct p-value for multiple testing
# Benjamini Hochberg correction
# round adjusted p-value to 5 decimals

df_enrichment_phyper_DEG$p_adjust <- p.adjust(p = df_enrichment_phyper_DEG$p, method = "BH")



##### check enrichment for random gene list #####
# take 200 random genes that are covered by RNA-seq
# calculate enrichment for these random genes (hypergeometric test)
# adjust p-value based only for the six samples (don't combine with other enrichment analyses)
# bootstrap sample randomization 1000 times

### calculate parameters and run hypergeometric test (phyper)
# 1000 iterations

for (t in c("LD3", "LG3", "LG7", "SD3", "SG3", "SG7")) {
  assign(paste0("q_",t), c())
  assign(paste0("m_",t), c())
  assign(paste0("n_",t), c())
  assign(paste0("k_",t), c())
  assign(paste0("p_",t), c())
}

for (i in c(1:1000)) {
  g <- sample(unique(df_expression_combined$gene_id),200)
  for (t in c("LD3", "LG3", "LG7", "SD3", "SG3", "SG7")) {
    assign(paste0("q_",t),
           c(get(paste0("q_",t)),
             length(unique(df_expression_combined[which(
               df_expression_combined$sample == t &
                 df_expression_combined$expression == "upregulated" &
                 df_expression_combined$gene_id %in% g),"gene_id"]))))
    assign(paste0("m_",t),
           c(get(paste0("m_",t)),
             length(unique(df_expression_combined[which(
               df_expression_combined$sample == t &
                 df_expression_combined$gene_id %in% g),"gene_id"]))))
    assign(paste0("n_",t),
           c(get(paste0("n_",t)),
             length(unique(df_expression_combined[which(
               df_expression_combined$sample == t &
                 !(df_expression_combined$gene_id %in% g)),"gene_id"]))))
    assign(paste0("k_",t),
           c(get(paste0("k_",t)),
             length(unique(df_expression_combined[which(
               df_expression_combined$sample == t &
                 df_expression_combined$expression == "upregulated"),"gene_id"]))))
    
    assign(paste0("p_",t),
           c(get(paste0("p_",t)),
             phyper(q = get(paste0("q_",t))[i]-1,
                    m = get(paste0("m_",t))[i],
                    n = get(paste0("n_",t))[i],
                    k = get(paste0("k_",t))[i],
                    lower.tail = FALSE)))
  }
}

remove(i)
remove(t)
remove(g)


### make dataframe with average p values for each treatment

random_list_enrichment <- data.frame("list" = rep("random", 6))
random_list_enrichment$treatment <- c("LD3", "LG3", "LG7", "SD3", "SG3", "SG7")

q <- c()
m <- c()
n <- c()
k <- c()
p <- c()

for (t in c("LD3", "LG3", "LG7", "SD3", "SG3", "SG7")) {
  q <- c(q, mean(get(paste0("q_",t))))
  m <- c(m, mean(get(paste0("m_",t))))
  n <- c(n, mean(get(paste0("n_",t))))
  k <- c(k, mean(get(paste0("k_",t))))
  p <- c(p, mean(get(paste0("p_",t))))
}

random_list_enrichment$q <- q
random_list_enrichment$m <- m
random_list_enrichment$n <- n
random_list_enrichment$k <- k
random_list_enrichment$p <- p

# remove temporary variables
remove(t)
remove(list = ls(pattern = "\\bq_...\\b"))
remove(list = ls(pattern = "\\bm_...\\b"))
remove(list = ls(pattern = "\\bn_...\\b"))
remove(list = ls(pattern = "\\bk_...\\b"))
remove(list = ls(pattern = "\\bp_...\\b"))
remove(q)
remove(m)
remove(n)
remove(k)
remove(p)


### adjust p-values (Benjamini Hochberg)
# multiple treatments -> multiple tests

random_list_enrichment$p_adjust <- p.adjust(p = random_list_enrichment$p, method = "BH")


### add permutations to summary dataframe

df_enrichment_phyper_DEG <-
  rbind(df_enrichment_phyper_DEG, random_list_enrichment)

remove(random_list_enrichment)



##### hypergeometric test: methylation sensitive genes #####
# hypergeometric test (phyper in base R)


### create dataframe
# list: list of genes of interest (immune/CTA/TSG/onco)
# treatment

df_enrichment_phyper_methyl_sensitive <-
  data.frame("list" = (c(rep("Immune",2), rep("CTA",2), rep("TSG",2), rep("onco",2))))
df_enrichment_phyper_methyl_sensitive$treatment <- rep(c("LG7", "SG7"),4)


### define values for
# q: number of genes in the gene list that are methylation-sensitive
# m: number of genes in the gene list that lose DNA methylation (DNA methylation difference ≤ -25% as compared to untreated cells (adjusted p-value < 0.01))
# n: number of genes that lose DNA methylation but are not part of the gene list
# k: total number of genes that are methylation-sensitive

q <- c()
m <- c()
n <- c()
k <- c()

for (c in c("Immune", "CTA", "TSG", "onco")) {
  for (t in c("LG7", "SG7")) {
    q <- c(q,
           length(unique(methylation_expression_N_transcripts[which(
             methylation_expression_N_transcripts$sample == t &
               methylation_expression_N_transcripts$methyl_sensitive == "yes" &
               methylation_expression_N_transcripts$gene_id %in%
               genes_of_interest[which(genes_of_interest$Category == c),
                                 "Gene"]),"gene_id"])))
    m <- c(m,
           length(unique(methylation_expression_N_transcripts[which(
             methylation_expression_N_transcripts$sample == t &
               methylation_expression_N_transcripts$methylation == "loss" &
               methylation_expression_N_transcripts$gene_id %in%
               genes_of_interest[which(genes_of_interest$Category == c),
                                 "Gene"]),"gene_id"])))
    n <- c(n,
           length(unique(methylation_expression_N_transcripts[which(
             methylation_expression_N_transcripts$sample == t &
               methylation_expression_N_transcripts$methylation == "loss" &
               !(methylation_expression_N_transcripts$gene_id %in%
                   genes_of_interest[which(genes_of_interest$Category == c),
                                     "Gene"])),"gene_id"])))
    k <- c(k,
           length(unique(methylation_expression_N_transcripts[which(
             methylation_expression_N_transcripts$sample == t &
               methylation_expression_N_transcripts$methyl_sensitive == "yes"),"gene_id"])))
  }
}

df_enrichment_phyper_methyl_sensitive$q <- q
df_enrichment_phyper_methyl_sensitive$m <- m
df_enrichment_phyper_methyl_sensitive$n <- n
df_enrichment_phyper_methyl_sensitive$k <- k

remove(q)
remove(m)
remove(n)
remove(k)


### run test (calculate p-value)

p <- c()

for (i in 1:nrow(df_enrichment_phyper_methyl_sensitive)) {
  p <- c(p,
         phyper(q = df_enrichment_phyper_methyl_sensitive[i, "q"]-1,
                m = df_enrichment_phyper_methyl_sensitive[i, "m"],
                n = df_enrichment_phyper_methyl_sensitive[i, "n"],
                k = df_enrichment_phyper_methyl_sensitive[i, "k"],
                lower.tail = FALSE))
}

df_enrichment_phyper_methyl_sensitive$p <- p

remove(p)


### correct p-value for multiple testing
# Benjamini Hochberg correction
# round adjusted p-value to 5 decimals

df_enrichment_phyper_methyl_sensitive$p_adjust <- 
  p.adjust(p = df_enrichment_phyper_methyl_sensitive$p, method = "BH")



##### check enrichment for random gene list #####
# take 200 random genes that are covered by RNA-seq
# calculate enrichment for these random genes (hypergeometric test)
# adjust p-value based only for the six samples (don't combine with other enrichment analyses)
# bootstrap sample ranodmization 1000 times

### calculate parameters and run hypergeometric test (phyper)
# 1000 iterations

for (t in c("LG7", "SG7")) {
  assign(paste0("q_",t), c())
  assign(paste0("m_",t), c())
  assign(paste0("n_",t), c())
  assign(paste0("k_",t), c())
  assign(paste0("p_",t), c())
}

for (i in c(1:1000)) {
  g <- sample(unique(methylation_expression_N_transcripts$gene_id),200)
  for (t in c("LG7", "SG7")) {
    assign(paste0("q_",t),
           c(get(paste0("q_",t)),
             length(unique(methylation_expression_N_transcripts[which(
               methylation_expression_N_transcripts$sample == t &
                 methylation_expression_N_transcripts$expression == "upregulated" &
                 methylation_expression_N_transcripts$gene_id %in% g),"gene_id"]))))
    assign(paste0("m_",t),
           c(get(paste0("m_",t)),
             length(unique(methylation_expression_N_transcripts[which(
               methylation_expression_N_transcripts$sample == t &
                 methylation_expression_N_transcripts$gene_id %in% g),"gene_id"]))))
    assign(paste0("n_",t),
           c(get(paste0("n_",t)),
             length(unique(methylation_expression_N_transcripts[which(
               methylation_expression_N_transcripts$sample == t &
                 !(methylation_expression_N_transcripts$gene_id %in% g)),"gene_id"]))))
    assign(paste0("k_",t),
           c(get(paste0("k_",t)),
             length(unique(methylation_expression_N_transcripts[which(
               methylation_expression_N_transcripts$sample == t &
                 methylation_expression_N_transcripts$expression == "upregulated"),"gene_id"]))))
    
    assign(paste0("p_",t),
           c(get(paste0("p_",t)),
             phyper(q = get(paste0("q_",t))[i]-1,
                    m = get(paste0("m_",t))[i],
                    n = get(paste0("n_",t))[i],
                    k = get(paste0("k_",t))[i],
                    lower.tail = FALSE)))
  }
}

remove(i)
remove(t)
remove(g)


### make dataframe with average p values for each treatment

random_list_enrichment <- data.frame("list" = rep("random", 2))
random_list_enrichment$treatment <- c("LG7", "SG7")

q <- c()
m <- c()
n <- c()
k <- c()
p <- c()

for (t in c("LG7", "SG7")) {
  q <- c(q, mean(get(paste0("q_",t))))
  m <- c(m, mean(get(paste0("m_",t))))
  n <- c(n, mean(get(paste0("n_",t))))
  k <- c(k, mean(get(paste0("k_",t))))
  p <- c(p, mean(get(paste0("p_",t))))
}

random_list_enrichment$q <- q
random_list_enrichment$m <- m
random_list_enrichment$n <- n
random_list_enrichment$k <- k
random_list_enrichment$p <- p

# remove temporary variables
remove(t)
remove(list = ls(pattern = "\\bq_...\\b"))
remove(list = ls(pattern = "\\bm_...\\b"))
remove(list = ls(pattern = "\\bn_...\\b"))
remove(list = ls(pattern = "\\bk_...\\b"))
remove(list = ls(pattern = "\\bp_...\\b"))
remove(q)
remove(m)
remove(n)
remove(k)
remove(p)


### adjust p-values (Benjamini Hochberg)
# multiple treatments -> multiple tests

random_list_enrichment$p_adjust <- p.adjust(p = random_list_enrichment$p, method = "BH")


### add permutations to summary dataframe

df_enrichment_phyper_methyl_sensitive <-
  rbind(df_enrichment_phyper_methyl_sensitive, random_list_enrichment)

remove(random_list_enrichment)



##### plot adjusted p-value #####

# Extended Data Fig. 7a


### define color palettes
npg_colors <- pal_npg()(10)
show_col (npg_colors)


### plot
# plot -log10 of the adjusted p-value
# plot line for cutoff for significance (adjusted p 0.05)
# plot separate plots for diffferentially expressed and methylation sensitive genes

### fill colour defined by gene list for which enrichment was tested
# Immune genes: npg [1] (red), CTAs: npg [2] (light blue), TSG: npg [7] (light green), oncogenes: npg [4] (dark blue)

plot <- ggarrange(
  # differentially expressed genes
  ggplot(data = df_enrichment_phyper_DEG, aes(x = list, y = -log10(p_adjust), fill = list)) +
    geom_bar(stat = "identity", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    facet_wrap(~treatment, nrow = 1) +
    scale_fill_manual(values = c("Immune" = npg_colors[1], "CTA" = npg_colors[2],
                                  "TSG" = npg_colors[7], "onco" = npg_colors[4],
                                 "random" = "grey")) +
    scale_x_discrete(limits = c("Immune", "CTA", "TSG", "onco", "random")) +
    scale_y_continuous(limits = c(0,27)) +
    labs (title = "differentially expressed genes") +
    theme_bw(),
  # methylation sensitive genes
  ggplot(data = df_enrichment_phyper_methyl_sensitive, aes(x = list, y = -log10(p_adjust), fill = list)) +
    geom_bar(stat = "identity", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    facet_wrap(~treatment, nrow = 1) +
    scale_fill_manual(values = c("Immune" = npg_colors[1], "CTA" = npg_colors[2],
                                 "TSG" = npg_colors[7], "onco" = npg_colors[4],
                                 "random" = "grey")) +
    scale_x_discrete(limits = c("Immune", "CTA", "TSG", "onco", "random")) +
    scale_y_continuous(limits = c(0,27)) +
    labs (title = "methylation sensitive genes") +
    theme_bw(),
  nrow = 2, ncol = 1, common.legend = TRUE
)

# export plot
pdf("R_exports/plots/Fig_S7a_DMR_p_value_enrichment_DEG_methylSensitive_bar.pdf")
plot(plot)
dev.off()

remove(plot)