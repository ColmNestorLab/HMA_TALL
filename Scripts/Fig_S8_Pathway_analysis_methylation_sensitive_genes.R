# Supplemental Figure S8

# Author: Sandra Hellberg
# Email: sandra.hellberg@liu.se

# Date: 2024-02-13
# Last modified: 2024-08-21

rm(list=ls()) # remove all entries in the global environment 
set.seed(311) # seed for reproducibility

#-------------------------------------------------------------------------------

### Set directory structure and load packages---- 

#-------------------------------------------------------------------------------
setwd("Working directoru")
rna_dir <- "Data directory"
fig_dir <- "Output folder for figures"


lib_path <- .libPaths("C:/Rpackages/") # specifing the package library
#update.packages(lib.loc = lib_path, ask = FALSE, checkBuilt = TRUE, dependencies = TRUE)
#BiocManager::install(lib.loc = lib_path, update = TRUE, ask = FALSE)


pack_R <- c("TissueEnrich", "tidyr", "ggpubr", "stringr", "clusterProfiler")


for (i in 1:length(pack_R)) {
  library(pack_R[i], character.only = TRUE)
}

#-------------------------------------------------------------------------------

                  ### Read in and format input data---- 

#-------------------------------------------------------------------------------
gene_input <- read.delim(paste0(rna_dir,"Supplementary_Table_13.csv"), sep = "," )

# LOUCY
LOUCY <- gene_input[grep("LG7", gene_input$sample),]
LOUCY_MSG <- unique(LOUCY$gene_id[LOUCY$methyl_sensitive == "yes"]) # select methylation-sensitive genes

LOUCY_background <- LOUCY[LOUCY$meth.diff <= -25 & LOUCY$qvalue <0.01,]
LOUCY_background <- unique(LOUCY_background$gene_id)

# SUPT1
SUPT1 <- gene_input[grep("SG7", gene_input$sample),]
SUPT1_MSG <- unique(SUPT1$gene_id[SUPT1$methyl_sensitive == "yes"]) # select methylation-sensitive genes

SUPT1_background <- SUPT1[SUPT1$meth.diff <= -25 & SUPT1$qvalue <0.01,]
SUPT1_background <- unique(SUPT1_background$gene_id)


#-------------------------------------------------------------------------------                      

                    ### KEGG/GO pathway analysis----

#-------------------------------------------------------------------------------

# convert symbol to entrez
LOUCY_background <- bitr(LOUCY_background, fromType = "SYMBOL",toType = "ENTREZID", OrgDb = "org.Hs.eg.db" ) # genes to be used as background
LOUCY_MSG <- bitr(LOUCY_MSG, fromType = "SYMBOL",toType = "ENTREZID", OrgDb = "org.Hs.eg.db" ) 

SUPT1_background <- bitr(SUPT1_background, fromType = "SYMBOL",toType = "ENTREZID", OrgDb = "org.Hs.eg.db" ) # genes to be used as background
SUPT1_MSG <- bitr(SUPT1_MSG, fromType = "SYMBOL",toType = "ENTREZID", OrgDb = "org.Hs.eg.db" ) 

# LOUCY
LOUCY_kegg <- enrichKEGG(
  as.character(LOUCY_MSG[,2]),
  organism = "hsa",
  keyType = "kegg",
  universe = LOUCY_background$ENTREZID,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  pAdjustMethod = "BH")

dotplot(LOUCY_kegg, showCategory=20, x="count")
LOUCY_kegg <- setReadable(LOUCY_kegg, OrgDb = "org.Hs.eg.db", keyType="ENTREZID")
LOUCY_kegg <- LOUCY_kegg@result

L_go <- enrichGO(gene          = as.character(LOUCY_MSG[,2]),
                 universe      = LOUCY_background$ENTREZID,
                 OrgDb         = "org.Hs.eg.db",
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)

results_L_go <- L_go@result

L_go_plot <- dotplot(L_go, showCategory=15)

# SUPT1
SUPT1_kegg <- enrichKEGG(
  as.character(SUPT1_MSG[,2]),
  organism = "hsa",
  keyType = "kegg",
  universe = SUPT1_background$ENTREZID,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  pAdjustMethod = "BH")

dotplot(SUPT1_kegg, showCategory=20, x="count")
SUPT1_kegg <- setReadable(SUPT1_kegg, OrgDb = "org.Hs.eg.db", keyType="ENTREZID")
SUPT1_kegg <- SUPT1_kegg@result

S_go <- enrichGO(gene          = as.character(SUPT1_MSG[,2]),
                 universe      = SUPT1_background$ENTREZID,
                 OrgDb         = "org.Hs.eg.db",
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)

results_S_go <- S_go@result
S_go_plot <- dotplot(S_go, showCategory=15)


#-------------------------------------------------------------------------------                      

                          ### Export data----

#-------------------------------------------------------------------------------

plot <- ggarrange(L_go_plot,S_go_plot, ncol=2, nrow=1) # arrange plots
pdf(
  paste0(fig_dir, "GO_MSG_240903.pdf"),
  width = 15,
  height = 5
)
plot(plot)
dev.off()

names <- list('L_KEGG' = LOUCY_kegg, 'L_GO' = results_L_go, 'S_KEGG' = SUPT1_kegg, 'S_GO' = results_S_go)
write.xlsx(names, file = 'Pathway_analysis_MSG_240903.xlsx')

# End of script
