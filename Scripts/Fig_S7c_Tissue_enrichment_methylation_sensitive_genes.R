# Supplemental Figure 7c

# Author: Sandra Hellberg
# Email: sandra.hellberg@liu.se

# Date: 2024-02-13
# Last modified: 2024-08-21

rm(list=ls()) # remove all entries in the global environment 
set.seed(311) # seed for reproducibility

#-------------------------------------------------------------------------------

              ### Set directory structure and load packages---- 

#-------------------------------------------------------------------------------
setwd("Set your working directory here")
rna_dir <- "input directory"
fig_dir <- "output directory for figures"


lib_path <- .libPaths("C:/Rpackages/") # specifing the package library
#update.packages(lib.loc = lib_path, ask = FALSE, checkBuilt = TRUE, dependencies = TRUE)
#BiocManager::install(lib.loc = lib_path, update = TRUE, ask = FALSE)


pack_R <- c("TissueEnrich", "tidyr", "ggpubr", "stringr") # libraries to load

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

               ### Tissue enrichment with Tissue Enrich----

#-------------------------------------------------------------------------------

## LOUCY----
gs_background_L <- LOUCY_background
gs_background_L <- GeneSet(geneIds=gs_background_L,organism="Homo Sapiens",geneIdType=SymbolIdentifier())

gs_LOUCY <- GeneSet(geneIds=LOUCY_MSG,organism="Homo Sapiens",geneIdType=SymbolIdentifier())

LOUCY_tissue <- teEnrichment(inputGenes = gs_LOUCY , rnaSeqDataset = 1,
                           tissueSpecificGeneType = 2, multiHypoCorrection = TRUE,
                           backgroundGenes = gs_background_L)

LOUCY_seEnrichmentOutput <- LOUCY_tissue[[1]]
LOUCY_enrichmentOutput <- setNames(data.frame(assay(LOUCY_seEnrichmentOutput),row.names = rowData(LOUCY_seEnrichmentOutput)[,1]), colData(LOUCY_seEnrichmentOutput)[,1])
LOUCY_enrichmentOutput$Tissue<-row.names(LOUCY_enrichmentOutput)
head(LOUCY_enrichmentOutput)

LOUCY_genes <- LOUCY_enrichmentOutput

LOUCY_genes <- LOUCY_genes[!LOUCY_genes$Tissue.Specific.Genes == 0,]
LOUCY_genes <- LOUCY_genes[!LOUCY_genes$Log10PValue == 0,] # remove tissues where pvalue=1, i.e. -log10Pvalue=0
LOUCY_genes$Sample <- "LOUCY" 


LOUCY_tissue_plot <- ggplot(LOUCY_genes,aes(x=Tissue,-Log10PValue,y=Log10PValue,label = Tissue.Specific.Genes,fill = Tissue))+
  geom_bar(stat = 'identity')+
  labs(x='', y = '-LOG10(P-Adjusted)')+
  ylim(0, 10)+
  theme_bw()+
  theme(legend.position="none")+
  theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),panel.grid.major= element_blank(),panel.grid.minor = element_blank()) +
  geom_hline(yintercept = -log10(0.05))

# Retrieve the genes that are significantly enriched
LOUCY_genes_sig <- LOUCY_genes[LOUCY_genes$Log10PValue > -log10(0.05),]
print(LOUCY_genes_sig$Tissue) # "Testis" 

# Testis
LOUCY_Testis <- LOUCY_tissue[[3]][["Testis"]]
LOUCY_Testis <- data.frame(assay(LOUCY_Testis))
LOUCY_Testis$Sample <- "LOUCY"
LOUCY_Testis$Tissue <- "Testis"

## SUPT1----
gs_background_S <- SUPT1_background
gs_background_S <- GeneSet(geneIds=gs_background_S,organism="Homo Sapiens",geneIdType=SymbolIdentifier())

gs_SUPT1 <- GeneSet(geneIds=SUPT1_MSG,organism="Homo Sapiens",geneIdType=SymbolIdentifier())

SUPT1_tissue <- teEnrichment(inputGenes = gs_SUPT1 , rnaSeqDataset = 1,
                             tissueSpecificGeneType = 2, multiHypoCorrection = TRUE,
                             backgroundGenes = gs_background_S)

SUPT1_seEnrichmentOutput <- SUPT1_tissue[[1]]
SUPT1_enrichmentOutput <- setNames(data.frame(assay(SUPT1_seEnrichmentOutput),row.names = rowData(SUPT1_seEnrichmentOutput)[,1]), colData(SUPT1_seEnrichmentOutput)[,1])
SUPT1_enrichmentOutput$Tissue<-row.names(SUPT1_enrichmentOutput)
head(SUPT1_enrichmentOutput)

SUPT1_genes <- SUPT1_enrichmentOutput

SUPT1_genes <- SUPT1_genes[!SUPT1_genes$Tissue.Specific.Genes == 0,]
SUPT1_genes <- SUPT1_genes[!SUPT1_genes$Log10PValue == 0,] # remove tissues where pvalue=1, i.e. -log10Pvalue=0
SUPT1_genes$Sample <- "SUPT1" 


SUPT1_tissue_plot <- ggplot(SUPT1_genes,aes(x=Tissue,-Log10PValue,y=Log10PValue,label = Tissue.Specific.Genes,fill = Tissue))+
  geom_bar(stat = 'identity')+
  labs(x='', y = '-LOG10(P-Adjusted)')+
  ylim(0, 10)+
  theme_bw()+
  theme(legend.position="none")+
  theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),panel.grid.major= element_blank(),panel.grid.minor = element_blank()) +
  geom_hline(yintercept = -log10(0.05))

# Retrieve the genes that are significantly enriched
SUPT1_genes_sig <- SUPT1_genes[SUPT1_genes$Log10PValue > -log10(0.05),]
print(SUPT1_genes_sig$Tissue) # "Testis" 

# Testis
SUPT1_Testis <- SUPT1_tissue[[3]][["Testis"]]
SUPT1_Testis <- data.frame(assay(SUPT1_Testis))
SUPT1_Testis$Sample <- "SUPT1"
SUPT1_Testis$Tissue <- "Testis"

# output combined plot
plot <- ggarrange(LOUCY_tissue_plot, SUPT1_tissue_plot, nrow=1, ncol=2)

pdf(paste0(fig_dir,"Tissue_enrichment_MSGs.pdf"),
    width = 20,
    height = 10)
plot(plot)
dev.off()


# Script ends here

