# Supplemental Figure 7b

# Author: Sandra Hellberg
# Email: sandra.hellberg@liu.se

# Date: 2024-02-13
# Last modified: 2024-08-16

# Description: 

rm(list=ls()) # remove all entries in the global environment 
set.seed(311) # seed for reproducibility

#-------------------------------------------------------------------------------

              ### Set directory structure and load packages---- 

#-------------------------------------------------------------------------------
setwd("Set your working directory here")
rna_dir <- "input directory for RNA-seq data"
fig_dir <- "output directory for plots"

lib_path <- .libPaths("C:/Rpackages/") # specifing the package library
#update.packages(lib.loc = lib_path, ask = FALSE, checkBuilt = TRUE, dependencies = TRUE)
#BiocManager::install(lib.loc = lib_path, update = TRUE, ask = FALSE)


pack_R <- c("TissueEnrich", "tidyr", "ggpubr", "stringr",
            "openxlsx", "dplyr") # libraries to load

for (i in 1:length(pack_R)) {
  library(pack_R[i], character.only = TRUE)
}

#-------------------------------------------------------------------------------

                  ### Read in and format input data---- 

#-------------------------------------------------------------------------------
gene_input <- read.delim(paste0(rna_dir,"Supplementary_Table_12.csv"), sep = "," )
loucy_input <- gene_input[grepl("L", gene_input$sample), ] # subset LOUCY
supt1_input <- gene_input[grepl("S", gene_input$sample), ] # subset SUPT1

background_loucy <- gene_input[grepl("LG7", gene_input$sample), ] # select one of the samples for background, should be same for all samples
background_loucy <- as.character(background_loucy$gene_id)

background_supt1 <- gene_input[grepl("SG7", gene_input$sample), ]
background_supt1 <- as.character(background_supt1$gene_id)

filter_rna <- function(data, sample_id) {
  data %>%
    filter(grepl(sample_id, sample),        
           logFC > 1, 
           padjust < 0.05) %>%
    pull('gene_id') 
}

### Combine TEs
# LOUCY
LD3 <- filter_rna(loucy_input, "LD3")
LG3 <- filter_rna(loucy_input, "LG3")
LG7 <- filter_rna(loucy_input, "LG7")

# SUPT1
SD3 <- filter_rna(supt1_input, "SD3")
SG3 <- filter_rna(supt1_input, "SG3")
SG7 <- filter_rna(supt1_input, "SG7")



#-------------------------------------------------------------------------------                      

                ### Tissue enrichment with TissueEnrich----

#-------------------------------------------------------------------------------

## LOUCY----

# background
gs_background_L <- background_loucy
gs_background_L <- GeneSet(geneIds=gs_background_L,organism="Homo Sapiens",geneIdType=SymbolIdentifier())

# LD3
gs_LD3 <- GeneSet(geneIds=LD3,organism="Homo Sapiens",geneIdType=SymbolIdentifier())

LD3_tissue <- teEnrichment(inputGenes = gs_LD3 , rnaSeqDataset = 1,
                           tissueSpecificGeneType = 2, multiHypoCorrection = TRUE,
                           backgroundGenes = gs_background_L)

LD3_seEnrichmentOutput <- LD3_tissue[[1]]
LD3_enrichmentOutput <- setNames(data.frame(assay(LD3_seEnrichmentOutput),row.names = rowData(LD3_seEnrichmentOutput)[,1]), colData(LD3_seEnrichmentOutput)[,1])
LD3_enrichmentOutput$Tissue<-row.names(LD3_enrichmentOutput)
head(LD3_enrichmentOutput)

LD3_genes <- LD3_enrichmentOutput

LD3_genes <- LD3_genes[!LD3_genes$Tissue.Specific.Genes == 0,]
LD3_genes <- LD3_genes[!LD3_genes$Log10PValue == 0,] # remove tissues where pvalue=1, i.e. -log10Pvalue=0
LD3_genes$Sample <- "LD3" 


LD3_tissue_plot <- ggplot(LD3_genes,aes(x=Tissue, y=Log10PValue,label = Tissue.Specific.Genes,fill = Tissue))+
  geom_bar(stat = 'identity')+
  labs(x='', y = '-LOG10(P-Adjusted)')+
  ylim(0, 32)+
  theme_bw()+
  theme(legend.position="none")+
  theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),panel.grid.major= element_blank(),panel.grid.minor = element_blank()) +
  geom_hline(yintercept = -log10(0.05))

# Not significantly enriched!

# LG3
gs_LG3 <- GeneSet(geneIds=LG3,organism="Homo Sapiens",geneIdType=SymbolIdentifier())

LG3_tissue <- teEnrichment(inputGenes = gs_LG3 , rnaSeqDataset = 1,
                           tissueSpecificGeneType = 2, multiHypoCorrection = TRUE,
                           backgroundGenes = gs_background_L)

LG3_seEnrichmentOutput <- LG3_tissue[[1]]
LG3_enrichmentOutput <- setNames(data.frame(assay(LG3_seEnrichmentOutput),row.names = rowData(LG3_seEnrichmentOutput)[,1]), colData(LG3_seEnrichmentOutput)[,1])
LG3_enrichmentOutput$Tissue<-row.names(LG3_enrichmentOutput)
head(LG3_enrichmentOutput)

LG3_genes <- LG3_enrichmentOutput

LG3_genes <- LG3_genes[!LG3_genes$Tissue.Specific.Genes == 0,]
LG3_genes <- LG3_genes[!LG3_genes$Log10PValue == 0,] # remove tissues where pvalue=1, i.e. -log10Pvalue=0
LG3_genes$Sample <- "LG3" 

LG3_tissue_plot <- ggplot(LG3_genes,aes(x=Tissue, y=Log10PValue,label = Tissue.Specific.Genes,fill = Tissue))+
  geom_bar(stat = 'identity')+
  labs(x='', y = '-LOG10(P-Adjusted)')+
  ylim(0, 32)+ 
  theme_bw()+
  theme(legend.position="none")+
  theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),panel.grid.major= element_blank(),panel.grid.minor = element_blank()) +
  geom_hline(yintercept = -log10(0.05))

# Retrieve the genes that are significantly enriched
LG3_genes_sig <- LG3_genes[LG3_genes$Log10PValue > -log10(0.05),]
print(LG3_genes_sig$Tissue) # "Testis" 

# Testis
LG3_Testis <- LG3_tissue[[3]][["Testis"]]
LG3_Testis <- data.frame(assay(LG3_Testis))
LG3_Testis$Sample <- "LG3"
LG3_Testis$Tissue <- "Testis"


# LG7
gs_LG7 <- GeneSet(geneIds=LG7, organism="Homo Sapiens",geneIdType=SymbolIdentifier())

LG7_tissue <- teEnrichment(inputGenes = gs_LG7 , rnaSeqDataset = 1,
                           tissueSpecificGeneType = 2, multiHypoCorrection = TRUE,
                           backgroundGenes = gs_background_L)

LG7_seEnrichmentOutput <- LG7_tissue[[1]]
LG7_enrichmentOutput <- setNames(data.frame(assay(LG7_seEnrichmentOutput),row.names = rowData(LG7_seEnrichmentOutput)[,1]), colData(LG7_seEnrichmentOutput)[,1])
LG7_enrichmentOutput$Tissue<-row.names(LG7_enrichmentOutput)
head(LG7_enrichmentOutput)

LG7_genes <- LG7_enrichmentOutput

LG7_genes <- LG7_genes[!LG7_genes$Tissue.Specific.Genes == 0,]
LG7_genes <- LG7_genes[!LG7_genes$Log10PValue == 0,] # remove tissues where pvalue=1, i.e. -log10Pvalue=0
LG7_genes$Sample <- "LG7" 


LG7_tissue_plot <- ggplot(LG7_genes,aes(x=Tissue, y=Log10PValue,label = Tissue.Specific.Genes,fill = Tissue))+
  geom_bar(stat = 'identity')+
  labs(x='', y = '-LOG10(P-Adjusted)')+
  ylim(0, 50)+ 
  theme_bw()+
  theme(legend.position="none")+
  theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),panel.grid.major= element_blank(),panel.grid.minor = element_blank()) +
  geom_hline(yintercept = -log10(0.05))

# Retrieve the genes that are significantly enriched
LG7_genes_sig <- LG7_genes[LG7_genes$Log10PValue > -log10(0.05),]
print(LG7_genes_sig$Tissue) # "Cerebral Cortex" "Kidney"          "Placenta"        "Skeletal Muscle" "Testis" 

# Cerebral cortex
LG7_CC <- LG7_tissue[[3]][["Cerebral Cortex"]]
LG7_CC <- data.frame(assay(LG7_CC))
LG7_CC$Sample <- "LG7"
LG7_CC$Tissue <- "Cerebral Cortex"

# Kidney
LG7_K <- LG7_tissue[[3]][["Kidney"]]
LG7_K <- data.frame(assay(LG7_K))
LG7_K$Sample <- "LG7"
LG7_K$Tissue <- "Kidney"

# Placenta
LG7_P <- LG7_tissue[[3]][["Placenta"]]
LG7_P <- data.frame(assay(LG7_P))
LG7_P$Sample <- "LG7"
LG7_P$Tissue <- "Placenta"

# Skeletal Muscle
LG7_SM <- LG7_tissue[[3]][["Skeletal Muscle"]]
LG7_SM <- data.frame(assay(LG7_SM))
LG7_SM$Sample <- "LG7"
LG7_SM$Tissue <- "Skeletal Muscle"

# Testis
LG7_Testis <- LG7_tissue[[3]][["Testis"]]
LG7_Testis <- data.frame(assay(LG7_Testis))
LG7_Testis$Sample <- "LG7"
LG7_Testis$Tissue <- "Testis"

Loucy_tissue <- data.frame(rbind(LG3_Testis,LG7_CC, LG7_K, LG7_P, LG7_SM, LG7_Testis)) 

# Combine all data and make one plot
LOUCY_tissue_combined <- rbind(LG3_genes, LG7_genes)
LOUCY_tissue_combined$Tissue <- rownames(LOUCY_tissue_combined)

ordered_tissues <- LOUCY_tissue_combined$Tissue # Adjust this list according to your desired order
LOUCY_tissue_combined$Tissue <- factor(LOUCY_tissue_combined$Tissue, levels = ordered_tissues)

LOUCY_tissue_combined_plot <- ggplot(LOUCY_tissue_combined,aes(x=Tissue, y=Log10PValue,label = Tissue.Specific.Genes,fill = Tissue))+
  geom_bar(stat = 'identity')+
  labs(x='', y = '-LOG10(P-Adjusted)')+
  ylim(0, 50)+ 
  theme_bw()+
  theme(legend.position="none")+
  theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),panel.grid.major= element_blank(),panel.grid.minor = element_blank()) +
  geom_hline(yintercept = -log10(0.05))


## SUP-T1----

# background
gs_background_S <- background_supt1
gs_background_S <- GeneSet(geneIds=gs_background_S,organism="Homo Sapiens",geneIdType=SymbolIdentifier())

# SG3
gs_SG3 <- GeneSet(geneIds=SG3,organism="Homo Sapiens",geneIdType=SymbolIdentifier())

SG3_tissue <- teEnrichment(inputGenes = gs_SG3 , rnaSeqDataset = 1,
                           tissueSpecificGeneType = 2, multiHypoCorrection = TRUE,
                           backgroundGenes = gs_background_S)

SG3_seEnrichmentOutput <- SG3_tissue[[1]]
SG3_enrichmentOutput <- setNames(data.frame(assay(SG3_seEnrichmentOutput),row.names = rowData(SG3_seEnrichmentOutput)[,1]), colData(SG3_seEnrichmentOutput)[,1])
SG3_enrichmentOutput$Tissue<-row.names(SG3_enrichmentOutput)
head(SG3_enrichmentOutput)

SG3_genes <- SG3_enrichmentOutput

SG3_genes <- SG3_genes[!SG3_genes$Tissue.Specific.Genes == 0,] # remove tissues with no enriched genes
SG3_genes <- SG3_genes[!SG3_genes$Log10PValue == 0,] # remove tissues where pvalue=1, i.e. -log10Pvalue=0
SG3_genes$Sample <- "SG3" 

SG3_tissue_plot <- ggplot(SG3_genes,aes(x=Tissue, y=Log10PValue,label = Tissue.Specific.Genes,fill = Tissue))+
  geom_bar(stat = 'identity')+
  labs(x='', y = '-LOG10(P-Adjusted)')+
  ylim(0, 50)+
  theme_bw()+
  theme(legend.position="none")+
  theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),panel.grid.major= element_blank(),panel.grid.minor = element_blank()) +
  geom_hline(yintercept = -log10(0.05))

# Retrieve the genes that are significantly enriched
SG3_genes_sig <- SG3_genes[SG3_genes$Log10PValue > -log10(0.05),]
print(SG3_genes_sig$Tissue) # "Testis" 

# Testis
SG3_Testis <- SG3_tissue[[3]][["Testis"]]
SG3_Testis <- data.frame(assay(SG3_Testis))
SG3_Testis$Sample <- "SG3"
SG3_Testis$Tissue <- "Testis"


# SG7
gs_SG7 <- GeneSet(geneIds=SG7, organism="Homo Sapiens",geneIdType=SymbolIdentifier())

SG7_tissue <- teEnrichment(inputGenes = gs_SG7 , rnaSeqDataset = 1,
                           tissueSpecificGeneType = 2, multiHypoCorrection = TRUE,
                           backgroundGenes = gs_background_S)

SG7_seEnrichmentOutput <- SG7_tissue[[1]]
SG7_enrichmentOutput <- setNames(data.frame(assay(SG7_seEnrichmentOutput),row.names = rowData(SG7_seEnrichmentOutput)[,1]), colData(SG7_seEnrichmentOutput)[,1])
SG7_enrichmentOutput$Tissue<-row.names(SG7_enrichmentOutput)
head(SG7_enrichmentOutput)

SG7_genes <- SG7_enrichmentOutput

SG7_genes <- SG7_genes[!SG7_genes$Tissue.Specific.Genes == 0,] # remove tissues with no enriched genes
SG7_genes <- SG7_genes[!SG7_genes$Log10PValue == 0,] # remove tissues where pvalue=1, i.e. -log10Pvalue=0
SG7_genes$Sample <- "SG7" 

SG7_tissue_plot <- ggplot(SG7_genes,aes(x=Tissue, y=Log10PValue,label = Tissue.Specific.Genes,fill = Tissue))+
  geom_bar(stat = 'identity')+
  labs(x='', y = '-LOG10(P-Adjusted)')+
  ylim(0, 50)+
  theme_bw()+
  theme(legend.position="none")+
  theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),panel.grid.major= element_blank(),panel.grid.minor = element_blank()) +
  geom_hline(yintercept = -log10(0.05))

# Retrieve the genes that are significantly enriched
SG7_genes_sig <- SG7_genes[SG7_genes$Log10PValue > -log10(0.05),]
print(SG7_genes_sig$Tissue) # "Placenta" "Testis" 

# Placenta
SG7_Placenta <- SG7_tissue[[3]][["Placenta"]]
SG7_Placenta <- data.frame(assay(SG7_Placenta))
SG7_Placenta$Sample <- "SG7"
SG7_Placenta$Tissue <- "Placenta"

# Testis
SG7_Testis <- SG7_tissue[[3]][["Testis"]]
SG7_Testis <- data.frame(assay(SG7_Testis))
SG7_Testis$Sample <- "SG7"
SG7_Testis$Tissue <- "Testis"

SUPT1_tissue <- data.frame(rbind(SG3_Testis, SG7_Placenta, SG7_Testis)) 

# Combine all data and make one plot
SUPT1_tissue_combined <- rbind(SG3_genes, SG7_genes)
SUPT1_tissue_combined$Tissue <- rownames(SUPT1_tissue_combined)

ordered_tissues <- SUPT1_tissue_combined$Tissue # Adjust this list according to your desired order
SUPT1_tissue_combined$Tissue <- factor(SUPT1_tissue_combined$Tissue, levels = ordered_tissues)

SUPT1_tissue_combined_plot <- ggplot(SUPT1_tissue_combined,aes(x=Tissue, y=Log10PValue,label = Tissue.Specific.Genes,fill = Tissue))+
  geom_bar(stat = 'identity')+
  labs(x='', y = '-LOG10(P-Adjusted)')+
  ylim(0, 50)+ 
  theme_bw()+
  theme(legend.position="none")+
  theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),panel.grid.major= element_blank(),panel.grid.minor = element_blank()) +
  geom_hline(yintercept = -log10(0.05))


# output combined plot
plot <- ggarrange(LOUCY_tissue_combined_plot, SUPT1_tissue_combined_plot, nrow=1, ncol=2)

pdf(paste0(fig_dir,"Tissue_enrichment_DEGs.pdf"),
    width = 20,
    height = 10)
plot(plot)
dev.off()


### END OF SCRIPT


