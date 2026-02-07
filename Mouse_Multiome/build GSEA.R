library(Seurat)
options(Seurat.object.assay.version = 'v5')
library(Signac)
library(ggplot2)
library(ggraph)
library(ggrepel)
library(cowplot)
library(patchwork)
library(dplyr)
library(tidyverse)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(tidyverse)
library(TFBSTools)
library(dplyr)
library(Matrix)
library(fgsea)
library(msigdbr)
library(tibble)
library(reshape2)
library(scales)
library(clusterProfiler)
library(org.Mm.eg.db)
library("AnnotationHub")
library(ComplexHeatmap)

##require objects built from "merged tissue data processing.R"
##this script apply fgsea and msigdbr to compute GSEA using RNA assay.

##make GSEA sets####
m_df_H<- msigdbr(species = "Mus musculus", category = "H")
m_df_H<- rbind(msigdbr(species = "Mus musculus", category = "C2"), m_df_H)
m_df_H<- rbind(msigdbr(species = "Mus musculus", category = "C7"), m_df_H)
fgsea_sets<- m_df_H %>% split(x = .$gene_symbol, f = .$gs_name)

##identify GSEA for each CD8 subset####
#repeat for other subsets that are shared across tissues, such as TEM, TEFF, gdT, and proliferating
DefaultAssay(Combined) <- "RNA.clean"; Idents(Combined) <- "tissue_celltype"

Tscm <- subset(Combined, idents = c("Kidney_Tcm", "Kidney_Tscm", "LN_Tscm", "spleen_Tscm"))
Tscm <- NormalizeData(Tscm, normalization.method = 'LogNormalize')
Tscm <- ScaleData(Tscm, verbose = TRUE)

#identify GSEA in kidney vs rLN and spleen
Idents(Tscm) <- "orig.ident"; DefaultAssay(Tscm) <- "RNA.clean"
DEG.GSEA.Tscm <- FindMarkers(Tscm, ident.1 = c("Kidney"),
                             only.pos = FALSE, logfc.threshold = 0, min.pct = 0, test.use = "wilcox")
Rank.DEG.GSEA.Tscm <- DEG.GSEA.Tscm$avg_log2FC
names(Rank.DEG.GSEA.Tscm) <- rownames(DEG.GSEA.Tscm)
fgsea_DEG.GSEA.Tscm <- fgsea(pathways = fgsea_sets, 
                             stats = Rank.DEG.GSEA.Tscm,
                             eps   = 0.0, minSize=15, maxSize=500, scoreType = "pos")

#after repeating for each cell, we generate fgsea for each major celltype and for combined data
#first make a gene list of significantly enriched gene sets of interest, for example:
gene_set_list <- c("REACTOME_TCR_SIGNALING",
                   "PID_CD8_TCR_PATHWAY",
                   "PID_CD8_TCR_DOWNSTREAM_PATHWAY",
                   "HALLMARK_INFLAMMATORY_RESPONSE",
                   "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                   "PID_AP1_PATHWAY",
                   "HALLMARK_TNFA_SIGNALING_VIA_NFKB")

filtered_fgsea_DEG.GSEA.Combined <- fgsea_DEG.GSEA.Combined[fgsea_DEG.GSEA.Combined$pathway %in% gene_set_list, ]
filtered_fgsea_DEG.GSEA.Tscm <- fgsea_DEG.GSEA.Tscm[fgsea_DEG.GSEA.Tscm$pathway %in% gene_set_list, ]
filtered_fgsea_DEG.GSEA.Differentiating <- fgsea_DEG.GSEA.Differentiating[fgsea_DEG.GSEA.Differentiating$pathway %in% gene_set_list, ]
filtered_fgsea_DEG.GSEA.gdT <- fgsea_DEG.GSEA.gdT[fgsea_DEG.GSEA.gdT$pathway %in% gene_set_list, ]
filtered_fgsea_DEG.GSEA.Proliferating <- fgsea_DEG.GSEA.Proliferating[fgsea_DEG.GSEA.Proliferating$pathway %in% gene_set_list, ]
filtered_fgsea_DEG.GSEA.TEFF <- fgsea_DEG.GSEA.TEFF[fgsea_DEG.GSEA.TEFF$pathway %in% gene_set_list, ]
filtered_fgsea_DEG.GSEA.TEM <- fgsea_DEG.GSEA.TEM[fgsea_DEG.GSEA.TEM$pathway %in% gene_set_list, ]

filtered_fgsea_DEG.GSEA.Combined$celltype <- "Combined"
filtered_fgsea_DEG.GSEA.Tscm$celltype <- "Tscm"
filtered_fgsea_DEG.GSEA.Differentiating$celltype <- "Differentiating"
filtered_fgsea_DEG.GSEA.gdT$celltype <- "gdT"
filtered_fgsea_DEG.GSEA.Proliferating$celltype <- "Proliferating"
filtered_fgsea_DEG.GSEA.TEFF$celltype <- "TEFF"
filtered_fgsea_DEG.GSEA.TEM$celltype <- "TEM"

filtered_fgsea <- rbind(filtered_fgsea_DEG.GSEA.Combined,
                        filtered_fgsea_DEG.GSEA.Tscm,
                        filtered_fgsea_DEG.GSEA.Differentiating,
                        filtered_fgsea_DEG.GSEA.gdT,
                        filtered_fgsea_DEG.GSEA.Proliferating,
                        filtered_fgsea_DEG.GSEA.TEFF,
                        filtered_fgsea_DEG.GSEA.TEM)

filtered_fgsea$log_adj <- -log10(filtered_fgsea$padj)




##identify leading edge genes from GSEA sets of interest####
gene_set_list <- c("REACTOME_TCR_SIGNALING",
                   "PID_CD8_TCR_PATHWAY",
                   "PID_CD8_TCR_DOWNSTREAM_PATHWAY",
                   "HALLMARK_INFLAMMATORY_RESPONSE",
                   "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                   "PID_AP1_PATHWAY",
                   "HALLMARK_TNFA_SIGNALING_VIA_NFKB")

LEG.REACTOME_TCR_SIGNALING <- filtered_fgsea_DEG.GSEA.Combined[which(filtered_fgsea_DEG.GSEA.Combined$pathway == 
                                                                       'REACTOME_TCR_SIGNALING'),]$leadingEdge
LEG.PID_CD8_TCR_PATHWAY <- filtered_fgsea_DEG.GSEA.Combined[which(filtered_fgsea_DEG.GSEA.Combined$pathway == 
                                                                    'PID_CD8_TCR_PATHWAY'),]$leadingEdge
LEG.PID_CD8_TCR_DOWNSTREAM_PATHWAY <- filtered_fgsea_DEG.GSEA.Combined[which(filtered_fgsea_DEG.GSEA.Combined$pathway == 
                                                                               'PID_CD8_TCR_DOWNSTREAM_PATHWAY'),]$leadingEdge
LEG.HALLMARK_INFLAMMATORY_RESPONSE <- filtered_fgsea_DEG.GSEA.Combined[which(filtered_fgsea_DEG.GSEA.Combined$pathway == 
                                                                               'HALLMARK_INFLAMMATORY_RESPONSE'),]$leadingEdge
LEG.HALLMARK_INTERFERON_GAMMA_RESPONSE <- filtered_fgsea_DEG.GSEA.Combined[which(filtered_fgsea_DEG.GSEA.Combined$pathway == 
                                                                                   'HALLMARK_INTERFERON_GAMMA_RESPONSE'),]$leadingEdge
LEG.PID_AP1_PATHWAY <- filtered_fgsea_DEG.GSEA.Combined[which(filtered_fgsea_DEG.GSEA.Combined$pathway == 
                                                                'PID_AP1_PATHWAY'),]$leadingEdge
LEG.HALLMARK_TNFA_SIGNALING_VIA_NFKB <- filtered_fgsea_DEG.GSEA.Combined[which(filtered_fgsea_DEG.GSEA.Combined$pathway == 
                                                                                 'HALLMARK_TNFA_SIGNALING_VIA_NFKB'),]$leadingEdge

LEG.REACTOME_TCR_SIGNALING <- unlist(LEG.REACTOME_TCR_SIGNALING)
LEG.PID_CD8_TCR_PATHWAY <- unlist(LEG.PID_CD8_TCR_PATHWAY)
LEG.PID_CD8_TCR_DOWNSTREAM_PATHWAY <- unlist(LEG.PID_CD8_TCR_DOWNSTREAM_PATHWAY)
LEG.HALLMARK_INFLAMMATORY_RESPONSE <- unlist(LEG.HALLMARK_INFLAMMATORY_RESPONSE)
LEG.HALLMARK_INTERFERON_GAMMA_RESPONSE <- unlist(LEG.HALLMARK_INTERFERON_GAMMA_RESPONSE)
LEG.PID_AP1_PATHWAY <- unlist(LEG.PID_AP1_PATHWAY)
LEG.HALLMARK_TNFA_SIGNALING_VIA_NFKB <- unlist(LEG.HALLMARK_TNFA_SIGNALING_VIA_NFKB)

