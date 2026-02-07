library(Seurat)
options(Seurat.object.assay.version = 'v5')
library(Signac)
library(ggplot2)
library(ggraph)
library(cowplot)
library(patchwork)
library(dplyr)
library(tidyverse)
library(dsb)
library(GenomicRanges)
library(future)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(tidyverse)
library(GenomeInfoDb)
library(TFBSTools)
library(JASPAR2022)
library(harmony)
library(clustree)
library(scSHC)
library(stringr)
library(scDblFinder)
library(BiocParallel)
library(loomR)
library(pheatmap)
library(SeuratDisk)
library(org.Mm.eg.db)
library(SeuratDisk)
library(data.table)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(loomR)
library(pheatmap)

##Figure 3A####
#require objects from 'merged tissue data preprocessing.R'
DimPlot(Combined_CD8, group.by = "orig.ident", reduction = "umap.wnn", label = TRUE, label.size = 6)

##Figure 3B####
#require objects from 'merged tissue data preprocessing.R'
DimPlot(Combined_CD8, group.by = "tissue_celltype", reduction = "umap.wnn", label = TRUE, label.size = 3) 

##Figure 3C####
#require objects from 'build GSEA.R'
DefaultAssay(Combined) <- "RNA.clean"; Idents(Combined) <- "tissue_celltype"
AVG.Combined.RNA <- AverageExpression(Combined, assays = "RNA.clean",return.seurat = TRUE, 
                                      group.by = 'tissue_celltype', layer = "data")
AVG.Combined.RNA$orig.ident <- colnames(AVG.Combined.RNA)

p1 <- DoHeatmap(AVG.Combined.RNA, features = LEG.REACTOME_TCR_SIGNALING, group.by = "orig.ident",
                disp.min = -2, disp.max = 2, draw.lines = F) + ggtitle("LEG.REACTOME_TCR_SIGNALING")
p2 <- DoHeatmap(AVG.Combined.RNA, features = LEG.PID_CD8_TCR_PATHWAY, group.by = "orig.ident",
                disp.min = -2, disp.max = 2, draw.lines = F) + ggtitle("LEG.PID_CD8_TCR_PATHWAY")
p3 <- DoHeatmap(AVG.Combined.RNA, features = LEG.PID_CD8_TCR_DOWNSTREAM_PATHWAY, group.by = "orig.ident",
                disp.min = -2, disp.max = 2, draw.lines = F) + ggtitle("LEG.PID_CD8_TCR_DOWNSTREAM_PATHWAY")
p4 <- DoHeatmap(AVG.Combined.RNA, features = LEG.HALLMARK_INFLAMMATORY_RESPONSE, group.by = "orig.ident",
                disp.min = -2, disp.max = 2, draw.lines = F) + ggtitle("LEG.HALLMARK_INFLAMMATORY_RESPONSE")
p5 <- DoHeatmap(AVG.Combined.RNA, features = LEG.HALLMARK_INTERFERON_GAMMA_RESPONSE, group.by = "orig.ident",
                disp.min = -2, disp.max = 2, draw.lines = F) + ggtitle("LEG.HALLMARK_INTERFERON_GAMMA_RESPONSE")
p6 <- DoHeatmap(AVG.Combined.RNA, features = LEG.PID_AP1_PATHWAY, group.by = "orig.ident",
                disp.min = -2, disp.max = 2, draw.lines = F) + ggtitle("LEG.PID_AP1_PATHWAY")
p7 <- DoHeatmap(AVG.Combined.RNA, features = LEG.HALLMARK_TNFA_SIGNALING_VIA_NFKB, group.by = "orig.ident",
                disp.min = -2, disp.max = 2, draw.lines = F) + ggtitle("LEG.HALLMARK_TNFA_SIGNALING_VIA_NFKB")

p11 <- p1[[1]]$data %>%  group_by(Cell) %>%  summarize(Mean = mean(Expression))
p12 <- p2[[1]]$data %>%  group_by(Cell) %>%  summarize(Mean = mean(Expression))
p13 <- p3[[1]]$data %>%  group_by(Cell) %>%  summarize(Mean = mean(Expression))
p14 <- p4[[1]]$data %>%  group_by(Cell) %>%  summarize(Mean = mean(Expression))
p15 <- p5[[1]]$data %>%  group_by(Cell) %>%  summarize(Mean = mean(Expression))
p16 <- p6[[1]]$data %>%  group_by(Cell) %>%  summarize(Mean = mean(Expression))
p17 <- p7[[1]]$data %>%  group_by(Cell) %>%  summarize(Mean = mean(Expression))

p <- cbind(p12, p11$Mean, p13$Mean, p14$Mean, p15$Mean, p16$Mean, p17$Mean)
rownames(p) <- p$Cell; p$Cell <- NULL
colnames(p) <- c("PID_CD8_TCR_PATHWAY", "REACTOME_TCR_SIGNALING", "PID_CD8_TCR_DOWNSTREAM_PATHWAY",
                 "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                 "PID_AP1_PATHWAY", "HALLMARK_TNFA_SIGNALING_VIA_NFKB")
p <- as.matrix(p); p <- t(p)
p <- as.data.frame(p)

p_order <- c("LN-Tscm", "LN-IFN-induced", "LN-Differentiating", "LN-Tem", "LN-Intermediate", 
             "LN-Terminal teff-exhausted", "LN-Proliferating",
             "spleen-Tscm", "spleen-Differentiating", "spleen-TEM", "spleen-Gzmk CD8", "spleen-Terminal teff-exhausted",
             "spleen-TRM regulatory", "spleen-gdT", "spleen-Proliferating",
             "Kidney-Stem-like I", "Kidney-Stem-like II", "Kidney-Activated", "Kidney-Innate-like", "Kidney-Trm",
             "Kidney-ISG", "Kidney-Cytotoxic I", "Kidney-Cytotoxic II",
             "Kidney-Terminally differentiated", "Kidney-Proliferating")
p <- p[, p_order]

ComplexHeatmap::Heatmap(p, row_names_gp = grid::gpar(fontsize = 10), cluster_columns = FALSE, cluster_rows = FALSE)

##Figure 3D####
#require objects from 'build GSEA.R'
DefaultAssay(Combined) <- "RNA.clean"; Idents(Combined) <- "tissue_celltype"
AVG.Combined.RNA <- AverageExpression(Combined, assays = "RNA.clean",return.seurat = TRUE, 
                                      group.by = 'tissue_celltype', layer = "data")
AVG.Combined.RNA$orig.ident <- colnames(AVG.Combined.RNA)

p_order <- c("LN-Tscm", "LN-IFN-induced", "LN-Differentiating", "LN-Tem", "LN-Intermediate", 
             "LN-Terminal teff-exhausted", "LN-Proliferating",
             "spleen-Tscm", "spleen-Differentiating", "spleen-TEM", "spleen-Gzmk CD8", "spleen-Terminal teff-exhausted",
             "spleen-TRM regulatory", "spleen-gdT", "spleen-Proliferating",
             "Kidney-Stem-like I", "Kidney-Stem-like II", "Kidney-Activated", "Kidney-Innate-like", "Kidney-Trm",
             "Kidney-ISG", "Kidney-Cytotoxic I", "Kidney-Cytotoxic II",
             "Kidney-Terminally differentiated", "Kidney-Proliferating")

DefaultAssay(AVG.Combined.RNA) <- "RNA.clean"; Idents(AVG.Combined.RNA) <- "orig.ident"
Idents(AVG.Combined.RNA) <- factor(Idents(AVG.Combined.RNA), levels = p_order)

DoHeatmap(AVG.Combined.RNA, features = LEG.REACTOME_TCR_SIGNALING,
          disp.min = -2, disp.max = 2, draw.lines = F) + ggtitle("LEG.REACTOME_TCR_SIGNALING")
DoHeatmap(AVG.Combined.RNA, features = LEG.PID_CD8_TCR_PATHWAY, 
          disp.min = -2, disp.max = 2, draw.lines = F) + ggtitle("LEG.PID_CD8_TCR_PATHWAY")
DoHeatmap(AVG.Combined.RNA, features = LEG.PID_CD8_TCR_DOWNSTREAM_PATHWAY, 
          disp.min = -2, disp.max = 2, draw.lines = F) + ggtitle("LEG.PID_CD8_TCR_DOWNSTREAM_PATHWAY")
DoHeatmap(AVG.Combined.RNA, features = LEG.HALLMARK_INFLAMMATORY_RESPONSE, 
          disp.min = -2, disp.max = 2, draw.lines = F) + ggtitle("LEG.HALLMARK_INFLAMMATORY_RESPONSE")
DoHeatmap(AVG.Combined.RNA, features = LEG.HALLMARK_INTERFERON_GAMMA_RESPONSE, 
          disp.min = -2, disp.max = 2, draw.lines = F) + ggtitle("LEG.HALLMARK_INTERFERON_GAMMA_RESPONSE")
DoHeatmap(AVG.Combined.RNA, features = LEG.PID_AP1_PATHWAY, 
          disp.min = -2, disp.max = 2, draw.lines = F) + ggtitle("LEG.PID_AP1_PATHWAY")
DoHeatmap(AVG.Combined.RNA, features = LEG.HALLMARK_TNFA_SIGNALING_VIA_NFKB, 
          disp.min = -2, disp.max = 2, draw.lines = F) + ggtitle("LEG.HALLMARK_TNFA_SIGNALING_VIA_NFKB")

##Figure 3E-F####
#require objects from 'build GRN.R'
exh.module <- c("Lag3","Itgav","E130308A19Rik", "Tank", "Ly75", "Rgs16",
                "Havcr2","Trerf1","Ikzf2","Penk","St6galnac3","Tox",
                "Trps1","Maf","Neb","Itga4","Gpr65","Gm13684",
                "Kif13b","Entpd1","Ccl5","Aopep","Ier3","Tnip3",
                "Nr4a2","Tmcc3","Cdk14","Wdfy2","Themis","Ttn",
                "Ptprk","Rabgap1l","Rasgef1b","Spire1","Osr2","Casp3",
                "Ctla4","Grhl1","Klhdc10")   

activation.module <- c("Id2","Cxcr6","Zfp36l2", "Slamf7","Ccl4","Junb","Actb","Tmsb4x","Ms4a4b",
                       "Ccl3","Dusp1","Klf6","Ifng","Gem","Jun","Fos","Hspa1a","Hspa1b",
                       "Hsp90ab1", "S100a4","Lgals1","Atf7ip","Fosb","Dnajb1","Ddx60","S100a6","Actg1",
                       "Fchsd2","Adgre5","Dnajc1","Vim","Gpr137c","Pcca")    

corKidney.CD8_corSpleen.CD8 <- list(corKidney.CD8, corSpleen.CD8)
splitcorheatmap(corKidney.CD8_corSpleen.CD8, exh.module, 
                "exh_module_Kidney_Spleen") 
splitcorheatmap(corKidney.CD8_corSpleen.CD8, activation.module, 
                "act_module_Kidney_Spleen") 

corKidney.CD8_corLN.CD8 <- list(corKidney.CD8, corLN.CD8)
splitcorheatmap(corKidney.CD8_corLN.CD8, exh.module, 
                "exh_module_Kidney_LN") 
splitcorheatmap(corKidney.CD8_corLN.CD8, activation.module, 
                "act_module_Kidney_LN") 

corSpleen.CD8_corLN.CD8 <- list(corSpleen.CD8, corLN.CD8)
splitcorheatmap(corSpleen.CD8_corLN.CD8, exh.module, 
                "exh_module_Spleen_LN") 
splitcorheatmap(corSpleen.CD8_corLN.CD8, activation.module, 
                "act_module_Spleen_LN") 

##Figure 3G####
#require objects from 'build GRN.R'
rankDrivers(figR.d.CD8,rankBy = "meanScore",interactive = FALSE)

##Figure 3H####
#require objects from 'build GRN.R'
#determine Runx3-driven/upregulated genes
Runx3_pos_edges <- figR.d.CD8 %>%
  mutate(
    Motif = trimws(Motif),                              # guard against stray spaces
    Enrichment.Q = p.adjust(Enrichment.P, "BH"),
    Corr.Q       = p.adjust(Corr.P,       "BH")
  ) %>%
  dplyr::filter(
    Motif == "Runx3",                                     # exact uppercase, no "::"
    Enrichment.Z > 0,              # motif evidence (FDR)
    Corr > 0, Corr.P < 0.05,
    Score > 0                 # positive regulation (FDR)
  ) %>%
  arrange(Corr.Q, desc(Corr))

Runx3_pos.reg_genes <- Runx3_pos_edges %>% distinct(DORC) %>% pull(DORC)

DefaultAssay(Kidney.CD8) <- "RNA.clean"
Kidney.CD8 <- AddModuleScore(Kidney.CD8, features = list(Runx3_pos.reg_genes), nbin = 10,
                             name = 'Runx3_pos.reg_genes')

#featureplot of Runx3 TF accessibility
DefaultAssay(Kidney.CD8) <- "chromvar"
FeaturePlot(Kidney.CD8, feature = c("MA0684.2"),  reduction = "umap.wnn", pt.size = 0.5,
            min.cutoff = 'q10', max.cutoff = 'q90')
#featureplot of Runx3-driven gene module score
FeaturePlot(Kidney.CD8, feature = c("Runx3_pos.reg_genes1"),  reduction = "umap.wnn", pt.size = 0.5,
            min.cutoff = 'q10', max.cutoff = 'q90')
