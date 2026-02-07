library(Seurat)
options(Seurat.object.assay.version = 'v5')
library(Signac)
library(ggplot2)
library(scales)
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
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(monocle3)

#these scripts apply monocle3 (https://cole-trapnell-lab.github.io/monocle-release/docs/) to build ATAC-based differentiation trajectory pseudotime.
#require processed objects from "per tissue data preprocessing.R"

##CD8 buid pseudotime####
#convert to as.cell_data_set
Kidney.CD8@reductions$UMAP <- Kidney.CD8@reductions$umap.wnn
Kidney.CD8.cds <- as.cell_data_set(Kidney.CD8)
Kidney.CD8.cds <- cluster_cells(cds = Kidney.CD8.cds, reduction_method = "UMAP", cluster_method = 'louvain')
Kidney.CD8.cds <- learn_graph(Kidney.CD8.cds, use_partition = TRUE)

#order cells
#note: here we define startpoint from Stem-like I  
Kidney.CD8.cds <- order_cells(Kidney.CD8.cds, reduction_method = "UMAP")

#pick and subset trajectory branch
#note: here we define 'starting node' from Stem-like I, and 'ending node' at Terminally differentiated, then press 'Connect nodes'
Kidney.CD8.cds_tscm.texh <- choose_graph_segments(Kidney.CD8.cds)

#add pseudotime score to seurat object
pData(Kidney.CD8.cds)$Pseudotime <- pseudotime(Kidney.CD8.cds)
Kidney.CD8$pseudotime <- pData(Kidney.CD8.cds)$Pseudotime

#subset seurat object by trajectory branch
Idents(Kidney.CD8) <- colnames(Kidney.CD8)
Kidney.CD8_tscm.texh <- subset(Kidney.CD8, idents = colnames(Kidney.CD8.cds_tscm.texh))

##differential feature expression by pseudotime bins from stem-like to terminally differentiated CD8 T in kidney####
#bin across 10 incremental group for selected trajectory branch 
Kidney.CD8_tscm.texh$Pseudotime_bin10 <- cut(Kidney.CD8_tscm.texh@meta.data$pseudotime, 10)

#identify differential TFs for each trajectory bins
DefaultAssay(Kidney.CD8_tscm.texh) <- "chromvar"; Idents(Kidney.CD8_tscm.texh) <- "Pseudotime_bin10"
Kidney.CD8_tscm.texh_bin.TFs <- FindAllMarkers(Kidney.CD8_tscm.texh, only.pos = TRUE, min.pct = 0,
                                               logfc.threshold = 0.3, 
                                               mean.fxn = rowMeans, fc.name = "avg_diff")
Kidney.CD8_tscm.texh_bin.TFs$BH <- p.adjust(Kidney.CD8_tscm.texh_bin.TFs$p_val, 
                                            method = "BH", n = nrow(Kidney.CD8_tscm.texh_bin.TFs))
Kidney.CD8_tscm.texh_bin.TFs$motifs <- ConvertMotifID(motif.names, id = rownames(Kidney.CD8_tscm.texh_bin.TFs))
Kidney.CD8_tscm.texh_bin.TFs <- na.omit(Kidney.CD8_tscm.texh_bin.TFs)
Kidney.CD8_tscm.texh_bin.TFs <- Kidney.CD8_tscm.texh_bin.TFs[!duplicated(Kidney.CD8_tscm.texh_bin.TFs$motifs),]
Kidney.CD8_tscm.texh_bin.TFs %>%
  group_by(cluster) %>%
  dplyr::filter(avg_diff > 0) %>%
  slice_head(n = 40) %>%
  ungroup() -> top10.Kidney.CD8_tscm.texh_bin.TFs

#identify differential RNA for each trajectory bins
DefaultAssay(Kidney.CD8_tscm.texh) <- "RNA.clean"; Idents(Kidney.CD8_tscm.texh) <- "Pseudotime_bin10"
Kidney.CD8_tscm.texh_bin.RNA  <- FindAllMarkers(Kidney.CD8_tscm.texh, only.pos = TRUE, min.pct = 0.3,
                                                logfc.threshold = 0.25)
Kidney.CD8_tscm.texh_bin.RNA$BH <- p.adjust(Kidney.CD8_tscm.texh_bin.RNA $p_val, 
                                            method = "BH", n = nrow(Kidney.CD8_tscm.texh_bin.RNA ))
Kidney.CD8_tscm.texh_bin.RNA %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0) %>%
  slice_head(n = 30) %>%
  ungroup() -> top10.Kidney.CD8_tscm.texh_bin.RNA
top10.Kidney.CD8_tscm.texh_bin.RNA <- top10.Kidney.CD8_tscm.texh_bin.RNA[!duplicated(top10.Kidney.CD8_tscm.texh_bin.RNA$gene),]

#identify differential gene accessibility for each trajectory bins
DefaultAssay(Kidney.CD8_tscm.texh) <- "activities"; Idents(Kidney.CD8_tscm.texh) <- "Pseudotime_bin10"
Kidney.CD8_tscm.texh_bin.activity <- FindAllMarkers(Kidney.CD8_tscm.texh, only.pos = TRUE, min.pct = 0.5,
                                                    logfc.threshold = 0.5, test.use = 'LR', latent.vars = 'nFeature_MACS2ATAC')
Kidney.CD8_tscm.texh_bin.activity$BH <- p.adjust(Kidney.CD8_tscm.texh_bin.activity $p_val, 
                                                 method = "BH", n = nrow(Kidney.CD8_tscm.texh_bin.activity ))
Kidney.CD8_tscm.texh_bin.activity <- Kidney.CD8_tscm.texh_bin.activity[!duplicated(Kidney.CD8_tscm.texh_bin.activity$gene), ]
top10.Kidney.CD8_tscm.texh_bin.activity <- Kidney.CD8_tscm.texh_bin.activity[
  Kidney.CD8_tscm.texh_bin.activity$gene %in% Kidney.CD8_tscm.texh_bin.RNA$gene, ]


