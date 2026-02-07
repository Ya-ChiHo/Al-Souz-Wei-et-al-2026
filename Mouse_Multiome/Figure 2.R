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
#library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
library(TFBSTools)
library(JASPAR2022)
library(harmony)
library(clustree)
library(scSHC)
library(stringr)
library(scDblFinder)
library(BiocParallel)
library(monocle3)
library(Matrix)
library(SeuratWrappers)
library(ArchR)
library(motifmatchr)
library(qlcMatrix)
library(WGCNA)
library(flashClust)
library(ade4)
library(glmGamPoi)
library(CelliD)
library(scriabin)
library(tidyverse)
library(ComplexHeatmap)
library(cowplot)
library(conflicted)
library(scales)
library(BuenColors)
library(FigR)
library(cisTopic)
library(parallel)
library(FNN)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(networkD3)
library(Hmisc)
library(qlcMatrix)
library(reshape2)
library(Matrix)
library(topGO)
library(fgsea)
library(msigdbr)
library(tibble)

##Figure 2A####
#require objects from 'per tissue data preprocessing.R'
DimPlot(Kidney.CD8, reduction = "umap.wnn", 
        label = TRUE, label.size = 10, pt.size = 1.5, repel = TRUE)

##Figure 2B####
#require objects from 'per tissue data preprocessing.R'
Kidney.CD8.RNA.markers <-c("Tcf7", "Slamf6", "Sell", "Lef1", "Ccr7", "Il7r", "Klf2",
                           "Ifng", "Tnf", "Jun", "Fos", 
                           "Trgv2", "Tcrg-C2", "Tcrg-C4", "Cxcr6",
                           "Il21r", "Il2rb", "Itga1",  "Itgae", 
                           "Il15ra","Ifit1", "Ifit3", "Isg15", "Irf7", "Klrk1", 
                           "Gzmb", "Prf1", "Il18r1", "Gzma", 
                           "Tox", "Entpd1", "Lag3", "Ctla4", "Pdcd1", "Havcr2", 
                           "Bag6", "Mki67", "Top2a")
DotPlot(Kidney.CD8, features = Kidney.CD8.RNA.markers, dot.scale = 30)

##Figure 2C####
#require objects from 'per tissue data preprocessing.R'
DefaultAssay(Kidney.CD8) <- "MACS2ATAC"; Idents(Kidney.CD8) <- "major_celltype"
Kidney.CD8 <- Footprint(Kidney.CD8, motif.name = c("TCF7"), assay = "MACS2ATAC", genome = BSgenome.Mmusculus.UCSC.mm10, in.peaks = TRUE)
PlotFootprint(Kidney.CD8, features = c("TCF7"))

##Figure 2D####
#require objects from 'build monocle2 trajectory.R'
Kidney.CD8$in_traj_tscm_Texh <- ifelse(colnames(Kidney.CD8) %in% colnames(Kidney.CD8_tscm.texh), TRUE, FALSE)
DimPlot(Kidney.CD8, reduction = "umap.wnn", group.by = "in_traj_tscm_Texh", label = TRUE, label.size = 4)

##Figure 2E####
#require objects from 'build monocle2 trajectory.R'
p <- top10.Kidney.CD8_tscm.texh_bin.TFs$gene
pt.matrix <- as.matrix(Kidney.CD8_tscm.texh[["chromvar"]]@data
                       [match(p, rownames(Kidney.CD8_tscm.texh[["chromvar"]]@data)),
                         order(Kidney.CD8_tscm.texh$pseudotime)])
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- ConvertMotifID(motif.names, id = rownames(pt.matrix))

ht <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  #km = 6,
  row_title_rot                = 0,
  cluster_rows                 = FALSE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  use_raster = TRUE)

##Figure 2F####
#require objects from 'build monocle2 trajectory.R'
p <- top10.Kidney.CD8_tscm.texh_bin.activity$gene

pt.matrix <- as.matrix(Kidney.CD8_tscm.texh[["activities"]]@data
                       [match(p, rownames(Kidney.CD8_tscm.texh[["activities"]]@data)),
                         order(Kidney.CD8_tscm.texh$pseudotime)])
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))

ht <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  #km = 6,
  row_title_rot                = 0,
  cluster_rows                 = FALSE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  use_raster = TRUE)

##Figure 2G####
#require objects from 'build monocle2 trajectory.R'
p <- top10.Kidney.CD8_tscm.texh_bin.RNA$gene

pt.matrix <- as.matrix(Kidney.CD8_tscm.texh[["RNA.clean"]]@data
                       [match(p, rownames(Kidney.CD8_tscm.texh[["RNA.clean"]]@data)),
                         order(Kidney.CD8_tscm.texh$pseudotime)])
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))

ht <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  #km = 6,
  row_title_rot                = 0,
  cluster_rows                 = FALSE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  use_raster = TRUE)

##Figure 2H####
#require objects from 'build monocle2 trajectory.R'
#for gene expression pt.matrix, see object made for Figure 2G 
#repeat for each gene
plot(pt.matrix[c("Tcf7"),], xlab = "cells ordered by pseudotime", ylab = "pseudotime", main = "Tcf7 gene expression dynamic change", col = "red") 

#for gene accessibility pt.matrix, see object made for Figure 2F 
#repeat for each gene
plot(pt.matrix[c("Tcf7"),], xlab = "cells ordered by pseudotime", ylab = "pseudotime", main = "Tcf7 activity expression dynamic change", col = "blue") 

##Figure 2I####
#require objects from 'build monocle2 trajectory.R'
#bin cells into 5 groups by pseudotime (early, early-mid, mid, mid-late, late) from stem-like to terminally differentiated states
Kidney.CD8_tscm.texh$Pseudotime_bin5 <- cut(Kidney.CD8_tscm.texh@meta.data$pseudotime, 5)
DefaultAssay(Kidney.CD8_tscm.texh) <- "MACS2ATAC"; Idents(Kidney.CD8_tscm.texh) <- "Pseudotime_bin5"

#identify significant peaks
da_peaks <- FindAllMarkers(object = Kidney.CD8_tscm.texh, only.pos = FALSE,
  test.use = 'LR', latent.vars = 'nCount_ATAC')
da_peaks$coord <- rownames(da_peaks)

#note: ranges.show is the significant peaks identified at each gene locus
#for example, for Tox: ranges.show <- StringToGRanges("chr4-6880614-6881595")
CoveragePlot(object = Kidney.CD8_tscm.texh,region = "Tcf7",
             ranges.title = "MACS2",extend.upstream = 2000,extend.downstream = 2000, ranges = ranges.show)
CoveragePlot(object = Kidney.CD8_tscm.texh, region = "Tox",
             ranges.title = "MACS2", extend.upstream = 2000, extend.downstream = 2000, ranges = ranges.show)
CoveragePlot(object = Kidney.CD8_tscm.texh, region = "Slamf6",
             ranges.title = "MACS2", extend.upstream = 2000, extend.downstream = 2000, ranges = ranges.show)
CoveragePlot(object = Kidney.CD8_tscm.texh, region = "Havcr2",
             ranges.title = "MACS2", extend.upstream = 2000, extend.downstream = 2000, ranges = ranges.show)
CoveragePlot(object = Kidney.CD8_tscm.texh, region = "Pdcd1",
             ranges.title = "MACS2", extend.upstream = 2000, extend.downstream = 2000, ranges = ranges.show)
