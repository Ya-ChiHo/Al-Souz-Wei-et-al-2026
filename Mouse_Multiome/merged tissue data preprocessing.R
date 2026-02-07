library(Seurat)
options(Seurat.object.assay.version = 'v5')
library(Signac)
library(ggplot2)
library(ggraph)
library(cowplot)
library(patchwork)
library(dplyr)
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
library(scDblFinder)

##require objects built from "per tissue data processing.R"
##this script describe steps to merge Seurat objects from each tissue (kidney, spleen, rLN), for CD4 T,  CD8 T, and combined T cell subsets.

##reduce each object to counts only and merge CD8 T or CD4 T cells between tissues####
#repeat below to build Combined_CD4 dataset
LN <- CreateSeuratObject(counts = LN.CD8[["RNA"]]$counts, assay =  "RNA")
LN[["ATAC"]] <- LN.CD8[["ATAC"]]
LN$nCount_ATAC <- LN.CD8$nCount_ATAC
LN$nFeature_ATAC <- LN.CD8$nFeature_ATAC
LN$percent.mt <- LN.CD8$percent.mt
LN$orig.ident <- LN.CD8$orig.ident 
LN$major_celltype <- LN.CD8$major_celltype 

spleen <- CreateSeuratObject(counts = spleen.CD8[["RNA"]]$counts, assay =  "RNA")
spleen[["ATAC"]] <- spleen.CD8[["ATAC"]]
spleen$nCount_ATAC <- spleen.CD8$nCount_ATAC
spleen$nFeature_ATAC <- spleen.CD8$nFeature_ATAC
spleen$percent.mt <- spleen.CD8$percent.mt
spleen$orig.ident <- spleen.CD8$orig.ident 
spleen$major_celltype <- spleen.CD8$major_celltype 

kidney <- CreateSeuratObject(counts = Kidney.CD8[["RNA"]]$counts, assay =  "RNA")
kidney[["ATAC"]] <- Kidney.CD8[["ATAC"]]
kidney$nCount_ATAC <- Kidney.CD8$nCount_ATAC
kidney$nFeature_ATAC <- Kidney.CD8$nFeature_ATAC
kidney$percent.mt <- Kidney.CD8$percent.mt
kidney$orig.ident <- Kidney.CD8$orig.ident 
kidney$major_celltype <- Kidney.CD8$major_celltype 

Combined_CD8 <- merge(x = LN, y = c(spleen, kidney))
Combined_CD8 <- JoinLayers(Combined_CD8)

##re-normalize RNA assay####
#repeat below for Combined_CD4 dataset
DefaultAssay(Combined_CD8) <- "RNA"
Combined_CD8 <- JoinLayers(Combined_CD8)
RNA.matrix <- Combined_CD8[["RNA"]]$counts
RNA.matrix <- RNA.matrix[!grepl('mt', rownames(RNA.matrix)),]
RNA.matrix <- RNA.matrix[!grepl('Igh', rownames(RNA.matrix)),]
RNA.matrix <- RNA.matrix[!grepl('Igk', rownames(RNA.matrix)),]
RNA.matrix <- RNA.matrix[!grepl('Igl', rownames(RNA.matrix)),]
RNA.matrix <- RNA.matrix[!grepl('Jchain', rownames(RNA.matrix)),]
RNA.matrix <- RNA.matrix[!grepl('Rpl', rownames(RNA.matrix)),]
RNA.matrix <- RNA.matrix[!grepl('Rps', rownames(RNA.matrix)),]
Combined_CD8[["RNA.clean"]] <- CreateAssayObject(counts = RNA.matrix )

DefaultAssay(Combined_CD8) <- "RNA.clean"
Combined_CD8 <- NormalizeData(Combined_CD8, normalization.method = 'LogNormalize')
Combined_CD8 <- FindVariableFeatures(Combined_CD8, selection.method = "vst", nfeatures = 2000)
Combined_CD8 <- ScaleData(Combined_CD8, verbose = TRUE)

##repeat MACS2 peak recall, annotation, and peak normalization####
#repeat below for Combined_CD4 dataset
DefaultAssay(Combined_CD8) <- "ATAC"
Combined_CD8 <- RunTFIDF(Combined_CD8)
Combined_CD8 <- FindTopFeatures(Combined_CD8, min.cutoff = 'q0')
Combined_CD8 <- RunSVD(Combined_CD8)

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
genome(annotation) <- "mm10"
Annotation(Combined_CD8) <- annotation

MACS2peaks <- CallPeaks(Combined_CD8)
MACS2peaks <- keepStandardChromosomes(MACS2peaks, pruning.mode = "coarse")
MACS2_counts <- FeatureMatrix(fragments = Fragments(Combined_CD8), features = MACS2peaks, cells = colnames(Combined_CD8))
Combined_CD8[["MACS2ATAC"]] <- CreateChromatinAssay(counts = MACS2_counts, fragments = Fragments(Combined_CD8), annotation = annotation)

DefaultAssay(Combined_CD8) <- "MACS2ATAC"
Combined_CD8 <- RunTFIDF(Combined_CD8)
Combined_CD8 <- FindTopFeatures(Combined_CD8, min.cutoff = 'q0')
Combined_CD8 <- RunSVD(Combined_CD8)

##ATAC integration (batch effect correction),  normalization, and UMAP generation####
#repeat below for Combined_CD4 dataset
#batch-effect correction
Combined_CD8.atac_anchors <- FindIntegrationAnchors(object.list = SplitObject(Combined_CD8_ATAC_integration, split.by = "orig.ident"),
                                                    anchor.features = rownames(subset(Combined_CD8_ATAC_integration, idents = c("spleen"))), 
                                                    reduction = "rlsi", dims = 2:20)
Combined_CD8.atac_integration <- IntegrateEmbeddings(anchorset = Combined_CD8.atac_anchors, 
                                                     reductions = Combined_CD8_ATAC_integration[["lsi"]],
                                                     new.reduction.name = "integrated_lsi", 
                                                     dims.to.integrate = 1:20, k.weight = 30) 
Combined_CD8@reductions$integrated_lsi <- Combined_CD8.atac_integration@reductions$integrated_lsi

#run ATAC UMAP
DefaultAssay(Combined_CD8) <- "MACS2ATAC"
Combined_CD8 <- RunUMAP(Combined_CD8, reduction = "integrated_lsi", dims = 2:40, 
                        seed.use = 42, reduction.name = "umap.integrated.atac", reduction.key = "atacIntegratedUMAP_")

##RNA integration (batch effect correction),  normalization, and UMAP generation####
#repeat below for Combined_CD4 dataset
DefaultAssay(Combined_CD8) <- "RNA.clean"
Combined_CD8 <- NormalizeData(Combined_CD8, normalization.method = 'LogNormalize') %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>% RunPCA(npcs = 50) %>% RunHarmony("orig.ident", plot_convergence = TRUE) %>%
  RunUMAP(reduction = "harmony", reduction.name = "umap.integrated.rna", reduction.key = "rnaIntegratedUMAP_", 
          dims = 1:40, assay = 'RNA.clean', seed.use = 42) %>%
  identity()

##make WNN UMAP####
#repeat below for Combined_CD4 dataset
DefaultAssay(Combined_CD8) <- "RNA.clean"
Combined_CD8 <- FindMultiModalNeighbors(Combined_CD8, reduction.list = list("harmony", "integrated_lsi"), 
                                        dims.list = list(1:40, 2:20), modality.weight.name = "WNN.weight")
Combined_CD8 <- RunUMAP(Combined_CD8, nn.name = "weighted.nn", reduction.name = "umap.wnn", 
                        reduction.key = "wnnUMAP_", seed.use = 42)

#T cell subset annotation identity for individual tissues were used to label cells in the combined data: 
Combined_CD8$tissue_celltype <- str_c(Combined_CD8$orig.ident, '_', Combined_CD8$major_celltype)









