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

##this script describe steps to build and process 10x Multiome RNA and ATAC assays for each tissue (kidney, spleen, rLN) using Cellranger-ARC outputs

##build RNA Seurat objects from raw feature_bc_matrix####
#kidney
Kidney.raw <- Read10X("~/Kidney_MAM_cellranger/raw_feature_bc_matrix")
Kidney <- CreateSeuratObject(counts = Kidney.raw$`Gene Expression`, assay =  "RNA")
Kidney$orig.ident <- 'Kidney'
Kidney[["percent.mt"]] <- PercentageFeatureSet(Kidney, pattern = "^mt-")
#spleen
spleen.raw <- Read10X("~/Spleen_Combined/outs/raw_feature_bc_matrix")
spleen <- CreateSeuratObject(counts = spleen.raw$`Gene Expression`, assay =  "RNA")
spleen$orig.ident <- 'spleen'
spleen[["percent.mt"]] <- PercentageFeatureSet(spleen, pattern = "^mt-")
#rLN
LN.raw <- Read10X("~/LN_MAM_cellranger/raw_feature_bc_matrix")
LN <- CreateSeuratObject(counts = LN.raw$`Gene Expression`, assay =  "RNA")
LN$orig.ident <- 'LN'
LN[["percent.mt"]] <- PercentageFeatureSet(LN, pattern = "^mt-")
#RNA QC filter 
kidney <- subset(kidney, nCount_RNA > 500 & nFeature_RNA > 240 & percent.mt < 7.5)
spleen <- subset(spleen, nCount_RNA > 500 & nFeature_RNA > 240 & percent.mt < 7.5)
LN <- subset(LN, nCount_RNA > 500 & nFeature_RNA > 240 & percent.mt < 7.5)

##create ATAC Seurat Objects from peaks and fragments####
#kidney
Kidney.peaks <- read.table(file = "~/Kidney_MAM_cellranger/atac_peaks.bed", col.names = c("chr", "start", "end"))
Kidney.peaks.gr <- makeGRangesFromDataFrame(Kidney.peaks)
Kidney.md <- read.table(file = "~/Kidney_MAM_cellranger/per_barcode_metrics.csv",
  stringsAsFactors = FALSE, sep = ",", header = TRUE, row.names = 1)[-1, ]
Kidney.frags <- CreateFragmentObject(path = "~/Kidney_MAM_cellranger/atac_fragments.tsv.gz", cells = rownames(Kidney.md))
Kidney.counts <- FeatureMatrix( fragments = Kidney.frags, features = Kidney.peaks.gr, cells = rownames(Kidney.md))
Kidney.assay <- CreateChromatinAssay(Kidney.counts, fragments = Kidney.frags, genome = "mm10", sep = c(":", "-"))
Kidney.ATAC <- CreateSeuratObject(Kidney.assay, assay = "ATAC", meta.data= Kidney.md)
Kidney.ATAC$orig.ident <- 'Kidney'

#spleen
spleen.peaks <- read.table(file = "~/Spleen_Combined/outs/atac_peaks.bed", col.names = c("chr", "start", "end"))
spleen.peaks.gr <- makeGRangesFromDataFrame(spleen.peaks)
spleen.md <- read.table(file = "~/Spleen_Combined/outs/per_barcode_metrics.csv",
  stringsAsFactors = FALSE, sep = ",", header = TRUE, row.names = 1)[-1, ]
spleen.frags <- CreateFragmentObject(path = "~/Spleen_Combined/outs/atac_fragments.tsv.gz", cells = rownames(spleen.md))
spleen.counts <- FeatureMatrix(fragments = spleen.frags, features = spleen.peaks.gr, cells = rownames(spleen.md))
spleen.assay <- CreateChromatinAssay(spleen.counts, fragments = spleen.frags, genome = "mm10", sep = c(":", "-"))
spleen.ATAC <- CreateSeuratObject(spleen.assay, assay = "ATAC", meta.data= spleen.md)
spleen.ATAC$orig.ident <- 'spleen'

#rLN
LN.peaks <- read.table(file = "~/LN_MAM_cellranger/atac_peaks.bed", col.names = c("chr", "start", "end"))
LN.peaks.gr <- makeGRangesFromDataFrame(LN.peaks)
LN.md <- read.table(file = "~/LN_MAM_cellranger/per_barcode_metrics.csv",
  stringsAsFactors = FALSE, sep = ",", header = TRUE, row.names = 1)[-1, ]
LN.frags <- CreateFragmentObject(path = "~/LN_MAM_cellranger/atac_fragments.tsv.gz", cells = rownames(LN.md))
LN.counts <- FeatureMatrix(fragments = LN.frags, features = LN.peaks.gr, cells = rownames(LN.md))
LN.assay <- CreateChromatinAssay(LN.counts, fragments = LN.frags, genome = "mm10", sep = c(":", "-"))
LN.ATAC <- CreateSeuratObject(LN.assay, assay = "ATAC", meta.data= LN.md)
LN.ATAC$orig.ident <- 'LN'

#add annotations
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
genome(annotation) <- "mm10"

Annotation(Kidney.ATAC) <- annotation
Annotation(spleen.ATAC) <- annotation
Annotation(LN.ATAC) <- annotation

#ATAC QC filter
Kidney.ATAC <- NucleosomeSignal(Kidney.ATAC)
Kidney.ATAC$nucleosome_group <- ifelse(Kidney.ATAC$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = Kidney.ATAC, group.by = 'nucleosome_group', region = 'chr1-1-10000000')
Kidney.ATAC <- TSSEnrichment(Kidney.ATAC)
Kidney.ATAC$high.tss <- ifelse(Kidney.ATAC$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(Kidney.ATAC, group.by = 'high.tss') + NoLegend()
Kidney.ATAC <- subset(Kidney.ATAC, nucleosome_signal < 1 & TSS.enrichment > 2 & nCount_ATAC > 350 
                      & nCount_ATAC < 30000 & nFeature_ATAC > 200)

spleen.ATAC <- NucleosomeSignal(spleen.ATAC)
spleen.ATAC$nucleosome_group <- ifelse(spleen.ATAC$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = spleen.ATAC, group.by = 'nucleosome_group', region = 'chr1-1-10000000')
spleen.ATAC <- TSSEnrichment(spleen.ATAC)
spleen.ATAC$high.tss <- ifelse(spleen.ATAC$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(spleen.ATAC, group.by = 'high.tss') + NoLegend()
spleen.ATAC <- subset(spleen.ATAC, nucleosome_signal < 1 & TSS.enrichment > 2 & nCount_ATAC > 350 
                      & nCount_ATAC < 30000 & nFeature_ATAC > 200)

LN.ATAC <- NucleosomeSignal(LN.ATAC)
LN.ATAC$nucleosome_group <- ifelse(LN.ATAC$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = LN.ATAC, group.by = 'nucleosome_group', region = 'chr1-1-10000000')
LN.ATAC <- TSSEnrichment(LN.ATAC)
LN.ATAC$high.tss <- ifelse(LN.ATAC$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(LN.ATAC, group.by = 'high.tss') + NoLegend()
LN.ATAC <- subset(LN.ATAC, nucleosome_signal < 1 & TSS.enrichment > 2 & nCount_ATAC > 350 
                  & nCount_ATAC < 30000 & nFeature_ATAC > 200)

##create Seurat Object with merged RNA and ATAC####
#kidney
joint.bcs<- intersect(colnames(Kidney), colnames(Kidney.ATAC))
joint.bcs <- as.data.frame(joint.bcs)
Kidney.ATAC$join_bc <- ifelse(colnames(Kidney.ATAC) %in% joint.bcs$joint.bcs, "TRUE", "FALSE")
Kidney.ATAC <- subset(Kidney.ATAC, idents = c("TRUE"))
Kidney$join_bc <- ifelse(colnames(Kidney) %in% joint.bcs$joint.bcs, "TRUE", "FALSE")
Kidney <- subset(Kidney, idents = c("TRUE"))
Kidney[["ATAC"]] <- Kidney.ATAC[["ATAC"]] 

#spleen
joint.bcs<- intersect(colnames(spleen), colnames(spleen.ATAC))
joint.bcs <- as.data.frame(joint.bcs)
spleen.ATAC$join_bc <- ifelse(colnames(spleen.ATAC) %in% joint.bcs$joint.bcs, "TRUE", "FALSE")
spleen.ATAC <- subset(spleen.ATAC, idents = c("TRUE"))
spleen$join_bc <- ifelse(colnames(spleen) %in% joint.bcs$joint.bcs, "TRUE", "FALSE")
spleen <- subset(spleen, idents = c("TRUE"))
spleen[["ATAC"]] <- spleen.ATAC[["ATAC"]] 

#rLN
joint.bcs<- intersect(colnames(LN), colnames(LN.ATAC))
joint.bcs <- as.data.frame(joint.bcs)
LN.ATAC$join_bc <- ifelse(colnames(LN.ATAC) %in% joint.bcs$joint.bcs, "TRUE", "FALSE")
LN.ATAC <- subset(LN.ATAC, idents = c("TRUE"))
LN$join_bc <- ifelse(colnames(LN) %in% joint.bcs$joint.bcs, "TRUE", "FALSE")
LN <- subset(LN, idents = c("TRUE"))
LN[["ATAC"]] <- LN.ATAC[["ATAC"]] 

##remove doublets by scdblfinder####
#kidney
sce_Kidney <- Kidney@assays$RNA@layers$counts
sce_Kidney <- scDblFinder(sce_Kidney, BPPARAM = bp)
Kidney$scDblFInder <- sce_Kidney$scDblFinder.class
Kidney <- subset(Kidney, idents = c("singlet"))
#spleen
sce_spleen <- spleen@assays$RNA@layers$counts
sce_spleen <- scDblFinder(sce_spleen, BPPARAM = bp)
spleen$scDblFInder <- sce_spleen$scDblFinder.class
spleen <- subset(spleen, idents = c("singlet"))
#rLN
sce_LN <- LN@assays$RNA@layers$counts
sce_LN <- scDblFinder(sce_LN, BPPARAM = bp)
LN$scDblFInder <- sce_LN$scDblFinder.class
LN <- subset(LN, idents = c("singlet"))

##MACS2 peak recall####
#kidney
DefaultAssay(Kidney) <- "ATAC"
MACS2peaks <- CallPeaks(Kidney)
MACS2peaks <- keepStandardChromosomes(MACS2peaks, pruning.mode = "coarse")
MACS2_counts <- FeatureMatrix(fragments = Fragments(Kidney), features = MACS2peaks, cells = colnames(Kidney))
Kidney[["MACS2ATAC"]] <- CreateChromatinAssay(counts = MACS2_counts, fragments = Fragments(Kidney), annotation = annotation)
#spleen
DefaultAssay(spleen) <- "ATAC"
MACS2peaks <- CallPeaks(spleen)
MACS2peaks <- keepStandardChromosomes(MACS2peaks, pruning.mode = "coarse")
MACS2_counts <- FeatureMatrix(fragments = Fragments(spleen), features = MACS2peaks, cells = colnames(spleen))
spleen[["MACS2ATAC"]] <- CreateChromatinAssay(counts = MACS2_counts, fragments = Fragments(spleen), annotation = annotation)
#rLN
DefaultAssay(LN) <- "ATAC"
MACS2peaks <- CallPeaks(LN)
MACS2peaks <- keepStandardChromosomes(MACS2peaks, pruning.mode = "coarse")
MACS2_counts <- FeatureMatrix(fragments = Fragments(LN), features = MACS2peaks, cells = colnames(LN))
LN[["MACS2ATAC"]] <- CreateChromatinAssay(counts = MACS2_counts, fragments = Fragments(LN), annotation = annotation)

##MACS2 ATAC normalization and generate ATAC_UMAP####
#kidney
DefaultAssay(Kidney) <- "MACS2ATAC"
Kidney <- RunTFIDF(Kidney)
Kidney <- FindTopFeatures(Kidney, min.cutoff = 'q0')
Kidney <- RunSVD(Kidney)
Kidney <- RunUMAP(Kidney, dims = 2:30, reduction = 'lsi', reduction.name = "umap.atac", reduction.key = "atacUMAP_", seed.use = 42)
#spleen
DefaultAssay(spleen) <- "MACS2ATAC"
spleen <- RunTFIDF(spleen)
spleen <- FindTopFeatures(spleen, min.cutoff = 'q0')
spleen <- RunSVD(spleen)
spleen <- RunUMAP(spleen, dims = 2:20, reduction = 'lsi', reduction.name = "umap.atac", reduction.key = "atacUMAP_", seed.use = 42)
#rLN
DefaultAssay(LN) <- "MACS2ATAC"
LN <- RunTFIDF(LN)
LN <- FindTopFeatures(LN, min.cutoff = 'q0')
LN <- RunSVD(LN)
LN <- RunUMAP(LN, dims = 2:20, reduction = 'lsi', reduction.name = "umap.atac", reduction.key = "atacUMAP_", seed.use = 42)

##RNA normalization and generate RNA_UMAP####
#kidney
DefaultAssay(Kidney) <- "RNA"
Kidney <- NormalizeData(Kidney, normalization.method = 'LogNormalize')
Kidney <- FindVariableFeatures(Kidney, selection.method = "vst", nfeatures = 2000)
Kidney <- ScaleData(Kidney, features = rownames(Kidney), verbose = TRUE)
Kidney <- RunPCA(Kidney, features = VariableFeatures(Kidney), verbose = TRUE, 
                 reduction.name = 'pca', assay = 'RNA')
Kidney <- RunUMAP(Kidney, reduction = "pca", seed.use = 42,  dims = 1:35,
                  reduction.name = "umap.rna", reduction.key = "rnaUMAP_")
#spleen
DefaultAssay(spleen) <- "RNA"
spleen <- NormalizeData(spleen, normalization.method = 'LogNormalize')
spleen <- FindVariableFeatures(spleen, selection.method = "vst", nfeatures = 2000)
spleen <- ScaleData(spleen, features = rownames(spleen), verbose = TRUE)
spleen <- RunPCA(spleen, features = VariableFeatures(spleen), verbose = TRUE, 
                 reduction.name = 'pca', assay = 'RNA')
spleen <- RunUMAP(spleen, reduction = "pca", seed.use = 42,  dims = 1:40,
                  reduction.name = "umap.rna", reduction.key = "rnaUMAP_")
#rLN
DefaultAssay(LN) <- "RNA"
LN <- NormalizeData(LN, normalization.method = 'LogNormalize')
LN <- FindVariableFeatures(LN, selection.method = "vst", nfeatures = 2000)
LN <- ScaleData(LN, features = rownames(LN), verbose = TRUE)
LN <- RunPCA(LN, features = VariableFeatures(LN), verbose = TRUE, 
             reduction.name = 'pca', assay = 'RNA')
LN <- RunUMAP(LN, reduction = "pca", seed.use = 42,  dims = 1:50,
              reduction.name = "umap.rna", reduction.key = "rnaUMAP_")

##generate WNN UMAP and subset out CD4 T and CD8 T####
#kidney
Kidney <- FindMultiModalNeighbors(Kidney, reduction.list = list("pca", "lsi"), k.nn = 10,
                                  dims.list = list(1:35, 2:30), modality.weight.name = "WNN.weight")
Kidney <- RunUMAP(Kidney, nn.name = "weighted.nn", reduction.name = "umap.wnn", reduction.key = "wnnUMAP_", seed.use = 42)
Kidney <- FindClusters(Kidney, graph.name = "wsnn", algorithm = 3, random.seed = 42, resolution = c(0.5), verbose = TRUE)
#spleen
spleen <- FindMultiModalNeighbors(spleen, reduction.list = list("pca", "lsi"), k.nn = 10,
                                  dims.list = list(1:40, 2:20), modality.weight.name = "WNN.weight")
spleen <- RunUMAP(spleen, nn.name = "weighted.nn", reduction.name = "umap.wnn", reduction.key = "wnnUMAP_", seed.use = 42)
spleen <- FindClusters(spleen, graph.name = "wsnn", algorithm = 3, random.seed = 42, resolution = c(15), verbose = TRUE)
#rLN
LN <- FindMultiModalNeighbors(LN, reduction.list = list("pca", "lsi"), k.nn = 10,
                              dims.list = list(1:50, 2:20), modality.weight.name = "WNN.weight")
LN <- RunUMAP(LN, nn.name = "weighted.nn", reduction.name = "umap.wnn", reduction.key = "wnnUMAP_", seed.use = 42)
LN <- FindClusters(LN, graph.name = "wsnn", algorithm = 3, random.seed = 42, resolution = c(0.5), verbose = TRUE)

#on cluster annotations:
#note that the UMAP algorithm is stochastic; therefore, results (including cluster structures and numbering) will vary between runs.
#As a result, a fully reproducible record of CD4 and CD8 T cell cluster annotation steps is not provided
#some clusters were further refined by subclustering, based on expression of key CD4 vs CD8 T cell defining markers.
#an example of identifying CD4 T vs CD8 T in kidney, based on Cd3e, Cd4, and Cd8a RNA expression is provided below.
#the same workflow is applied to spleen and rLN datasets to create spleen.CD4, spleen.CD8, LN.CD4, and LN.CD8 objects
Idents(Kidney) <- "seurat_clusters"
Kidney <- FindSubCluster(Kidney, cluster = "10", graph.name = "wsnn", subcluster.name = "sub.seurat_clusters", resolution = 1, algorithm = 3)
DimPlot(Kidney, reduction = "umap.wnn", group.by = "sub.seurat_clusters", label = TRUE, label.size = 3, raster=FALSE)
Idents(Kidney) <- "sub.seurat_clusters"
Kidney <-RenameIdents(Kidney, '10_2' = '10_0')
Kidney$seurat_clusters <- Idents(Kidney)

Kidney$major_celltype <- Kidney$seurat_clusters
Idents(Kidney) <- "major_celltype"
Kidney <-RenameIdents(Kidney, '11' = 'CD8', '1' = 'CD8', '4' = 'CD8', '8' = 'CD8', '7' = 'CD8', '10_0' = 'CD8',
                      '3' = 'CD8', '9' = 'CD8', '13' = 'other',
                      '10_1' = 'CD4', '5' = 'CD4', '0' = 'CD4', '2' = 'CD4', '6' = 'CD4', '12' = 'CD4')
Kidney$major_celltype <- Idents(Kidney)

Kidney.CD4 <- subset(Kidney, idents = c("CD4"))
Kidney.CD8 <- subset(Kidney, idents = c("CD8"))

##re-normalize CD4 T and CD8 T cell subsets####
#the same steps as below were repeated for kidney.CD8, and for spleen and rLN datasets:
RNA.matrix <- Kidney.CD4[["RNA"]]$counts
RNA.matrix <- RNA.matrix[!grepl('mt', rownames(RNA.matrix)),]
RNA.matrix <- RNA.matrix[!grepl('Igh', rownames(RNA.matrix)),]
RNA.matrix <- RNA.matrix[!grepl('Igk', rownames(RNA.matrix)),]
RNA.matrix <- RNA.matrix[!grepl('Igl', rownames(RNA.matrix)),]
RNA.matrix <- RNA.matrix[!grepl('Jchain', rownames(RNA.matrix)),]
RNA.matrix <- RNA.matrix[!grepl('Rpl', rownames(RNA.matrix)),]
RNA.matrix <- RNA.matrix[!grepl('Rps', rownames(RNA.matrix)),]
Kidney.CD4[["RNA.clean"]] <- CreateAssayObject(counts = RNA.matrix )
#normalize RNA assays
DefaultAssay(Kidney.CD4) <- "RNA.clean"
Kidney.CD4 <- NormalizeData(Kidney.CD4, normalization.method = 'LogNormalize')
Kidney.CD4 <- FindVariableFeatures(Kidney.CD4, selection.method = "vst", nfeatures = 2000)
Kidney.CD4 <- ScaleData(Kidney.CD4, verbose = TRUE)
#normalize MACS2_ATAC assays
DefaultAssay(Kidney.CD4) <- "MACS2ATAC"
Kidney.CD4 <- RunTFIDF(Kidney.CD4)
Kidney.CD4 <- FindTopFeatures(Kidney.CD4, min.cutoff = 'q0')
Kidney.CD4 <- RunSVD(Kidney.CD4)

##build gene accessibility assay from ATAC####
#the same steps as below were repeated for kidney.CD8, and for spleen and rLN datasets:
DefaultAssay(Kidney.CD4) <- "ATAC"
gene.activities <- GeneActivity(Kidney.CD4)
Kidney.CD4[['activities']] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(Kidney.CD4) <- "activities"
Kidney.CD4 <- NormalizeData(object = Kidney.CD4, assay = 'activities', 
                            normalization.method = 'LogNormalize', 
                            scale.factor = median(Kidney.CD4$nCount_activities))
Kidney.CD4 <- ScaleData(Kidney.CD4)

##add JASPAR2022 TF motifs and chromVAR scores####
#the same steps as below were repeated for kidney.CD8, and for spleen and rLN datasets:
pfm <- getMatrixSet(x = JASPAR2022, opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))
main.chroms <- standardChromosomes(BSgenome.Mmusculus.UCSC.mm10)

DefaultAssay(Kidney.CD4) <- "MACS2ATAC"
Kidney.CD4 <- AddMotifs(Kidney.CD4, genome = BSgenome.Mmusculus.UCSC.mm10, pfm = pfm)
motif.names <- Kidney.CD8@assays$MACS2ATAC@motifs@motif.names
Kidney.CD4 <- RunChromVAR(Kidney.CD4, genome = BSgenome.Mmusculus.UCSC.mm10)

##generate WNN UMAP and identify clusters for CD4 T and CD8 T cell subsets####
#kidney
DefaultAssay(Kidney.CD4) <- "MACS2ATAC"
Kidney.CD4 <- RunUMAP(Kidney.CD4, dims = 2:20, reduction = 'lsi',
                      reduction.name = "umap.atac", reduction.key = "atacUMAP_", seed.use = 42)
DefaultAssay(Kidney.CD4) <- "RNA"
Kidney.CD4 <- RunUMAP(Kidney.CD4, reduction = "pca", seed.use = 42,  dims = 1:35,
                      reduction.name = "umap.rna", reduction.key = "rnaUMAP_")
DefaultAssay(Kidney.CD8) <- "MACS2ATAC"
Kidney.CD8 <- RunUMAP(Kidney.CD8, dims = 2:10, reduction = 'lsi',
                      reduction.name = "umap.atac", reduction.key = "atacUMAP_", seed.use = 42)
DefaultAssay(Kidney.CD8) <- "RNA"
Kidney.CD8 <- RunUMAP(Kidney.CD8, reduction = "pca", seed.use = 42,  dims = 1:40,
                      reduction.name = "umap.rna", reduction.key = "rnaUMAP_")

Kidney.CD4 <- FindMultiModalNeighbors(Kidney.CD4, reduction.list = list("pca", "lsi"), k.nn = 10,
                                      dims.list = list(1:35, 2:20), modality.weight.name = "WNN.weight")
Kidney.CD4 <- RunUMAP(Kidney.CD4, nn.name = "weighted.nn", reduction.name = "umap.wnn", reduction.key = "wnnUMAP_", seed.use = 42)
Kidney.CD4 <- FindClusters(Kidney.CD4, graph.name = "wsnn", algorithm = 3, random.seed = 42, resolution = c(1.5), verbose = TRUE)

Kidney.CD8 <- FindMultiModalNeighbors(Kidney.CD8, reduction.list = list("pca", "lsi"), k.nn = 10,
                                      dims.list = list(1:35, 2:20), modality.weight.name = "WNN.weight")
Kidney.CD8 <- RunUMAP(Kidney.CD8, nn.name = "weighted.nn", reduction.name = "umap.wnn", reduction.key = "wnnUMAP_", seed.use = 42)
Kidney.CD8 <- FindClusters(Kidney.CD8, graph.name = "wsnn", algorithm = 3, random.seed = 42, resolution = c( 0.5), verbose = TRUE)

#spleen
DefaultAssay(spleen.CD4) <- "MACS2ATAC"
spleen.CD4 <- RunUMAP(spleen.CD4, dims = 2:20, reduction = 'lsi',
                      reduction.name = "umap.atac", reduction.key = "atacUMAP_", seed.use = 42)
DefaultAssay(spleen.CD4) <- "RNA"
spleen.CD4 <- RunUMAP(spleen.CD4, reduction = "pca", seed.use = 42,  dims = 1:40,
                      reduction.name = "umap.rna", reduction.key = "rnaUMAP_")
DefaultAssay(spleen.CD8) <- "MACS2ATAC"
spleen.CD8 <- RunUMAP(spleen.CD8, dims = 2:30, reduction = 'lsi',
                      reduction.name = "umap.atac", reduction.key = "atacUMAP_", seed.use = 42)
DefaultAssay(spleen.CD8) <- "RNA"
spleen.CD8 <- RunUMAP(spleen.CD8, reduction = "pca", seed.use = 42,  dims = 1:40,
                      reduction.name = "umap.rna", reduction.key = "rnaUMAP_")

spleen.CD4 <- FindMultiModalNeighbors(spleen.CD4, reduction.list = list("pca", "lsi"), k.nn = 10,
                                      dims.list = list(1:40, 2:20), modality.weight.name = "WNN.weight")
spleen.CD4 <- RunUMAP(spleen.CD4, nn.name = "weighted.nn", reduction.name = "umap.wnn", reduction.key = "wnnUMAP_", seed.use = 42)
spleen.CD4 <- FindClusters(spleen.CD4, graph.name = "wsnn", algorithm = 3, random.seed = 42, resolution = c(1), verbose = TRUE)

spleen.CD8 <- FindMultiModalNeighbors(spleen.CD8, reduction.list = list("pca", "lsi"), k.nn = 10,
                                      dims.list = list(1:40, 2:40), modality.weight.name = "WNN.weight")
spleen.CD8 <- RunUMAP(spleen.CD8, nn.name = "weighted.nn", 
                      reduction.name = "umap.wnn", reduction.key = "wnnUMAP_", seed.use = 42)
spleen.CD8 <- FindClusters(spleen.CD8, graph.name = "wsnn", algorithm = 3, random.seed = 42, resolution = c(1), verbose = TRUE)

#rLN
DefaultAssay(LN.CD4) <- "MACS2ATAC"
LN.CD4 <- RunUMAP(LN.CD4, dims = 2:20, reduction = 'lsi',
                  reduction.name = "umap.atac", reduction.key = "atacUMAP_", seed.use = 42)
DefaultAssay(LN.CD4) <- "RNA"
LN.CD4 <- RunUMAP(LN.CD4, reduction = "pca", seed.use = 42,  dims = 1:30,
                  reduction.name = "umap.rna", reduction.key = "rnaUMAP_")
DefaultAssay(LN.CD8) <- "MACS2ATAC"
LN.CD8 <- RunUMAP(LN.CD8, dims = 2:25, reduction = 'lsi',
                  reduction.name = "umap.atac", reduction.key = "atacUMAP_", seed.use = 42)
DefaultAssay(LN.CD8) <- "RNA"
LN.CD8 <- RunUMAP(LN.CD8, reduction = "pca", seed.use = 42,  dims = 1:40,
                  reduction.name = "umap.rna", reduction.key = "rnaUMAP_")

LN.CD4 <- FindMultiModalNeighbors(LN.CD4, reduction.list = list("pca", "lsi"), k.nn = 10,
                                  dims.list = list(1:30, 2:20), modality.weight.name = "WNN.weight")
LN.CD4 <- RunUMAP(LN.CD4, nn.name = "weighted.nn", 
                  reduction.name = "umap.wnn", reduction.key = "wnnUMAP_", seed.use = 42)
LN.CD4 <- FindClusters(LN.CD4, graph.name = "wsnn", algorithm = 3, random.seed = 42, resolution = c(1), verbose = TRUE)

LN.CD8 <- FindMultiModalNeighbors(LN.CD8, reduction.list = list("pca", "lsi"), k.nn = 10,
                                  dims.list = list(1:40, 2:25), modality.weight.name = "WNN.weight")
LN.CD8 <- RunUMAP(LN.CD8, nn.name = "weighted.nn",  reduction.name = "umap.wnn", reduction.key = "wnnUMAP_", seed.use = 42)
LN.CD8 <- FindClusters(LN.CD8, graph.name = "wsnn", algorithm = 3, random.seed = 42, resolution = c(1.2), verbose = TRUE)

#on cluster annotations:
#note that the UMAP algorithm is stochastic; therefore, results (including cluster structures and numbering) will vary between runs.
#As a result, a fully reproducible record of CD4 and CD8 T cell cluster annotation steps is not provided
#some clusters were further refined by subclustering, based on expression of key CD4 vs CD8 T cell defining markers.
#an example of identifying CD8 T cell subsets in kidney, based on RNA expression of celltype defining markers, is provided below.
#See Al Souz, et. al., Immunity 2026 for markers used.
#the same workflow is applied to spleen and rLN datasets to annotate spleen.CD4, spleen.CD8, LN.CD4, and LN.CD8 objects
Idents(Kidney.CD8) <- "seurat_clusters"
Kidney.CD8 <- FindSubCluster(Kidney.CD8, cluster = "5", graph.name = "wsnn", 
                             subcluster.name = "seurat_subclusters", resolution = 1, algorithm = 3)
DimPlot(Kidney.CD8, reduction = "umap.wnn", group.by = "seurat_subclusters", label = TRUE, label.size = 2.5)
Idents(Kidney.CD8) <- "seurat_subclusters"
Kidney.CD8 <-RenameIdents(Kidney.CD8, '5_4' = '12', '5_0' = '5', '5_1' = '5', '5_2' = '5', '5_3' = '5', 
                          '5_5' = '5', '5_6' = '5', '5_7' = '5', '5_8' = '5') 
Kidney.CD8$seurat_clusters <- Idents(Kidney.CD8)

Idents(Kidney.CD8) <- "seurat_clusters"
Kidney.CD8 <- FindSubCluster(Kidney.CD8, cluster = "0", graph.name = "wsnn", 
                             subcluster.name = "seurat_subclusters", resolution = 0.6, algorithm = 3)
Idents(Kidney.CD8) <- "seurat_subclusters"
Kidney.CD8 <-RenameIdents(Kidney.CD8, '0_1' = '11', '0_3' = '11', '0_5' = '3', 
                          '0_4' = '0', '0_2' = '0', '0_0' = '0', '0_6' = '0') 
Kidney.CD8$seurat_clusters <- Idents(Kidney.CD8)

Kidney.CD8$major_celltype <- Kidney.CD8$seurat_clusters
Idents(Kidney.CD8) <- "major_celltype"
Kidney.CD8 <-RenameIdents(Kidney.CD8, '8' = 'Proliferating', '9' = 'Tscm I', '0' = 'Trm', '12' = 'Trm', '11' = 'Innate-like',
                          '2' = 'Tscm II', '6' = 'Tem', '5' = 'Teff', '1' = 'Teff',
                          '3' = 'Activated', '7' = 'Activated', '10' = 'Activated', '4' = 'Terminal teff-exhausted')
Kidney.CD8$major_celltype <- Idents(Kidney.CD8)





