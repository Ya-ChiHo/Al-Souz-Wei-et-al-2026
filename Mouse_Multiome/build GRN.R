library(Seurat)
options(Seurat.object.assay.version = 'v5')
library(Signac)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(stringr)
library(cicero)
library(data.table)
library(Matrix)
library(GenomicRanges)
library(magrittr)
library(SummarizedExperiment)
library(optparse)
library(yaml)
library(Rcpp)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(chromVAR)
library(doParallel)
library(BuenColors)
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
library(scales)
library(Matrix)
library(tibble)


##require objects built from "per tissue data processing.R"
##this script apply cisTopic and FigR packages to generate the Gene Regulatory Networks using ATAC and RNA assays.
##shown are GRNs built for kidney CD8 T cells
##the same steps were repearted for spleen and rLN CD8 T cells

##make cisTopic####
Kidney.CD8$barcode <- colnames(Kidney.CD8)
Kidney.CD8.ATAC.counts <- Kidney.CD8[["MACS2ATAC"]]@counts
rownames(Kidney.CD8.ATAC.counts) <- str_replace(rownames(Kidney.CD8.ATAC.counts), "-", ":")

barcode <- colnames(Kidney.CD8)
barcode <- as.data.frame(barcode)
barcode$is__cell_barcode <- "1"
barcode$total <- Kidney.CD8$nCount_ATAC
barcode$TSS_fragments <- Kidney.CD8$atac_TSS_fragments

cisTopicObject.CD8 <- createcisTopicObject(Kidney.CD8.ATAC.counts, barcode, project.name='Kidney.CD8', 
                                           keepCountsMatrix = FALSE,
                                           min.cells = 1, min.regions = 1, is.acc = 1)
cisTopicObject.CD8 <- runCGSModels(cisTopicObject.CD8, topic=c(5, 10, 20, 30, 40), seed=987, 
                                   nCores=20, burnin = 120, iterations = 150, addModels=FALSE)
cisTopicObject.CD8 <- selectModel(cisTopicObject.CD8)
cisTopicObject.CD8 <- selectModel(cisTopicObject.CD8, select= 40)

cisTopicObject.CD8 <- runUmap(cisTopicObject.CD8, target='cell')
topic.mat.CD8 <- modelMatSelection(cisTopicObject.CD8, 'cell', 'Probability')
topic.mat.CD8 <- t(topic.mat.CD8)
topic.mat.CD8 <- as.matrix(topic.mat.CD8)

set.seed(123)
cellkNN.CD8 <- get.knn(topic.mat.CD8,k = 30)$nn.index
dim(cellkNN.CD8)
rownames(cellkNN.CD8) <- barcode$barcode

##make required ATAC and RNA formats####
CD8.ATAC <- Kidney.CD8[["MACS2ATAC"]]@counts
CD8.ATAC <- as.data.frame(as.matrix(CD8.ATAC))
CD8.ATAC$seqnames <- sapply(strsplit(rownames(CD8.ATAC),"-"), `[`, 1)
CD8.ATAC$start <- sapply(strsplit(rownames(CD8.ATAC),"-"), `[`, 2)
CD8.ATAC$end <- sapply(strsplit(rownames(CD8.ATAC),"-"), `[`, 3)
CD8.ATAC <- subset(CD8.ATAC, grepl('chr', rownames(CD8.ATAC)))

CD8.ATAC.se <- makeSummarizedExperimentFromDataFrame(CD8.ATAC)
counts(CD8.ATAC.se) <- assay(CD8.ATAC.se)
assay(CD8.ATAC.se) <- as(assay(CD8.ATAC.se), 'sparseMatrix')

CD8.RNA <- GetAssayData(Kidney.CD8, assay = "RNA.clean", slot = "data")
# Remove genes with zero expression across all cells
CD8.RNA <- CD8.RNA[Matrix::rowSums(CD8.RNA)!=0,]
CD8.RNA <- as.matrix(CD8.RNA)

##compute peak to gene associations and determine DORCs####
cisCorr.CD8 <- FigR::runGenePeakcorr(ATAC.se = CD8.ATAC.se,
                                     RNAmat = CD8.RNA,
                                     genome = "mm10", 
                                     windowPadSize = 25000,
                                     nCores = 4,
                                     p.cut = NULL, # Set this to NULL and we can filter later
                                     n_bg = 50)

cisCorr.CD8.filt <- cisCorr.CD8 %>% filter(pvalZ <= 0.05)
dorcGenes.CD8 <- dorcJPlot(dorcTab = cisCorr.CD8.filt,
                           cutoff = 3, # No. sig peaks needed to be called a DORC
                           labelTop = 20,
                           returnGeneList = TRUE, # Set this to FALSE for just the plot
                           force=2)
dorcMat.CD8 <- getDORCScores(ATAC.se = CD8.ATAC.se, # Has to be same SE as used in previous step
                             dorcTab = cisCorr.CD8.filt,
                             geneList = dorcGenes.CD8,
                             nCores = 4)

dorcMat.CD8.s <- smoothScoresNN(NNmat = cellkNN.CD8[,1:30], mat = dorcMat.CD8, nCores = 4)
RNAmat.CD8.s <- smoothScoresNN(NNmat = cellkNN.CD8[,1:30], mat = CD8.RNA, nCores = 10)
umap.CD8.d <- as.data.frame(Kidney.CD8@reductions[["umap.wnn"]]@cell.embeddings)

##compute TF-gene associations####
motif.names <- Kidney.CD8@assays$MACS2ATAC@motifs@motif.names
motif.names <- unlist(motif.names)
tf_names <- unique(motif.names[!grepl("::", motif.names)])
#change lower to upper case
norm <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))
tf_key  <- norm(tf_names)
rna_key <- norm(rownames(RNAmat.CD8.s))
#keep TF genes present in RNAmat.s
keep_tf <- rna_key %in% tf_key
RNAmat_TF.CD8 <- RNAmat.CD8.s[keep_tf, , drop = FALSE] 
stopifnot(identical(colnames(RNAmat_TF.CD8), colnames(dorcMat.CD8.s)))
#filter cisCorr to dorc genes
cisCorr.CD8.use <- subset(cisCorr.CD8.filt, Gene %in% rownames(dorcMat.CD8.s))
#reduce ATAC.se object, remove anything thats not counts
if (!"counts" %in% assayNames(CD8.ATAC.se)) {
  # if the first assay is counts but unnamed, rename it
  assayNames(CD8.ATAC.se)[1] <- "counts"
}

cnt <- assay(CD8.ATAC.se, "counts")
if (!inherits(cnt, "dgCMatrix")) cnt <- Matrix(as.matrix(cnt), sparse = TRUE)
assays(CD8.ATAC.se) <- S4Vectors::SimpleList(counts = cnt)

figR.d.CD8 <- runFigRGRN(ATAC.se = CD8.ATAC.se, # Must be the same input as used in runGenePeakcorr()
                         dorcTab = cisCorr.CD8.use, # Filtered peak-gene associations
                         genome = "mm10",
                         dorcMat = dorcMat.CD8.s,
                         rnaMat = RNAmat_TF.CD8, 
                         nCores = 20)




