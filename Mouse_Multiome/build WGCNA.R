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
library(WGCNA)
library(flashClust)
library(Hmisc)
library(dplyr)
library(openxlsx)
library(qlcMatrix)
library(reshape2)
library(scales)
library(Matrix)
library(tibble)

##require objects built from "merged tissue data processing.R"
##this script modifies WGCNA functions from Kazer et al. 2019, Nature Medicine 
##shown are WGCNA built for CD8 T cells in kidney, spleen, and rLN

##required WGCNA functions from Kazer et al. 2019, Nature Medicine####
# Choosing the appropriate power for generating the adjacency matrix
FindPower <- function(datExpr){
  #choose soft-threshold power
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  sft=pickSoftThreshold(datExpr,dataIsExpr = TRUE,powerVector = powers,corFnc = cor,
                        corOptions = list(use = 'p', method = "pearson"),networkType = "signed")
  
  # Plot the results
  par(mfrow = c(1,2));
  cex1 = 0.9;
  
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
  
  # Red line corresponds to using an R^2 cut-off
  abline(h=0.80,col="red")
  
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
}
# Generating the adjacency matrix and performing clustering
ClusterTOM <- function(datExpr, softPower){
  #dev.off()
  #Calclute the adjacency matrix
  adj= adjacency(datExpr,type = "signed", power = softPower);
  
  #Turn adjacency matrix into a topological overlap matrix (TOM) to minimize the effects of noise and spurious associations.
  TOM=TOMsimilarityFromExpr(datExpr,networkType = "signed", TOMType = "signed", power = softPower, corType="bicor");
  
  colnames(TOM) = rownames(TOM) = colnames(datExpr)
  dissTOM=1-TOM
  
  #Hierarchical clustering of the genes based on the TOM dissimilarity measure
  geneTree = flashClust(as.dist(dissTOM),method="complete");
  
  #Plot the resulting clustering tree (dendrogram)
  plot(geneTree, xlab="", sub="",cex=0.3);
  
  return(list(dissTOM = dissTOM, geneTree = geneTree)) #returns list with dissimilarity TOM, and the clustered gene tree.
}
# Cut the resulting clustering dendrogram using the "tree" method for cutreeDynamic. Minimum module size can be specified.
CutTOMTree <- function(datExpr, dissTOM, geneTree, minModuleSize = 10){
  #dev.off()
  # Module identification using dynamic tree cut, you can also choose the hybrid method
  dynamicMods = cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = minModuleSize);
  #dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, method="hybrid", deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize);
  
  #Get the module labels and the size of each module. Lable 0 is reserved for unassigned genes
  print(table(dynamicMods))
  
  #Plot the module assignment under the dendrogram; note: The grey color is reserved for unassigned genes
  dynamicColors = labels2colors(dynamicMods)
  table(dynamicColors)
  plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
  
  #Set the diagonal of the dissimilarity to NA 
  diag(dissTOM) = NA;
  
  #extract modules
  module_colors= setdiff(unique(dynamicColors), "grey")
  modules = lapply(module_colors, function(x){colnames(datExpr)[which(dynamicColors==x)]})
  names(modules) = module_colors
  return(list(dyanmicColors = dynamicColors, modules = modules)) #returns list with module colors, and the modules themselves
}
# Merge modules with low dissimilarity. Cutoff for dissimilarity merge can be specified
MergeSimilarModules <- function(datExpr, dynamicColors, geneTree, MEDissThres = 0.5){
  #cacluate eigengenes
  MEList = moduleEigengenes(datExpr, colors=dynamicColors)
  MEs = MEList$eigengenes
  
  # Calculate dissimilarity of module eigengenes
  MEDiss = 1-cor(MEs);
  
  # Cluster module eigengenes
  METree = hclust(as.dist(MEDiss), method = "average");
  
  # Plot the result
  #sizeGrWindow(7, 6)
  plot(METree, main = "Clustering of module eigengenes",
       xlab = "", sub = "")
  abline(h = MEDissThres, lwd=2, col="red")
  
  # Call an automatic merging function
  merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
  # The merged module colors
  mergedColors = merge$colors;
  # Eigengenes of the new merged modules:
  mergedMEs = merge$newMEs;
  #plot showing how merged modules exist
  plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                      c("Dynamic Tree Cut", "Merged dynamic"),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  
  #extract merged modules
  merged_module_colors= setdiff(unique(mergedColors), "grey")
  merged_modules = lapply(merged_module_colors, function(x){colnames(datExpr)[which(mergedColors==x)]})
  names(merged_modules) = merged_module_colors
  
  return(list(mergedColors = mergedColors, merged_modules = merged_modules)) #returns list with merged colors, and the merged modules themselves
}
# Test to determine if the genes within the module are truly the least dissimilar compared to randomly generated modules of the same size.
TestModuleSignificance <- function(mod, dissTOM, expr.data, n_perm = 10000, pval = 0.05, n.bin = 10){
  #vectorize the actual distribution of (dis)similarities, and remove zeros!
  true.diss = as.vector(dissTOM[mod,mod])
  true.diss = true.diss[-which(true.diss == 0)]
  
  #size of module for permutations
  mod.size = length(mod)
  
  #bin all genes by expression
  expr.avg = rowMeans(expr.data)
  expr.avg = expr.avg[order(expr.avg)]
  expr.avg.cut = as.numeric(x = cut2(x = expr.avg, m=round(length(expr.avg)/n.bin)))
  names(expr.avg.cut) = names(expr.avg)
  
  #create a table of binnings of all genes and of our module genes
  all.bin.table = table(expr.avg.cut)
  mod.bin.table = table(expr.avg.cut[mod])
  
  #randomly generate module with same expression binning structure and run t.test and record results
  test.results = data.frame(statistic = rep(NA, n_perm), pval = rep(NA, n_perm)) #container for results
  
  for (i in 1:n_perm){ #by permutation
    random.mod = list() #create an empty list we will fill with gene names (each element will be gene names by bin)
    
    #for each bin of the mod bin table, randomly select that number of genes from the full set with that expression
    for (j in names(mod.bin.table)){ 
      bin.genes = sample(names(expr.avg.cut)[which(expr.avg.cut == as.numeric(j))], mod.bin.table[j], replace = FALSE)
      random.mod[[as.numeric(j)]] = bin.genes #stick those genes into the random.mod list
    }
    #unlist and vectorize the distribution of (dis)similarities (remember to remove zeros)
    random.mod = unlist(random.mod)
    random.diss = as.vector(dissTOM[random.mod,random.mod])
    random.diss = random.diss[-which(random.diss == 0)]
    
    #perform a one-sided wilcox.test and record the statistic and p-value.
    #Note, IMPORTANT: here we perform the test asking if the true diss is LESS THAN the random diss, as we are trying to minimize dissimilarity
    test = wilcox.test(x = true.diss, y = random.diss, alternative = "less")
    test.results[i,] = c(test$statistic, test$p.value)
  }
  
  #correct for multiple hypothesis testing, and then report the proportion of bad tests
  test.results$FDR = p.adjust(test.results$pval, method = "fdr")
  num.failed.tests = sum(test.results$FDR > pval)
  print(paste(paste(num.failed.tests, n_perm, sep="/"), "permutations failed the Mann-Whitney test.", sep=" "))
  
  #is the percentage of failed tests less than or equal to the p-val specified?
  return(num.failed.tests/n_perm <= pval) #returns a vector of booleans indicating if each module was significant based on the specific p.val
}

##wrapper and plotting functions, from Wei et al. 2023, Immunity####
#function just wraps all the steps in one easy function w/ prompts
modid<-function(seuratobj, cluster, prefix="modules"){
  #initial ensuring we're scaled and PCA ready
  test_clus<-subset(seuratobj, idents=cluster)
  test_clus<-FindVariableFeatures(test_clus)
  test_clus<-ScaleData(test_clus)
  test_clus<-RunPCA(test_clus, reduction.key = "PCWGCNA_")
  print(ElbowPlot(test_clus))
  nPCS<-readline("how many NPCs?")
  nPCS<-as.integer(nPCS)
  #expanding the dim predicted list
  test_clus<-ProjectDim(test_clus)
  
  genes<-c()
  for (i in 1:nPCS){genelist<-TopFeatures(test_clus, dim = i, nfeatures = 50, balanced = T, projected = T )
  for (i in 1:length(genelist)){genes<-c(genes, genelist[[i]])}}
  genes<-unique(genes)
  test_clus<-as.matrix(Kidney.CD8[["RNA"]]$data[genes,])
  FindPower(datExpr=t(test_clus))
  softpower<-readline("what softpower?")
  softpower<-as.integer(softpower)
  #run all of sams functions
  test_clus_tom<-ClusterTOM(datExpr = t(test_clus), softPower = softpower)
  test_clus_mods<-CutTOMTree(datExpr = t(test_clus), geneTree = test_clus_tom$geneTree, dissTOM = test_clus_tom$dissTOM, minModuleSize = 8)
  test_clus_merge_mods<-MergeSimilarModules(datExpr=t(test_clus), dynamicColors = test_clus_mods$dyanmicColors, geneTree = test_clus_tom$geneTree, MEDissThres = 0.5)
  print(test_clus_merge_mods$merged_modules)                                          
  test_clus_merge_mods.isSig = sapply(test_clus_merge_mods$merged_modules, function(module){
    TestModuleSignificance(mod = module, dissTOM = test_clus_tom$dissTOM, expr.data = test_clus,
                           n_perm = 10000, pval = 0.05, n.bin = 10)
  })
  test_clus_merge_mods.isSig = test_clus_merge_mods$merged_modules[test_clus_merge_mods.isSig]
  print(test_clus_merge_mods.isSig)
  names(test_clus_merge_mods.isSig)<-paste(prefix,names(test_clus_merge_mods.isSig ), sep="_")
  return(test_clus_merge_mods.isSig)
}
#changes to function:
MergeSimilarModules <- function(datExpr, dynamicColors, geneTree, MEDissThres = 0.5){
  #cacluate eigengenes
  MEList = moduleEigengenes(datExpr, colors=dynamicColors)
  MEs = MEList$eigengenes
  
  # Calculate dissimilarity of module eigengenes
  MEDiss = 1-cor(MEs);
  
  # Cluster module eigengenes
  METree = hclust(as.dist(MEDiss), method = "average");
  
  # Plot the result
  #sizeGrWindow(7, 6)
  plot(METree, main = "Clustering of module eigengenes",
       xlab = "", sub = "")
  abline(h = MEDissThres, lwd=2, col="red")
  
  # Call an automatic merging function
  merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
  # The merged module colors
  mergedColors = merge$colors;
  # Eigengenes of the new merged modules:
  mergedMEs = merge$newMEs;
  #plot showing how merged modules exist
  plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                      c("Dynamic Tree Cut", "Merged dynamic"),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  
  #extract merged modules
  merged_module_colors= setdiff(unique(mergedColors), "grey")
  merged_modules = lapply(merged_module_colors, function(x){colnames(datExpr)[which(mergedColors==x)]})
  names(merged_modules) = merged_module_colors
  
  return(list(mergedColors = mergedColors, merged_modules = merged_modules)) #returns list with merged colors, and the merged modules themselves
}
extract.cor<-function(cormatrix, gene, ngene=50){
  test<-list()
  for (i in 1:length(cormatrix)){test[[i]]<-names(sort(cormatrix[[i]][gene,],decreasing = T)[1:ngene])}
  return(test)
}
corheatmap<-function(cor, genes, prefix, n_ext=50,exp=TRUE){
  samples<-names(cor)
  plot<-c()
  if(exp==TRUE){for(i in 1:length(genes)){
    plot<-unique(c(plot,extract.cor(cor,genes[[i]], ngene = n_ext)))
  }
    plot2<-c()
    for(i in 1:length(plot)){
      plot2<-c(plot2, plot[[i]])
    }
    plot2<-unique(plot2)}else{plot2<-genes}
  
  for (z in 1:length(cor)){
    corx<-cor[[z]][plot2,plot2]
    x<-setNames(melt(corx), c("x","y","weight"))
    x$x<-factor(x$x, levels = plot2)
    x$y<-factor(x$y, levels = plot2)
    x[is.na(x)]<-0
    p<-ggplot(x, aes(x=x,y, fill=weight))+geom_tile()+
      scale_fill_gradient2(low = "blue",high = "red", limits=c(-0.3, 0.3), oob=squish)+
      theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.4))
    savename<-paste(prefix,"_",samples[z],".pdf", sep = "")
    ggsave(filename = savename, plot = p, device = "pdf", width = 10, height = 9)
  }
}
splitcorheatmap<-function(cor, genes, prefix, n_ext=1,exp=TRUE){
  samples<-names(cor)
  plot<-c()
  if(exp==TRUE){for(i in 1:length(genes)){
    plot<-unique(c(plot,extract.cor(cor,genes[[i]], ngene = n_ext)))
  }
    plot2<-c()
    for(i in 1:length(plot)){
      plot2<-c(plot2, plot[[i]])
    }
    plot2<-unique(plot2)}else{plot2<-genes}
  
  corx<-cor[[1]][plot2,plot2]
  cory<-cor[[2]][plot2,plot2]
  corx[lower.tri(corx)]<-cory[lower.tri(cory)]
  x<-setNames(melt(corx), c("x","y","weight"))
  x[is.na(x)]<-0
  x$x<-factor(x$x, levels = plot2)
  x$y<-factor(x$y, levels = plot2)
  p<-ggplot(x, aes(x=x,y, fill=weight))+geom_tile()+
    scale_fill_gradient2(low = "blue",high = "red", limits=c(-0.3, 0.3), oob=squish)+
    theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.4))+theme(axis.ticks = element_blank())
  savename<-paste(prefix,"_",samples[1],"_",samples[2],"split.pdf", sep = "")
  ggsave(filename = savename, plot = p, device = "pdf", width = 10, height = 9)
}

##construct correlation data####
DefaultAssay(Kidney.CD8) <- "RNA.clean"
corKidney.CD8 <- corSparse(t(Kidney.CD8@assays$RNA.clean$data))
rownames(corKidney.CD8) <- rownames(Kidney.CD8)
colnames(corKidney.CD8) <- rownames(Kidney.CD8)

DefaultAssay(LN.CD8) <- "RNA.clean"
corLN.CD8 <- corSparse(t(LN.CD8@assays$RNA.clean$data))
rownames(corLN.CD8) <- rownames(LN.CD8)
colnames(corLN.CD8) <- rownames(LN.CD8)

DefaultAssay(spleen.CD8) <- "RNA.clean"
corSpleen.CD8 <- corSparse(t(spleen.CD8@assays$RNA.clean$data))
rownames(corSpleen.CD8) <- rownames(spleen.CD8)
colnames(corSpleen.CD8) <- rownames(spleen.CD8)

##identify WGCNA modules in Teff and Tscm cells from kidney####
Idents(Kidney.CD8) <- "major_celltype"; DefaultAssay(Kidney.CD8) <- "RNA.clean"
modid(Kidney.CD8, "Terminal teff-exhausted", prefix="WGCNA.teff.Kidney.CD8")
#picked pc = 8, power = 2

Idents(Kidney.CD8) <- "major_celltype"; DefaultAssay(Kidney.CD8) <- "RNA.clean"
modid(Kidney.CD8, "Tscm", prefix="WGCNA.Tscm.Kidney.CD8")
#pc = 8, power = 1



