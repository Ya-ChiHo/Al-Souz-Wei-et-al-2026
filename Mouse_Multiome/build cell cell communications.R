library(Seurat)
options(Seurat.object.assay.version = 'v5')
library(Signac)
library(ggplot2)
library(ggraph)
library(cowplot)
library(patchwork)
library(dplyr)
library(tidyverse)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(TFBSTools)
library(scSHC)
library(stringr)
library(scriabin)
library(pbapply)
library(scriabin)
library(ComplexHeatmap)
library(cowplot)
library(flashClust)
library(ade4)
library(glmGamPoi)
library(CelliD)

#these scripts apply scriabin package to identify ligand-receptor sender and receiver cells using RNA assay.
#require processed objects from "per tissue data preprocessing.R"

##modified functions from Scriabin####
IPFeaturePlot1 <- function (seu, ip, cols = c("grey90", "blue", "orangered3"), 
                            order = T) 
{
  seu$lig_feature <- seu[["IPligands"]]@data[ip, ]
  seu$rec_feature <- seu[["IPreceptors"]]@data[ip, ]
  p <- FeaturePlot(seu, features = c("lig_feature", "rec_feature"), 
                   blend = T, combine = F, cols = cols, order = order, reduction = 'umap.wnn')
  p[[3]] + NoLegend() + labs(x = NULL, y = NULL, title = NULL) + 
    theme(aspect.ratio = 1, axis.ticks = element_blank(), 
          axis.text = element_blank())
}
GenerateCellSignature1 <- function (seu, variant_genes, dq = 0.05) 
{
  seu <- RunMCA(seu, features = variant_genes, slot = "RNA.clean")
  ds2 <- t(CelliD:::GetCellGeneDistance(seu, reduction = "mca", 
                                        dims = 1:30))
  rq <- matrixStats::rowQuantiles(ds2, probs = dq)
  dsp <- ds2
  dsp[rq[row(dsp)] < dsp] <- 0
  dsp[dsp > 0] <- 1
  return(dsp)
}
TopLigandsByIdent1 <- function (seu, active_ligands = NULL, sender = NULL, receiver = NULL, 
                                group.by = "orig.ident", pearson.threshold = 0.1, ligands.display = 25, 
                                assay = "RNA.clean") 
  
{
  receiver_cells <- colnames(seu)[seu@meta.data[, group.by] == 
                                    receiver]
  ral_sub <- active_ligands[, receiver_cells]
  ral_means <- rowMeans(ral_sub)
  ral_pct <- ral_sub
  ral_pct[ral_pct > pearson.threshold] <- 1
  ral_pct[ral_pct < 1] <- 0
  ral_pct <- 100 * rowSums(ral_pct)/ncol(ral_sub)
  ral_res <- data.frame(ligand = names(ral_means), mean = ral_means, 
                        pct = ral_pct) %>% as_tibble() %>% top_n(wt = mean, n = ligands.display)
  avg_exprs <- AverageExpression(seu, features = unique(ral_res$ligand), 
                                 group.by = group.by, assay = assay)[[1]][, sender]
  avg_exprs <- data.frame(ligand = names(avg_exprs), exprs = avg_exprs)
  ral_res <- merge(ral_res, avg_exprs, by = "ligand") %>% arrange(-mean)
  ral_res$ligand <- factor(ral_res$ligand, levels = ral_res$ligand)
  ggplot(ral_res, aes(x = ligand, y = mean, size = pct, color = exprs)) + 
    geom_point() + scale_radius() + theme_cowplot() + scale_color_gradientn(colours = c("grey", 
                                                                                        "yellow", "red3")) + ggpubr::rotate_x_text() + labs(x = "Ligand", 
                                                                                                                                            y = paste0("Mean ligand activity in\n", receiver, " receiver cells"), 
                                                                                                                                            color = paste0("Expression by\n", sender, "\nsender cells"), 
                                                                                                                                            size = paste0("Percent receiver cells\nwith activity > ", 
                                                                                                                                                          pearson.threshold))
}


##generate scrainbin interaction programs####
#add back subset annotations from Kidney.CD8 and Kidney.CD4 to Kidney
Kidney.CD4$celltype <- "CD4"
Kidney.CD4$toCombine <- str_c(Kidney.CD4$celltype, '_', Kidney.CD4$major_celltype)
Kidney.CD8$celltype <- "CD8"
Kidney.CD8$toCombine <- str_c(Kidney.CD8$celltype, '_', Kidney.CD8$major_celltype)

CD8.cell_annotation <- as.data.frame(Kidney.CD8$toCombine)
CD8.cell_annotation$celltype <- CD8.cell_annotation$`Kidney.CD8$toCombine`
CD8.cell_annotation$`Kidney.CD8$toCombine` <- NULL
CD4.cell_annotation <- as.data.frame(Kidney.CD4$toCombine)
CD4.cell_annotation$celltype <- CD4.cell_annotation$`Kidney.CD4$toCombine`
CD4.cell_annotation$`Kidney.CD4$toCombine` <- NULL

T_cell_annotation <- rbind(CD4.cell_annotation, CD8.cell_annotation)
T_cell_annotation$BC <- rownames(T_cell_annotation)

Kidney$major_celltype <- ifelse(colnames(Kidney) %in% rownames(T_cell_annotation), 
                                    as.character(T_cell_annotation$celltype)[
                                      match(colnames(Kidney), rownames(T_cell_annotation))],FALSE)


#perform ALRA denoising for sparse data
Kidney <- SeuratWrappers::RunALRA(Kidney)

#compute interaction programs
Kidney.ip <- FindAllInteractionPrograms(Kidney, iterate.threshold = 300, 
                                            species = "mouse", database = "OmniPath",
                                            group.by = "major_celltype", assay = "RNA.clean", 
                                            sim_threshold = 0.4)

#identify significant interaction programs
Kidney.ip_sig <- InteractionProgramSignificance(Kidney.ip, n.replicate = 100)

#keep IP that are significant in atleast 1 condition type
Kidney.ip_pvals <- Kidney.ip_sig %>% as_tibble() %>%
  dplyr::select(name,ends_with("pval")) %>% unique() %>%
  pivot_longer(!name, names_to = "celltype", values_to = "pval") %>%
  group_by(name) %>% dplyr::mutate(min_p = min(pval)) %>%
  dplyr::select(name,min_p) %>% unique() %>% 
  dplyr::filter(min_p < 0.05) %>% pull(name)
Kidney.ip_sig <- Kidney.ip_sig %>% dplyr::filter(name %in% Kidney.ip_pvals)

#score cells by expression of interaction program
Kidney <- ScoreInteractionPrograms(Kidney, mods = Kidney.ip_sig)




##ranking ligand activity####
#convert nichnet object human gene name to mosue gene name
mouse_human_genes <- read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
mouse <- split.data.frame(mouse_human_genes,mouse_human_genes$Common.Organism.Name)[[2]]
human <- split.data.frame(mouse_human_genes,mouse_human_genes$Common.Organism.Name)[[1]]
mouse <- mouse[,c(1,4)]
human <- human[,c(1,4)]
mh_data <- merge.data.frame(mouse,human,by = "DB.Class.Key",all.y = TRUE) 
colnames(ligand_target_matrix) <- 
  ifelse(colnames(ligand_target_matrix) %in% mh_data$Symbol.y, 
         mh_data$Symbol.x[match(colnames(ligand_target_matrix), mh_data$Symbol.y)],FALSE)
rownames(ligand_target_matrix) <- 
  ifelse(rownames(ligand_target_matrix) %in% mh_data$Symbol.y,  
         mh_data$Symbol.x[match(rownames(ligand_target_matrix), mh_data$Symbol.y)],FALSE)
ligand_target_matrix <- as.matrix(ligand_target_matrix)
rownames(ligand_target_matrix) <- gsub("^(\\w)(\\w*)", "\\U\\1\\L\\2", rownames(ligand_target_matrix), perl = TRUE)
colnames(ligand_target_matrix) <- gsub("^(\\w)(\\w*)", "\\U\\1\\L\\2", colnames(ligand_target_matrix), perl = TRUE)

scriabin::load_nichenet_database()
data("lr_network")  
head(lr_network$from)
unique(lr_network$from)[1:20] 

Idents(Kidney) <- "major_celltype"
DefaultAssay(Kidney) <- "RNA.clean"
variant_genes <- IDVariantGenes(Kidney, group.by = "major_celltype", assay = "RNA.clean", n.gene = 5000)
gene_signature <- GenerateCellSignature1(Kidney, variant_genes = variant_genes)
active_ligands <- RankActiveLigands(Kidney, signature_matrix = gene_signature, 
                                    species = "mouse", assay = "RNA.clean")
TopLigandsByIdent1(Kidney, active_ligands = active_ligands, ligands.display = 10, 
                   sender = "CD4-Il10 CD4", receiver = "CD8-Tscm", group.by = "major_celltype")
