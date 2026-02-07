R Scripts related to the Al Souz et al 2026 manuscript.

Please see STAR Methods in Al Souz et al 2026 manuscript for software and package versions.

These scripts are minimally commented to be functional and reproducible, they are not tutorials. Please see the Methods section in the manuscript for additional details and non-R related analyses.

Scripts are located in 3 separate directories, depending on the type of single-cell data (mouse CITE-seq+TCR, mouse 10x Multiome, and human RNA-seq), as described in the manuscript.

'Mouse CITE-seq and TCR' and 'Human RNA-seq' directories contain functional scripts to generate the corresponding figures in the manuscript. 


In 'Mouse_Multiome":

The two 'preprocessing.R' files describe example steps taken to process Cellranger data to generate Seurat objects: Building and merging ATAC, RNA, and CITE Seurat objects, recall peaks by MACS2, doublet removal, low quality cell filter, 
data normalization and batch effect removal, building gene accessibility data, adding chromVAR transcription factor motif scores, and celltype annotation.

build GSEA.R - minimal example scripts to build required R objects for Gene Set Enrichment Analysis (by msigdbr and fgsea).

build GRN.R - minimal example scripts to build required R objects for Gene Regulatory Network analyses (by cisTopic and FigR).

build WGCNA.R - minimal example scripts to build required R objects for WGCNA analyses by methods from Kazer et al. (Nature Medicine, 2019).

build cell cell communication.R - minimal example scripts to build required R objects for ligand-receptor interaction analyses (by Scriabin and NicheNet).

build monocle2 trajectory.R - minimal example scripts to build required R objects for pseudotime trajectory analysis (by Monocle3).

In 'Mouse_Multiome', scripts in the 'Figure #.R' require processed Seurat objects generated from the above (as described in script comment). 


All 'Figure #.R' generate the plots as shown in manuscript main figures and, when relevant, the functions/data processing steps required to generate these plots. Scripts in 'Figures.R' may be modified to reproduce all supplemental figures. 

Raw reads and processed data (e.g., cellranger filtered feature bc matrices, hashtag keys, fragment.tsv etc.) are hosted on GEO:GSE318514. 

For all other information, see supplemental files and the STAR Methods in the associated manuscript.
If there are particular aspects of the analyses you would like to see that are not here, or have any other questions, please email Yulong.Wei@yale.edu or jafar.alsouz@yale.edu.

Comment on Seurat v5 and Signac for Seurat v5 (2026-02-07): Note that some functions have been removed/replaced, while default parameters changed in others (e.g., FindMarkers), and data structure different (e.g., layers), in Seurat v5 vs Seurat v4. 
Version difference may lead to some scripts breaking and may lead to slight graphical differences (e.g., volcanoplot because of change in FC in FindMarkers), but do not otherwise affect the main results (e.g., key differentially expressed genes).
