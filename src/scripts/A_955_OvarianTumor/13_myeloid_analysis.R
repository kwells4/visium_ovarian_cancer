library(ggplot2)
library(Seurat)
library(tidyr)
library(cowplot)
library(dplyr)
library(openxlsx)
library(here)
library(clustree)
library(AUCell)


# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

sample <- "A_955_OvarianTumor"
cell_types <- "SCT_celltype_combined"
celltype2 <- "Spatial_celltype_combined"
clusters <- "SCT_cluster"

seurat_assay <- "SCT"
vars.to.regress <- NULL
ffpe <- TRUE # Setting to TRUE will assume a probe set was used.
HTO <- FALSE
ADT <- FALSE


source(here("src/scripts/functions.R"))

base_dir_proj <- here(file.path("results", sample))

no_multi_dir <- file.path(base_dir_proj, "R_analysis_no_multi")

# Read in the data
seurat_data_no_multi <- readRDS(file.path(no_multi_dir, "rda_obj",
                                          "seurat_no_mult_processed.rds"))

cells <- seurat_data_no_multi$SCT_celltype_no_mult == "Monocyte"
counts <- GetAssayData(seurat_data_no_multi, slot = "counts")[ ,cells]
genes_use <- rowSums(counts) > 3

genes <- names(which(genes_use))


myeloid_seurat <- subset(seurat_data_no_multi,
                         cells = colnames(counts),
                         features = genes)

# Renormalize ------------------------------------------------------------------
myeloid_seurat <- SCTransform(myeloid_seurat, assay = "Spatial", verbose = FALSE,
                          vars.to.regress = vars.to.regress)

DefaultAssay(myeloid_seurat) <- "Spatial"

myeloid_seurat <- FindVariableFeatures(myeloid_seurat) %>%
  ScaleData(vars.to.regress = vars.to.regress)

# PCA --------------------------------------------------------------------------

# PCA of gene expression
myeloid_seurat <- PCA_dimRed(myeloid_seurat, assay = seurat_assay)

RNA_plots <- plot_PCA(HTO = FALSE, assay = seurat_assay,
                      sample_object = myeloid_seurat,
                      data_type = "spatial", ffpe = ffpe)

pdf(file.path(no_multi_dir, "images", "RNA_pca_myeloid.pdf"))
print(RNA_plots)
dev.off()

# UMAP -------------------------------------------------------------------------
RNA_pcs <- 8

set.seed(0)
# UMAP of gene expression
umap_data <- group_cells(myeloid_seurat, sample, nPCs = RNA_pcs,
                         resolution = 0.6, assay = seurat_assay, HTO = FALSE)

myeloid_seurat <- umap_data$object

# Test a range of resolutions
myeloid_seurat <- FindClusters(myeloid_seurat, resolution = c(0.1, 0.3, 0.5, 0.8,1))
clustree(myeloid_seurat)

# UMAP of gene expression with final resolution selelction
set.seed(0)
umap_data <- group_cells(myeloid_seurat, sample, no_multi_dir, nPCs = RNA_pcs,
                         resolution = 0.6, assay = seurat_assay, HTO = HTO)

myeloid_seurat <- umap_data$object

gene_plots <- umap_data$plots

plotDimRed(myeloid_seurat, col_by = "SCT_cluster", plot_type = "rnasct.umap",
           size = 1)


SpatialDimPlot(myeloid_seurat, label = TRUE, label.size = 3, cols = "Set1")

# M1/M2 analysis ---------------------------------------------------------------

# This list is taken from https://insight.jci.org/articles/view/126556/figure/3
# Single cell RNA sequencing identifies unique inflammatory airspace 
# macrophage subsets

gene_list <- list(M1_lung = c("AZIN1",
                              "CD38",
                              "CD86",
                              "CXCL10",
                              "CXCL9",
                              "FPR2",
                              "GPR18",
                              "IL12B",
                              "IL18",
                              "IL1B",
                              "IRF5",
                              "NFKBIZ",
                              "NOS2",
                              "PTGS2",
                              "SOCS3",
                              "TLR4",
                              "TNF"),
                  M2_lung = c("ALOX15",
                              "ARG1",
                              "CHIL3",
                              "CLEC7A",
                              "EGR2",
                              "IL10",
                              "IRF4",
                              "KLF4",
                              "MRC1",
                              "MYC",
                              "SOCS2",
                              "TGM2"))

# From M1 Macrophage and M1/M2 ratio defined by transcriptomic signatures 
# resemble only part of their conventional clinical characteristics in breast 
# cancer" and is copied from supplemental table 1 and supplemental table 2. 
# https://doi.org/10.1038/s41598-020-73624-w

list_path <- "/Users/wellskr/Documents/Analysis/Rachel_Friedman/Friedman_mertk/Rachel_F_nPOD_10X/files/M1_M2_human_signatures.xlsx"

list_names <- openxlsx::getSheetNames(list_path)

gene_list2 <- lapply(list_names, function(x){
  gene_file <- openxlsx::read.xlsx(list_path, sheet = x)
  gene_list <- gene_file$Gene
  return(gene_list)
})

names(gene_list2) <- paste(list_names, "breast_cancer", sep = "_")

all_lists <- c(gene_list, gene_list2)

all_lists <- lapply(all_lists, function(x){
  genes <- x[x %in% rownames(myeloid_seurat)]
  return(genes)
})

DefaultAssay(myeloid_seurat) <- "SCT"
exp_mat <- GetAssayData(myeloid_seurat, slot = "data")

cells_rankings <- AUCell_buildRankings(exp_mat,
                                       plotStats = TRUE)

cells_AUC <- AUCell_calcAUC(all_lists, cells_rankings)

add_metadata <- AUCell::getAUC(cells_AUC) %>%
  t() %>%
  data.frame

colnames(add_metadata) <- paste0("AUCell_", colnames(add_metadata))

add_metadata$AUCell_M1_M2_lung <- add_metadata$AUCell_M1_lung -
                                  add_metadata$AUCell_M2_lung

add_metadata$AUCell_M1_M2_breast_cancer <- add_metadata$AUCell_Human_M1_breast_cancer - 
                                           add_metadata$AUCell_Human_M2_breast_cancer

cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE)

#myeloid_seurat <- AddMetaData(myeloid_seurat, add_metadata)

# plot M1 vs M2
ggplot(add_metadata, aes(x = AUCell_M1_lung, y = AUCell_M2_lung,
                         color = AUCell_M1_M2_lung)) +
  geom_point() +
  scale_color_viridis(option = "magma")


ggplot(add_metadata, aes(x = AUCell_Human_M1_breast_cancer,
                         y = AUCell_Human_M2_breast_cancer,
                         color = AUCell_M1_M2_breast_cancer)) +
  geom_point() +
  scale_color_viridis(option = "magma")

ggplot(add_metadata, aes(x = AUCell_M1_lung, y = AUCell_M2_lung,
                         color = AUCell_M1_M2_breast_cancer)) +
  geom_point() +
  scale_color_viridis(option = "magma")


ggplot(add_metadata, aes(x = AUCell_Human_M1_breast_cancer,
                         y = AUCell_Human_M2_breast_cancer,
                         color = AUCell_M1_M2_lung)) +
  geom_point() +
  scale_color_viridis(option = "magma")

write.csv(add_metadata, file = file.path(no_multi_dir, "files",
                                         "myeloid_AUCell.csv"))

saveRDS(myeloid_seurat, file = file.path(no_multi_dir, "rda_obj",
                                      "seurat_no_mult_myeloid.rds"))
