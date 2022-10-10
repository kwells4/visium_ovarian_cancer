library(tidyverse)
library(Seurat)
library(harmony)

ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

# Variables
sample <- "OvarianTumor"
ffpe <- TRUE # Setting to TRUE will assume a probe set was used.
seurat_assay <- "SCT"

# Set data dirs
base_dir <- "/Users/wellskr/Documents/Analysis/Benjamin_Bitler/visium_ovary"

source(file.path(base_dir, "src", "scripts", "functions.R"))

results_dir <- file.path(base_dir, "results")

save_dir <- file.path(results_dir, sample, "R_analysis")

# Read in the data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_start.rds"))

# PCA --------------------------------------------------------------------------

# PCA of gene expression
seurat_data <- PCA_dimRed(seurat_data, assay = seurat_assay)

RNA_plots <- plot_PCA(HTO = FALSE, assay = seurat_assay,
                      sample_object = seurat_data,
                      data_type = "spatial", ffpe = ffpe)

# Run harmony, when combining all, think about adding in slide information?
seurat_data <- RunHarmony(seurat_data, "orig.ident",
                          plot_convergence = TRUE,
                          theta = 10, 
                          assay.use = "SCT",
                          reduction = "sctpca",
                          nclust = 5,
                          reduction.save = "harmony_clust_5")

seurat_data <- RunHarmony(seurat_data, "orig.ident",
                          plot_convergence = TRUE, lambda = 20,
                          theta = 10, sigma = 0.1,
                          assay.use = "SCT",
                          reduction = "sctpca",
                          reduction.save = "harmony_high_lambda")

seurat_data <- RunHarmony(seurat_data, "orig.ident",
                          plot_convergence = TRUE, 
                          theta = 10, 
                          assay.use = "SCT",
                          reduction = "sctpca",
                          reduction.save = "harmony_default")


harmony_plot <- plotDimRed(seurat_data, col_by = "orig.ident",
                           plot_type = "harmony")


pdf(file.path(save_dir, "images", "RNA_pca.pdf"))
print(RNA_plots)
print(harmony_plot)
dev.off()

# UMAP -------------------------------------------------------------------------
RNA_pcs <- 30

set.seed(0)
# UMAP of gene expression
umap_data <- group_cells(seurat_data, sample, save_dir, nPCs = RNA_pcs,
                         resolution = 0.6, assay = seurat_assay, HTO = FALSE,
                         reduction = "harmony_clust_5")

seurat_data <- umap_data$object

gene_plots <- umap_data$plots

plotDimRed(seurat_data, col_by = "orig.ident",
           plot_type = "harmony_clust_5.umap")

plotDimRed(seurat_data, col_by = "SCT_celltype",
           plot_type = "harmony_clust_5.umap")

# UMAP of gene expression
umap_data <- group_cells(seurat_data, sample, save_dir, nPCs = RNA_pcs,
                         resolution = 0.6, assay = seurat_assay, HTO = FALSE,
                         reduction = "harmony_high_lambda")

seurat_data <- umap_data$object

gene_plots <- umap_data$plots

plotDimRed(seurat_data, col_by = "orig.ident",
           plot_type = "harmony_high_lambda.umap")

plotDimRed(seurat_data, col_by = "SCT_celltype",
           plot_type = "harmony_high_lambda.umap")

# UMAP of gene expression
umap_data <- group_cells(seurat_data, sample, save_dir, nPCs = RNA_pcs,
                         resolution = 0.6, assay = seurat_assay, HTO = FALSE,
                         reduction = "harmony_default")

seurat_data <- umap_data$object

gene_plots <- umap_data$plots

plotDimRed(seurat_data, col_by = "orig.ident",
           plot_type = "harmony_default.umap")

plotDimRed(seurat_data, col_by = "SCT_celltype",
           plot_type = "harmony_default.umap")


plotDimRed(seurat_data, col_by = "CD3D", plot_type = "harmony.umap")


plotDimRed(seurat_data, col_by = "SCT_celltype", plot_type = "harmony.umap")

SpatialDimPlot(seurat_data, label = TRUE, label.size = 3, cols = "Set1")
#SpatialFeaturePlot(seurat_data, features = "CD3D", pt.size.factor = 1)
#SpatialFeaturePlot(seurat_data, features = "CD3D", pt.size.factor = 1,
#                   alpha = c(0.1, 1))

# Highlight the T cells
SpatialDimPlot(seurat_data,
               cells.highlight = CellsByIdentities(object = seurat_data,
                                                   idents = 3))

saveRDS(seurat_data, file = file.path(save_dir, "rda_obj", "seurat_processed.rds"))

