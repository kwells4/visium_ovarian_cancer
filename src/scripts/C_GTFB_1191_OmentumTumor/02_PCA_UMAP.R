library(tidyverse)
library(Seurat)

ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

# Variables
sample <- "C_GTFB_1191_OmentumTumor"
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

pdf(file.path(save_dir, "images", "RNA_pca.pdf"))
print(RNA_plots)
dev.off()

# UMAP -------------------------------------------------------------------------
RNA_pcs <- 19

set.seed(0)
# UMAP of gene expression
umap_data <- group_cells(seurat_data, sample, save_dir, nPCs = RNA_pcs,
                         resolution = 0.6, assay = seurat_assay, HTO = FALSE)

seurat_data <- umap_data$object

gene_plots <- umap_data$plots

#DefaultAssay(seurat_data) <- "Spatial"

#plotDimRed(seurat_data, col_by = "CD4", plot_type = "rnasct.umap")
#plotDimRed(seurat_data, col_by = "CD8A", plot_type = "rnasct.umap")
#plotDimRed(seurat_data, col_by = "CD3D", plot_type = "rnasct.umap")
#plotDimRed(seurat_data, col_by = "ITLN1", plot_type = "rnasct.umap")
#plotDimRed(seurat_data, col_by = "IGHG2", plot_type = "rnasct.umap")
#plotDimRed(seurat_data, col_by = "PRDM1", plot_type = "rnasct.umap")
#plotDimRed(seurat_data, col_by = "SUGCT", plot_type = "rnasct.umap")
#plotDimRed(seurat_data, col_by = "BAMBI", plot_type = "rnasct.umap")
#plotDimRed(seurat_data, col_by = "CCL21", plot_type = "rnasct.umap")
#plotDimRed(seurat_data, col_by = "MKI67", plot_type = "rnasct.umap")


SpatialDimPlot(seurat_data, label = TRUE, label.size = 3, cols = "Set1")
#SpatialFeaturePlot(seurat_data, features = "CD3D", pt.size.factor = 1)
#SpatialFeaturePlot(seurat_data, features = "CD3D", pt.size.factor = 1,
#                   alpha = c(0.1, 1))

# Highlight the T cells
#SpatialDimPlot(seurat_data,
#               cells.highlight = CellsByIdentities(object = seurat_data,
#                                                   idents = 3))

saveRDS(seurat_data, file = file.path(save_dir, "rda_obj", "seurat_processed.rds"))

