library(ggplot2)
library(Seurat)
library(tidyr)
library(cowplot)
library(dplyr)
library(openxlsx)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

sample <- "C_GTFB_1191_OmentumTumor"
cell_types <- "SCT_celltype_combined"
celltype2 <- "Spatial_celltype_combined"
clusters <- "SCT_cluster"


spatial <- TRUE

if(normalization_method == "SCT"){
  SCT <- TRUE
  seurat_assay <- "SCT"
} else {
  SCT <- FALSE
  seurat_assay <- "RNA"
}

source(file.path(base_dir, "src", "scripts", "functions.R"))

base_dir_proj <- file.path(base_dir, "results", sample)

save_dir <- file.path(base_dir_proj, "R_analysis_no_multi")

# Read in the data
seurat_data <- readRDS(file.path(save_dir, "rda_obj",
                                 "seurat_no_mult_processed.rds"))

gene <- "HP"

plotDimRed(seurat_data, plot_type = "rnasct.umap",
           col_by = gene)

plot(featDistPlot(seurat_data, geneset = gene))
