library(ggplot2)
library(Seurat)
library(tidyr)
library(cowplot)
library(dplyr)
library(clustifyr)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

sample <- "A_955_OvarianTumor"

normalization_method <- "log" # can be SCT or log

HTO <- FALSE
ADT <- FALSE

if(normalization_method == "SCT"){
  SCT <- TRUE
  seurat_assay <- "SCT"
  clusters <- "SCT_cluster"
} else {
  SCT <- FALSE
  seurat_assay <- "Spatial"
  clusters <- "SCT_cluster"
}

# Set data dirs
base_dir <- "/Users/wellskr/Documents/Analysis/Benjamin_Bitler/visium_ovary"

source(file.path(base_dir, "src", "scripts", "functions.R"))

base_dir_proj <- file.path(base_dir, "results", sample)

save_dir <- file.path(base_dir_proj, "R_analysis")

# Read in the data
seurat_data <-  readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))
DefaultAssay(seurat_data) <- seurat_assay

# Published data ---------------------------------------------------------------
# Information for cell mapping
data_dir <-
  "/Users/wellskr/Documents/Analysis/Benjamin_Bitler/visium_ovary/resources"

ref_mat <- read.csv(file.path(data_dir,
                           "average_cluster_expression_my_clusters.csv"),
                    header = TRUE, row.names = 1)
# Run with SCT
seurat_res_list <- name_clusters(seurat_object = seurat_data,
                                 ref_mat = ref_mat,
                                 save_name = "celltype",
                                 ADT = ADT,
                                 nfeatures = 500,
                                 clusters = clusters,
                                 assay = "SCT",
                                 plot_type = "rnasct.umap")

seurat_data <- seurat_res_list$object

plotDimRed(seurat_data, "SCT_celltype", plot_type = "rnasct.umap")

# Repeat with Spatial
seurat_res_list <- name_clusters(seurat_object = seurat_data,
                                 ref_mat = ref_mat,
                                 save_name = "celltype",
                                 ADT = ADT,
                                 nfeatures = 500,
                                 clusters = clusters,
                                 assay = "Spatial",
                                 plot_type = "rnasct.umap")

seurat_data <- seurat_res_list$object

plotDimRed(seurat_data, "Spatial_celltype", plot_type = "rnasct.umap")


# Pan cancer -------------------------------------------------------------------
# Information for cell mapping
data_dir <-
  "/Users/wellskr/Documents/Analysis/Benjamin_Bitler/visium_ovary/resources/pan_cancer"

ref_mat <- read.csv(file.path(data_dir,
                              "average_cluster_expression_clusters.csv"),
                    header = TRUE, row.names = 1)
# Run with SCT
seurat_res_list <- name_clusters(seurat_object = seurat_data,
                                 ref_mat = ref_mat,
                                 save_name = "celltype_pc",
                                 ADT = ADT,
                                 nfeatures = 1000,
                                 clusters = clusters,
                                 assay = "SCT",
                                 plot_type = "rnasct.umap")

seurat_data <- seurat_res_list$object

plotDimRed(seurat_data, "SCT_celltype_pc", plot_type = "rnasct.umap")

# Repeat with Spatial
seurat_res_list <- name_clusters(seurat_object = seurat_data,
                                 ref_mat = ref_mat,
                                 save_name = "celltype_pc",
                                 ADT = ADT,
                                 nfeatures = 1000,
                                 clusters = clusters,
                                 assay = "Spatial",
                                 plot_type = "rnasct.umap")

seurat_data <- seurat_res_list$object

plotDimRed(seurat_data, "Spatial_celltype_pc", plot_type = "rnasct.umap")


saveRDS(seurat_data, file.path(save_dir, "rda_obj", "seurat_processed.rds"))
