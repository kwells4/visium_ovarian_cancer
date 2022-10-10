library(ggplot2)
library(Seurat)
library(tidyr)
library(cowplot)
library(dplyr)
library(clustifyr)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

sample <- "D_GTFB_1230_OvarianTumor"

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

save_dir <- file.path(base_dir_proj, "R_analysis_no_multi")

# Read in the data
seurat_data <-  readRDS(file.path(save_dir, "rda_obj",
                                  "seurat_no_mult_processed.rds"))
DefaultAssay(seurat_data) <- seurat_assay
# Information for cell mapping
data_dir <-
  "/Users/wellskr/Documents/Analysis/Benjamin_Bitler/visium_ovary/resources"

ref_mat_1 <- read.csv(file.path(data_dir,
                              "average_cluster_expression_my_clusters.csv"),
                    header = TRUE, row.names = 1)
# Run with SCT
seurat_res_list <- name_clusters(seurat_object = seurat_data,
                                 ref_mat = ref_mat_1,
                                 save_name = "celltype_no_mult",
                                 ADT = ADT,
                                 nfeatures = 500,
                                 clusters = clusters,
                                 assay = "SCT",
                                 plot_type = "rnasct.umap")

seurat_data <- seurat_res_list$object

plotDimRed(seurat_data, "SCT_celltype_no_mult", plot_type = "rnasct.umap")
plotDimRed(seurat_data, "SCT_celltype", plot_type = "rnasct.umap")
# Repeat with Spatial
seurat_res_list <- name_clusters(seurat_object = seurat_data,
                                 ref_mat = ref_mat_1,
                                 save_name = "celltype_no_mult",
                                 ADT = ADT,
                                 nfeatures = 500,
                                 clusters = clusters,
                                 assay = "Spatial",
                                 plot_type = "rnasct.umap")

seurat_data <- seurat_res_list$object

plotDimRed(seurat_data, "Spatial_celltype_no_mult", plot_type = "rnasct.umap")

# Pan cancer -------------------------------------------------------------------
# Information for cell mapping
data_dir <-
  "/Users/wellskr/Documents/Analysis/Benjamin_Bitler/visium_ovary/resources/pan_cancer"

ref_mat_2 <- read.csv(file.path(data_dir,
                              "average_cluster_expression_all_info.csv"),
                    header = TRUE, row.names = 1)
# Run with SCT
seurat_res_list <- name_clusters(seurat_object = seurat_data,
                                 ref_mat = ref_mat_2,
                                 save_name = "celltype_pc",
                                 ADT = ADT,
                                 nfeatures = 500,
                                 clusters = clusters,
                                 assay = "SCT",
                                 plot_type = "rnasct.umap")

seurat_data <- seurat_res_list$object

plotDimRed(seurat_data, "SCT_celltype_pc", plot_type = "rnasct.umap")

# Repeat with Spatial
seurat_res_list <- name_clusters(seurat_object = seurat_data,
                                 ref_mat = ref_mat_2,
                                 save_name = "celltype_pc",
                                 ADT = ADT,
                                 nfeatures = 500,
                                 clusters = clusters,
                                 assay = "Spatial",
                                 plot_type = "rnasct.umap")

seurat_data <- seurat_res_list$object

plotDimRed(seurat_data, "Spatial_celltype_pc", plot_type = "rnasct.umap")


# combined ---------------------------------------------------------------------
# Information for cell mapping
combined_ref <- merge(ref_mat_1, ref_mat_2, by = "row.names",
                      all.x = FALSE, all.y = FALSE)

rownames(combined_ref) <- combined_ref$Row.names

combined_ref$Row.names <- NULL

# Run with SCT
seurat_res_list <- name_clusters(seurat_object = seurat_data,
                                 ref_mat = combined_ref,
                                 save_name = "celltype_combined",
                                 ADT = ADT,
                                 nfeatures = 500,
                                 clusters = clusters,
                                 assay = "SCT",
                                 plot_type = "rnasct.umap")

seurat_data <- seurat_res_list$object

plotDimRed(seurat_data, "SCT_celltype", plot_type = "rnasct.umap")
plotDimRed(seurat_data, "SCT_celltype_combined", plot_type = "rnasct.umap")

# Repeat with Spatial
seurat_res_list <- name_clusters(seurat_object = seurat_data,
                                 ref_mat = combined_ref,
                                 save_name = "celltype_combined",
                                 ADT = ADT,
                                 nfeatures = 500,
                                 clusters = clusters,
                                 assay = "Spatial",
                                 plot_type = "rnasct.umap")

seurat_data <- seurat_res_list$object

plotDimRed(seurat_data, "SCT_celltype", plot_type = "rnasct.umap")
plotDimRed(seurat_data, "Spatial_celltype_combined", plot_type = "rnasct.umap")


saveRDS(seurat_data, file.path(save_dir, "rda_obj",
                               "seurat_no_mult_processed.rds"))

gene_list <- c("EHMT1",
               "EHMT2",
               "CLDN4",
               "ATF6",
               "IL6",
               #"GOLM2", #CASC4
               "CXCL10",
               "DUSP1",
               "PTP4A3",
               "PAX8",
               "CRABP2",
               "MUC16",
               "POLQ",
               "BRCA1",
               "BRCA2",
               "TMEM173",
               "CBX2",
               "CPT1A",
               "IER3",
               "WNT4",
               "WNT3A",
               "FZD2",
               "CTNNB1",
               #"T53BP1", 
               "XRCC1",
               "MAPK14",
               "EGFR",
               "ERBB3",
               "OLFML3",
               "RB1",
               "NOTCH3",
               "MUC1",
               "THBS1",
               "RXRA",
               "IGFBP2")

plots <- plotDimRed(seurat_data,
                    gene_list, plot_type = "rnasct.umap")
