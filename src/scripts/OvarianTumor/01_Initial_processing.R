library(tidyverse)
library(Seurat)

ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

# Variables
samples <- c("A_955_OvarianTumor", "D_GTFB_1230_OvarianTumor")
save_sample <- "OvarianTumor"
ffpe <- TRUE # Setting to TRUE will assume a probe set was used.

vars.to.regress <- NULL

# Set data dirs
base_dir <- "/Users/wellskr/Documents/Analysis/Benjamin_Bitler/visium_ovary"

source(file.path(base_dir, "src", "scripts", "functions.R"))

results_dir <- file.path(base_dir, "results")

save_dir_one <- file.path(results_dir, samples[[1]], "R_analysis")

save_dir_two <- file.path(results_dir, samples[[2]], "R_analysis")

save_dir <- file.path(results_dir, save_sample, "R_analysis")

seurat_one <- readRDS(file.path(save_dir_one, "rda_obj",
                                "seurat_processed.rds"))

seurat_two <- readRDS(file.path(save_dir_two, "rda_obj",
                                "seurat_processed.rds"))

# Make output directories
ifelse(!dir.exists(file.path(results_dir, save_sample)),
       dir.create(file.path(results_dir, save_sample)), FALSE)

ifelse(!dir.exists(file.path(save_dir)),
       dir.create(file.path(save_dir)), FALSE)

ifelse(!dir.exists(file.path(save_dir, "images")),
       dir.create(file.path(save_dir, "images")), FALSE)

ifelse(!dir.exists(file.path(save_dir, "files")),
       dir.create(file.path(save_dir, "files")), FALSE)

ifelse(!dir.exists(file.path(save_dir, "rda_obj")),
       dir.create(file.path(save_dir, "rda_obj")), FALSE)

# Merge files
seurat_obj <- merge(seurat_one, seurat_two)

DefaultAssay(seurat_obj) <- "Spatial"
DefaultAssay(seurat_one) <- "Spatial"
DefaultAssay(seurat_two) <- "Spatial"

VariableFeatures(seurat_obj) <- intersect(VariableFeatures(seurat_one),
                                         VariableFeatures(seurat_two))

DefaultAssay(seurat_obj) <- "SCT"
DefaultAssay(seurat_one) <- "SCT"
DefaultAssay(seurat_two) <- "SCT"

VariableFeatures(seurat_obj) <- intersect(VariableFeatures(seurat_one),
                                         VariableFeatures(seurat_two))
#SpatiallyVariableFeatures(seurat_obj) <- 
#  unique(c(SpatiallyVariableFeatures(seurat_one),
#           SpatiallyVariableFeatures(seurat_two)))

# Integrate files
#object_list <- list(seurat_one, seurat_two)

#spatial_features <- SelectIntegrationFeatures(object.list = object_list,
#                                              nfeatures = 3000)

#object_list <- PrepSCTIntegration(object.list = object_list,
#                                  anchor.features = spatial_features)

#spatial_anchors <- FindIntegrationAnchors(object.list = object_list,
#                                          normalization.method = "SCT",
#                                          anchor.features = spatial_features)

#seurat_obj <- IntegrateData(anchorset = spatial_anchors,
#                            normalization.method = "SCT")

saveRDS(seurat_obj, file = file.path(save_dir, "rda_obj", "seurat_start.rds"))
