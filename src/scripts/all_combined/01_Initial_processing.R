library(tidyverse)
library(Seurat)

ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

# Variables
samples <- c("A_955_OvarianTumor", "A_GTFB1154_OvarianCancerTumor",
             "B_1180_Omentum", "B_GTFB1191_OvarianCancerTumor",
             "C_GTFB1170_SmallCellOvarianCancer", "C_GTFB_1191_OmentumTumor",
             "D_GTFB1170_SmallCellOvarianCancer", "D_GTFB_1230_OvarianTumor")
save_sample <- "all_combined"
ffpe <- TRUE # Setting to TRUE will assume a probe set was used.

vars.to.regress <- NULL

# Set data dirs
base_dir <- "/Users/wellskr/Documents/Analysis/Benjamin_Bitler/visium_ovary"

source(file.path(base_dir, "src", "scripts", "functions.R"))

results_dir <- file.path(base_dir, "results")

save_dir <- file.path(results_dir, save_sample, "R_analysis")

seurat_obj_list <- lapply(samples, function(x){
        save_dir_obj <- file.path(results_dir, x, "R_analysis")
        seurat_obj <- readRDS(file.path(save_dir_obj, "rda_obj",
                                        "seurat_processed.rds"))
        return(seurat_obj)
        
})

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
seurat_obj <- merge(x = seurat_obj_list[[1]], y = seurat_obj_list[2:8],
                    merge.data = TRUE)

variable_features_spatial <- lapply(seurat_obj_list, function(x){
        DefaultAssay(x) <- "Spatial"
        return(VariableFeatures(x))
})

variable_features_spatial <- unique(unlist(variable_features_spatial))

variable_features_sct <- rownames(seurat_obj[["SCT"]]@scale.data)

DefaultAssay(seurat_obj) <- "Spatial"

VariableFeatures(seurat_obj) <- variable_features_spatial

DefaultAssay(seurat_obj) <- "SCT"


VariableFeatures(seurat_obj) <- variable_features_sct

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
