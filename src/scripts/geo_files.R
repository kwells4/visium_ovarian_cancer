library(tidyverse)
library(Seurat)
library(here)

ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

source(here("src", "scripts", "functions.R"))


# All samples
samples <- c("A_955_OvarianTumor", "A_GTFB1154_OvarianCancerTumor",
             "B_1180_Omentum", "B_GTFB1191_OvarianCancerTumor",
             "C_GTFB1170_SmallCellOvarianCancer", "C_GTFB_1191_OmentumTumor",
             "D_GTFB1170_SmallCellOvarianCancer", "D_GTFB_1230_OvarianTumor")

save_dir_all <- here("results", "geo_files")
ifelse(!dir.exists(save_dir_all), dir.create(save_dir_all), FALSE)

save_info <- function(sample_name){
  base_dir_proj <- here("results", sample_name)
  
  save_dir <- file.path(base_dir_proj, "R_analysis")
  
  spatial_dir <- file.path(base_dir_proj, "outs", "spatial")
  
  new_save_dir <- file.path(save_dir_all, sample_name)
  ifelse(!dir.exists(new_save_dir), dir.create(new_save_dir), FALSE)
  
  # Copy all spatial files
  spatial_files <- file.path(spatial_dir, list.files(spatial_dir))
  
  new_files <- file.path(new_save_dir, list.files(spatial_dir))
  
  file.copy(spatial_files, new_files)
  
  # Read in the data
  seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds")) 
  
  # Pull out sct normalized counts
  DefaultAssay(seurat_data) <- "SCT"
  
  assay_data <- GetAssayData(seurat_data) %>%
    data.frame()
  
  write.csv(assay_data,
            file.path(save_dir_all, sample_name, "sct_normalized_counts.csv"))
  
  
  assay_data_raw <- GetAssayData(seurat_data, slot = "counts") %>%
    data.frame()
  
  write.csv(assay_data_raw,
            file.path(save_dir_all, sample_name, "raw_counts.csv"))
  
  
  # Pull out meta data
  meta_data <- seurat_data[[]] %>%
    dplyr::select(orig.ident, nCount_Spatial, nFeature_Spatial,
                  S.Score, G2M.Score, Phase, nCount_SCT, nFeature_SCT,
                  SCT_cluster, SCT_celltype, in_tissue, row, col, imagerow,
                  imagecol, spot.idx, high_res_cluster)

  write.csv(meta_data,
            file.path(save_dir_all, sample_name, "metadata.csv")) 
  
  umap_coords <- Embeddings(seurat_data, reduction = "rnasct.umap")
  colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
  
  write.csv(umap_coords,
            file.path(save_dir_all, sample_name, "umap_coords.csv")) 
}

invisible(lapply(samples, save_info))

# Once this is made, tar all files