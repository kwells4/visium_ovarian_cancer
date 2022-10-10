library(ggplot2)
library(Seurat)
library(tidyr)
library(cowplot)
library(dplyr)
library(openxlsx)
library(STutility)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

sample <- "A_955_OvarianTumor"
cell_types <- "SCT_celltype_combined"
celltype2 <- "Spatial_celltype_combined"
clusters <- "SCT_cluster"

DE_assay <- "Spatial" # Keep spatial or RNA for now

# Set directories
base_dir <- "/Users/wellskr/Documents/Analysis/Benjamin_Bitler/visium_ovary"


normalization_method <- "SCT" # can be SCT or log

HTO <- FALSE
ADT <- FALSE # Set to true if you want to run DE on ADT (not enough ADT here)
spatial <- TRUE
pval <- 0.05
logfc <- 0.5

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

save_dir_full <- file.path(base_dir_proj, "R_analysis")

# Read in the data
seurat_data <- readRDS(file.path(save_dir, "rda_obj",
                                 "seurat_no_mult_processed.rds"))

# Read in the data (not subset)
seurat_data_full <- readRDS(file.path(save_dir_full, "rda_obj",
                                      "seurat_processed.rds"))

seurat_data_full <- AddMetaData(seurat_data_full,
                                seurat_data$SCT_celltype_combined,
                                col.name = "no_multi_celltype")


# Add celltypes from subset to full
seurat_data_full$no_multi_celltype[is.na(seurat_data_full$no_multi_celltype)] <-
  "multiple_celltypes"

SpatialDimPlot(seurat_data_full, group.by = "no_multi_celltype")

SpatialDimPlot(seurat_data, group.by = "SCT_celltype_combined")

# Read in data for STutility
info_table <- read.table(file.path(base_dir, "Files",
                                   paste0(sample, "_infoTable.csv")),
                         sep = ",", header = TRUE)

se <- InputFromTable(infotable = info_table,
                     platform =  "Visium")

# Add cell types
metadata <- seurat_data_full[[]]

rownames(metadata) <- paste0(rownames(metadata), "_1")

se <- AddMetaData(se, metadata = metadata)

se <- SetIdent(se, value = "no_multi_celltype")


se <- LoadImages(se, time.resolve = FALSE, verbose = TRUE)

# Find neighbors
se <- RegionNeighbours(se,
                       id = "Myeloid_Ovarium_Tumor",
                       verbose = TRUE,
                       keep.within.id = FALSE)


FeatureOverlay(se, features = "no_multi_celltype", sampleids = 1)

FeatureOverlay(se, features = "SCT_celltype", sampleids = 1)


FeatureOverlay(se, features = "nbs_Myeloid_Ovarium_Tumor", sampleids = 1)

# Add neighbors back to original seurat object
nbs_meta <- se[[]] %>%
  dplyr::select(nbs_Myeloid_Ovarium_Tumor)

rownames(nbs_meta) <- sub("_1", "", rownames(nbs_meta))

seurat_data_full <- AddMetaData(seurat_data_full,
                                nbs_meta)

Idents(seurat_data_full) <- "nbs_Myeloid_Ovarium_Tumor"

# Find markers
nbs_myeloid.markers <- FindMarkers(seurat_data_full,
                             ident.1 = "Myeloid_Ovarium_Tumor",
                             ident.2 = "nbs_Myeloid_Ovarium_Tumor")

nbs_myeloid.markers %>%
  dplyr::filter(p_val_adj < 0.05)

SpatialDimPlot(seurat_data_full, group.by = "no_multi_celltype",
               image.alpha = 0, cols = "Set1")

SpatialDimPlot(seurat_data_full, group.by = "nbs_Myeloid_Ovarium_Tumor",
               image.alpha = 0, cols = "Set1")


seurat_data_full$neighbors_full <- paste0(seurat_data_full$no_multi_celltype,
                                          "_",
                                          seurat_data_full$nbs_Myeloid_Ovarium_Tumor)

seurat_data_full$neighbors_full[is.na(seurat_data_full$nbs_Myeloid_Ovarium_Tumor)] <- NA

SpatialDimPlot(seurat_data_full, group.by = "neighbors_full",
               image.alpha = 0, cols = "Set1")
