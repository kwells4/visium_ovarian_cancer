library(tidyverse)
library(Seurat)
library(clustifyr)

ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

# Variables
sample <- "pan_cancer_ref"
ffpe <- TRUE # Setting to TRUE will assume a probe set was used.

vars.to.regress <- NULL

# Set data dirs
base_dir <- "/Users/wellskr/Documents/Analysis/Benjamin_Bitler/visium_ovary/"

source(file.path(base_dir, "src", "scripts", "functions.R"))

count_path <- file.path(base_dir, "resources/pan_cancer")

seurat_obj <- create_seurat_object(count_path = count_path,
                                   sample = sample,
                                   ADT = FALSE,
                                   hashtag = FALSE)


meta_data <- read.csv(file.path(count_path, "2101-Ovariancancer_metadata.csv"),
                      row.names = 1)

# Subset to only the cells in the meta data
seurat_obj <- subset(seurat_obj, cells = rownames(meta_data))

# Add meta data to object
seurat_obj <- AddMetaData(seurat_obj, metadata = meta_data)

seurat_obj <- NormalizeData(seurat_obj) %>% 
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(npcs = 50) %>%
  FindNeighbors(dims = 1:50) %>%
  RunUMAP(dims = 1:50) %>%
  FindClusters(resolution = 0.3)

seurat_obj$CellFromTumorWord <- ifelse(seurat_obj$CellFromTumor == 1,
                                       "Tumor", "NotTumor")

seurat_obj$tumorsite_celltype <- paste0(seurat_obj$TumorSite, "_",
                                    seurat_obj$CellType)

seurat_obj$tumorsite_from_tumor <- paste0(seurat_obj$TumorSite, "_",
                                          seurat_obj$CellFromTumorWord)

seurat_obj$celltype_from_tumor <- paste0(seurat_obj$CellType, "_",
                                         seurat_obj$CellFromTumorWord)

seurat_obj$all_info <- paste0(seurat_obj$CellType, "_",
                              seurat_obj$TumorSite, "_",
                              seurat_obj$CellFromTumorWord)

# Make reference ---------------------------------------------------------------
ref_counts <- GetAssayData(seurat_obj, "data")

ref_metadata <- seurat_obj[[]]

## CellType --------------------------------------------------------------------
meta_cluster <- "CellType"

ref_mat_celltype <- average_clusters(
  mat = ref_counts,
  metadata = ref_metadata,
  cluster_col = meta_cluster,
  if_log = TRUE
)

write.csv(ref_mat_celltype,
          file.path(count_path, "average_cluster_expression_celltype.csv"))

## CellType tumor site ---------------------------------------------------------
meta_cluster <- "tumorsite_celltype"

ref_mat_celltype_tumorsite <- average_clusters(
  mat = ref_counts,
  metadata = ref_metadata,
  cluster_col = meta_cluster,
  if_log = TRUE
)

write.csv(ref_mat_celltype_tumorsite,
          file.path(count_path,
                    "average_cluster_expression_celltype_tumorsite.csv"))

## CellType From Tumor ---------------------------------------------------------
meta_cluster <- "celltype_from_tumor"

ref_mat_celltype_from_tumor <- average_clusters(
  mat = ref_counts,
  metadata = ref_metadata,
  cluster_col = meta_cluster,
  if_log = TRUE
)

write.csv(ref_mat_celltype_from_tumor,
          file.path(count_path,
                    "average_cluster_expression_celltype_from_tumor.csv"))

## All Info --------------------------------------------------------------------
meta_cluster <- "all_info"

ref_mat_all_info <- average_clusters(
  mat = ref_counts,
  metadata = ref_metadata,
  cluster_col = meta_cluster,
  if_log = TRUE
)

write.csv(ref_mat_all_info,
          file.path(count_path,
                    "average_cluster_expression_all_info.csv"))

saveRDS(seurat_obj, file = file.path(count_path, "seurat_obj.rds"))


seurat_obj$PatientNumber <- factor(seurat_obj$PatientNumber)

plotDimRed(seurat_obj, col_by = "celltype_from_tumor", plot_type = "umap")
plotDimRed(seurat_obj, col_by = "TumorType", plot_type = "umap")
plotDimRed(seurat_obj, col_by = "CellType", plot_type = "umap")
plotDimRed(seurat_obj, col_by = "CellFromTumorWord", plot_type = "umap")
plotDimRed(seurat_obj, col_by = "PatientNumber", plot_type = "umap")
plotDimRed(seurat_obj, col_by = "seurat_clusters", plot_type = "umap")
plotDimRed(seurat_obj, col_by = "CD3D", plot_type = "umap")
plotDimRed(seurat_obj, col_by = "CD8A", plot_type = "umap")

seurat_obj$high_res_cluster <- paste0(seurat_obj$CellType, "_",
                                      seurat_obj$seurat_clusters)


## clusters --------------------------------------------------------------------
meta_cluster <- "high_res_cluster"
ref_metadata <- seurat_obj[[]]

ref_mat_cluster <- average_clusters(
  mat = ref_counts,
  metadata = ref_metadata,
  cluster_col = meta_cluster,
  if_log = TRUE
)

# Accurate clusters have more than 10 cells
cluster_counts <- data.frame(table(seurat_obj$high_res_cluster)) %>%
  dplyr::filter(Freq > 20)

ref_mat_cluster <- ref_mat_cluster %>%
  data.frame %>%
  select(cluster_counts$Var1)

write.csv(ref_mat_cluster,
          file.path(count_path,
                    "average_cluster_expression_clusters.csv"))


## Markers ---------------------------------------------------------------------
Idents(seurat_obj) <- "high_res_cluster"
marker_list <- FindAllMarkers(seurat_obj, only.pos = TRUE)

marker_list <- marker_list %>%
  dplyr::filter(p_val_adj < 0.05)

write.csv(marker_list,
          file.path(count_path, "cluster_makers.csv"))

# T cell markers
t_markers <- marker_list[grepl("T_cell", marker_list$cluster),] %>%
  dplyr::filter(p_val_adj < 0.05)
