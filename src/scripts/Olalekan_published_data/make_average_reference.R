library(Seurat)
library(tidyverse)
library(harmony)
library(clustifyr)

sample = "published_data"

ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

# Set data dirs
base_dir <- "/Users/wellskr/Documents/Analysis/Benjamin_Bitler/visium_ovary/"

source(file.path(base_dir, "src", "scripts", "functions.R"))

data_dir <-
  "/Users/wellskr/Documents/Analysis/Benjamin_Bitler/visium_ovary/resources"

seurat_merge <- readRDS(file.path(data_dir, "seurat_processed.rds"))

#Identify cell types by cluster
cm <- confusionMatrix(seurat_merge$seurat_clusters,
                      seurat_merge$assignedCelltype)

assignments <- colnames(cm)[apply(cm, 1 , which.max)]
assignments <- cbind(assignments, rownames(cm)) #Assignments

new_celltypes <- assignments[,1]
names(new_celltypes) <- assignments[,2]

# Make meta cluster with new assignments
seurat_merge$new_assignments <- new_celltypes[
  as.character(seurat_merge$seurat_clusters)]

plotDimRed(seurat_merge, "new_assignments", plot_type = "harmony.umap")

plotDimRed(seurat_merge, "seurat_clusters", plot_type = "harmony.umap")

plotDimRed(seurat_merge, "assignedCelltype", plot_type = "harmony.umap")

# Make average assignments

ref_counts <- GetAssayData(seurat_merge, "data")

ref_metadata <- seurat_merge[[]]

meta_cluster_pub <- "assignedCelltype"

meta_cluster_mine <- "new_assignments"

# Everything looks good so now to make the reference for clustifyr
ref_mat_pub <- average_clusters(
  mat = ref_counts,
  metadata = ref_metadata,
  cluster_col = meta_cluster_pub,
  if_log = TRUE
)

ref_mat_mine <- average_clusters(
  mat = ref_counts,
  metadata = ref_metadata,
  cluster_col = meta_cluster_pub,
  if_log = TRUE
)

write.csv(ref_mat_pub,
          file.path(data_dir, "average_cluster_expression_publisehd.csv"))

write.csv(ref_mat_mine,
          file.path(data_dir, "average_cluster_expression_my_clusters.csv"))
