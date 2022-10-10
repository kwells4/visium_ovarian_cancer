library(ggplot2)
library(Seurat)
library(tidyr)
library(cowplot)
library(dplyr)
library(clustifyr)
library(SingleCellExperiment)
library(BayesSpace)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

sample <- "B_1180_Omentum"

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


# Add in position info
position_info <- read.table(file.path(base_dir_proj,
                                      "outs/spatial/tissue_positions_list.csv"),
                            sep = ",", row.names = 1)

colnames(position_info) <- c("in_tissue", "row", "col",
                             "imagerow", "imagecol")

# Read in the data
seurat_data <-  readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))

seurat_data <- AddMetaData(seurat_data, position_info)

# Convert to sce
seurat_sce <- as.SingleCellExperiment(seurat_data)

# Add pca
seurat_pca <- Embeddings(seurat_data, reduction = "sctpca")

npcs <- 12

seurat_pca <- seurat_pca[ , 1:npcs]

colnames(seurat_pca) <- sub("_",  "", colnames(seurat_pca))

reducedDim(seurat_sce, "PCA") <- seurat_pca

# Set up for bayes space
set.seed(102)
seurat_sce <- spatialPreprocess(seurat_sce, platform="Visium", 
                              n.PCs=7, n.HVGs=2000, log.normalize=FALSE,
                              skip.PCA = TRUE)

# Determine number of clusters
seurat_sce <- qTune(seurat_sce, qs=seq(2, 10), platform="Visium", d=npcs)
qPlot(seurat_sce)

clusters <- 6

set.seed(149)
seurat_sce <- spatialCluster(seurat_sce, q=clusters, platform="Visium", d=npcs,
                           init.method="mclust", model="t", gamma=2,
                           nrep=10000, burn.in=100,
                           save.chain=TRUE)

clusterPlot(seurat_sce, palette = RColorBrewer::brewer.pal(name = "Set1",
                                                           n = clusters))


seurat_enhanced <- spatialEnhance(seurat_sce, q=clusters, platform="Visium", d=npcs,
                                    model="t", gamma=2,
                                    jitter_prior=0.3, jitter_scale=3.5,
                                    nrep=1000, burn.in=100,
                                    save.chain=TRUE)

clusterPlot(seurat_enhanced, palette = RColorBrewer::brewer.pal(name = "Set1",
                                                           n = clusters))

spot_data <- colData(seurat_enhanced)

data_info <- spot_data %>%
  data.frame %>%
  dplyr::count(spot.idx, spatial.cluster)

# any times the "n" is 6, all subspots agree on cluster
agree_spots <- data_info %>%
  dplyr::filter(n == 6) %>%
  dplyr::mutate(high_res_cluster = spatial.cluster)

diff_spots <- data_info %>%
  dplyr::filter(n < 6) %>%
  dplyr::mutate(high_res_cluster = "multiple")

diff_spots_new <- diff_spots %>%
  dplyr::distinct(spot.idx, high_res_cluster)

agree_spots_new <- agree_spots %>%
  dplyr::select(spot.idx, high_res_cluster)

all_spots <- rbind(agree_spots_new, diff_spots_new)

all_spots <- all_spots[order(all_spots$spot.idx),]

rownames(all_spots) <- rownames(seurat_data[[]])

seurat_data <- AddMetaData(seurat_data, all_spots)

SpatialDimPlot(seurat_data, group.by = "high_res_cluster",
               label.size = 3, cols = "Set1")


table(paste0(seurat_data$high_res_cluster, "_", seurat_data$SCT_celltype))

# Compare spatial cluster to seurat
spatial_data <- colData(seurat_sce) %>%
  data.frame %>%
  select(spatial.cluster)

seurat_data <- AddMetaData(seurat_data, spatial_data)

SpatialDimPlot(seurat_data, group.by = "spatial.cluster",
               label.size = 3, cols = "Set1")


table(paste0(seurat_data$spatial.cluster, "_", seurat_data$SCT_celltype))

plotDimRed(seurat_data, col_by = "high_res_cluster", plot_type = "rnasct.umap")


seurat_data$spatial.cluster <- factor(seurat_data$spatial.cluster)
plotDimRed(seurat_data, col_by = "spatial.cluster", plot_type = "rnasct.umap")

# Checking gene expression

# T cells - All are pretty similar, monocytes also high in seurat
plot(featDistPlot(seurat_data, geneset = "CD3D", sep_by = "high_res_cluster"))
plot(featDistPlot(seurat_data, geneset = "CD3D", sep_by = "spatial.cluster"))
plot(featDistPlot(seurat_data, geneset = "CD3D", sep_by = "SCT_cluster"))
plot(featDistPlot(seurat_data, geneset = "CD3D", sep_by = "SCT_celltype"))

# Monocytes - the spatial looks best here although seurat based on clusters is fine
plot(featDistPlot(seurat_data, geneset = "CD14", sep_by = "high_res_cluster"))
plot(featDistPlot(seurat_data, geneset = "CD14", sep_by = "spatial.cluster"))
plot(featDistPlot(seurat_data, geneset = "CD14", sep_by = "SCT_cluster"))
plot(featDistPlot(seurat_data, geneset = "CD14", sep_by = "SCT_celltype"))

# Fibroblasts - spatial looks better 
plot(featDistPlot(seurat_data, geneset = "COL1A1", sep_by = "high_res_cluster"))
plot(featDistPlot(seurat_data, geneset = "COL1A1", sep_by = "spatial.cluster"))
plot(featDistPlot(seurat_data, geneset = "COL1A1", sep_by = "SCT_cluster"))
plot(featDistPlot(seurat_data, geneset = "COL1A1", sep_by = "SCT_celltype"))

# Epithelial cells - None look good
plot(featDistPlot(seurat_data, geneset = "CDH1", sep_by = "high_res_cluster"))
plot(featDistPlot(seurat_data, geneset = "CDH1", sep_by = "spatial.cluster"))
plot(featDistPlot(seurat_data, geneset = "CDH1", sep_by = "SCT_cluster"))
plot(featDistPlot(seurat_data, geneset = "CDH1", sep_by = "SCT_celltype"))


## Go further once Ben sends marker genes
saveRDS(seurat_data, file.path(save_dir, "rda_obj", "seurat_processed.rds"))
saveRDS(seurat_enhanced, file.path(save_dir, "rda_obj", "enhanced_clusters.rds"))
saveRDS(seurat_sce, file.path(save_dir, "rda_obj", "sce_object.rds"))