library(ggplot2)
library(Seurat)
library(tidyr)
library(cowplot)
library(dplyr)
library(openxlsx)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

sample <- "B_1180_Omentum"
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

gene_lists <- NULL

# Read in the data
seurat_data <- readRDS(file.path(save_dir, "rda_obj",
                                 "seurat_no_mult_processed.rds"))

# Cell type DE -----------------------------------------------------------------

marker_list <- find_write_markers(seurat_object = seurat_data,
                                  meta_col = cell_types,
                                  pval = pval,
                                  logfc = logfc,
                                  assay = DE_assay,
                                  save_dir = save_dir)

if(ADT){
  marker_list <- find_write_markers(seurat_object = seurat_data,
                                    meta_col = cell_types,
                                    pval = pval,
                                    logfc = logfc,
                                    assay = "ADT",
                                    save_dir = save_dir)
}

# Cell type DE 2 ---------------------------------------------------------------

marker_list <- find_write_markers(seurat_object = seurat_data,
                                  meta_col = celltype2,
                                  pval = pval,
                                  logfc = logfc,
                                  assay = DE_assay,
                                  save_dir = save_dir)

if(ADT){
  marker_list <- find_write_markers(seurat_object = seurat_data,
                                    meta_col = celltype2,
                                    pval = pval,
                                    logfc = logfc,
                                    assay = "ADT",
                                    save_dir = save_dir)
}

# RNA cluster DE ---------------------------------------------------------------

marker_list <- find_write_markers(seurat_object = seurat_data,
                                  meta_col = clusters,
                                  pval = pval,
                                  logfc = logfc,
                                  assay = DE_assay,
                                  save_dir = save_dir)

if(ADT){
  marker_list <- find_write_markers(seurat_object = seurat_data,
                                    meta_col = clusters,
                                    pval = pval,
                                    logfc = logfc,
                                    assay = "ADT",
                                    save_dir = save_dir)
}

# Can also subset and get cells that way --> ex find DE between T cells in
# Different locations. First subset to T cells, then highlight those areas
# to get the area specific T cells, then do DE.
#cortex <- subset(seurat_data,
#                 slice1_imagerow > 200 | slice1_imagecol < 300,
#                 invert = TRUE)

#p2 <- SpatialDimPlot(cortex, crop = FALSE, label = TRUE,
#                     pt.size.factor = 1, label.size = 3)


saveRDS(seurat_data, file.path(save_dir, "rda_obj",
                               "seurat_no_mult_processed.rds"))
