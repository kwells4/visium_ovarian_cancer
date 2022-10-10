library(ggplot2)
library(Seurat)
library(tidyr)
library(cowplot)
library(dplyr)
library(openxlsx)
library(here)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

sample <- "A_955_OvarianTumor"
cell_types <- "SCT_celltype_combined"
celltype2 <- "Spatial_celltype_combined"
clusters <- "SCT_cluster"

source(here("src/scripts/functions.R"))

base_dir_proj <- here(file.path("results", sample))

no_multi_dir <- file.path(base_dir_proj, "R_analysis_no_multi")
all_dir <- file.path(base_dir_proj, "R_analysis")

# Read in the data
seurat_data_no_multi <- readRDS(file.path(no_multi_dir, "rda_obj",
                                 "seurat_no_mult_processed.rds"))

# Read in the data (not subset)
seurat_data_full <- readRDS(file.path(all_dir, "rda_obj",
                                      "seurat_processed.rds"))

annotations <- read.csv(file.path(base_dir_proj, "outs/Annotated.csv"))

rownames(annotations) <- annotations$Barcode
annotations$Barcode <- NULL

seurat_data_no_multi <- AddMetaData(seurat_data_no_multi,
                                    annotations)

seurat_data_full <- AddMetaData(seurat_data_full,
                                annotations)

seurat_data_full$Annotated[seurat_data_full$Annotated == ""] <- "not_annotated"

seurat_data_no_multi$Annotated[seurat_data_no_multi$Annotated == ""] <- "not_annotated"


confusionMatrix(seurat_data_full$SCT_celltype, seurat_data_full$Annotated)

celltypes1 <- confusionMatrix(seurat_data_no_multi$SCT_celltype_pc,
                seurat_data_no_multi$Annotated)

celltypes2 <- confusionMatrix(seurat_data_no_multi$SCT_celltype_no_mult,
                seurat_data_no_multi$Annotated)

celltypes3 <- confusionMatrix(seurat_data_no_multi$SCT_celltype_combined,
                seurat_data_no_multi$Annotated)

cols <- RColorBrewer::brewer.pal(5, "Set1")
cols <- c("#D3D3D3", cols)

names(cols) <- unique(seurat_data_no_multi$Annotated)

pdf(file.path(no_multi_dir, "images", "annotated_plot.pdf"))
print(SpatialDimPlot(seurat_data_no_multi, group.by = "Annotated", cols = cols))
print(plotDimRed(seurat_data_no_multi, col_by = "Annotated", color = cols,
           plot_type = "rnasct.umap")[[1]])

dev.off()

confusion_matrix <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb = confusion_matrix, sheet = "celltype_pan_cancer")
openxlsx::writeData(wb = confusion_matrix, sheet = "celltype_pan_cancer",
                    x = data.frame(celltypes1), rowNames = TRUE)

openxlsx::addWorksheet(wb = confusion_matrix, sheet = "celltype_sc_ref")
openxlsx::writeData(wb = confusion_matrix, sheet = "celltype_sc_ref",
                    x = data.frame(celltypes2), rowNames = TRUE)


openxlsx::addWorksheet(wb = confusion_matrix, sheet = "celltype_combined_ref")
openxlsx::writeData(wb = confusion_matrix, sheet = "celltype_combined_ref",
                    x = data.frame(celltypes3), rowNames = TRUE)
                    
openxlsx::saveWorkbook(wb = confusion_matrix,
                       file = file.path(no_multi_dir, "files",
                                        "annotated_vs_reference_mapping.xlsx"),
                       overwrite = TRUE)


# Read in the data
saveRDS(seurat_data_no_multi, file.path(no_multi_dir, "rda_obj",
                                          "seurat_no_mult_processed.rds"))

# Read in the data (not subset)
saveRDS(seurat_data_full, file.path(all_dir, "rda_obj",
                                      "seurat_processed.rds"))
