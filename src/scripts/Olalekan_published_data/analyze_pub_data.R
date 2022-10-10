library(Seurat)
library(tidyverse)
library(harmony)

sample = "published_data"

ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

# Set data dirs
base_dir <- "/Users/wellskr/Documents/Analysis/Benjamin_Bitler/visium_ovary/"

source(file.path(base_dir, "src", "scripts", "functions.R"))

data_dir <-
  "/Users/wellskr/Documents/Analysis/Benjamin_Bitler/visium_ovary/resources"

data_files <- list.files(path = data_dir, pattern = "GSM.*csv")

meta_data <- file.path(data_dir, "meta.cells.csv")

meta_df <- read.csv(meta_data)
meta_df$Barcode <- sub("_", "", meta_df$Barcode)

# Pull out idents to match meta df
idents <- sub(".*-([0-9]*)\\.csv", "\\1", data_files)

# Fix the idents that don't match
idents[idents == "3232"] <- "3233"
idents[idents == "3401"] <- "401"
idents[idents == "5150"] <- "5150frozen"

idents <- paste0("omentum", idents)

names(data_files) <- idents

vars.to.regress <- NULL

plot_list <- c("nCount_RNA", "nFeature_RNA",
               "percent.mt")

seurat_assay <- "SCT"

# Make objects -----------------------------------------------------------------
seurat_obj_list <- lapply(names(data_files), function(x){
  data_df <- read.csv(file.path(data_dir, data_files[[x]]),
                      row.names = 1)
  seurat_object <- CreateSeuratObject(counts = data_df,
                                      project = x)
  
  seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object,
                                                        pattern = "^MT\\.")
  
  plot1 <- VlnPlot(seurat_object, features = plot_list, pt.size = 0.1) + 
    NoLegend()
  
  # Remove outliers (percent mt is very high)
  seurat_object <- subset(x = seurat_object, subset = percent.mt < 15 &
                            nFeature_RNA > 500 & nFeature_RNA < 5000)
  
  # Single Cell Transform normalization
  seurat_object <- SCTransform(seurat_object, vars.to.regress = vars.to.regress,
                               verbose = FALSE)
  # Default normalization
  DefaultAssay(seurat_object) <- "RNA"
  seurat_object <- NormalizeData(seurat_object) %>% 
    FindVariableFeatures() %>%
    ScaleData(vars.to.regress = vars.to.regress)
  
})

seurat_merge <- merge(x = seurat_obj_list[[1]],
                      y = seurat_obj_list[2:length(data_files)],
                      add.cell.ids = idents)

# Subset to only idents in the other file
seurat_merge <- subset(seurat_merge, cells = meta_df$Barcode)

# Add in meta data
rownames(meta_df) <- meta_df$Barcode
seurat_merge <- AddMetaData(seurat_merge, meta_df)

seurat_merge <- SCTransform(seurat_merge, vars.to.regress = vars.to.regress,
                             verbose = FALSE)
# PCA --------------------------------------------------------------------------

# PCA of gene expression
seurat_merge <- PCA_dimRed(seurat_merge, assay = seurat_assay)

RNA_plots <- plot_PCA(HTO = FALSE, assay = seurat_assay,
                      sample_object = seurat_merge)

# Need to run harmony to integrate samples based on PCA
seurat_merge <- RunHarmony(seurat_merge, "orig.ident",
                          plot_convergence = TRUE, lambda = 1,
                          theta = 8,
                          reduction = "sctpca")

harmony_plot <- plotDimRed(seurat_merge, col_by = "orig.ident",
                           plot_type = "harmony")
# UMAP -------------------------------------------------------------------------
RNA_pcs <- 30

set.seed(0)
# UMAP of gene expression
umap_data <- group_cells(seurat_merge, sample, data_dir, nPCs = RNA_pcs,
                         resolution = 0.6, assay = seurat_assay, HTO = FALSE,
                         reduction = "harmony")

seurat_merge <- umap_data$object

gene_plots <- umap_data$plots


plotDimRed(seurat_merge, col_by = "assignedCelltype", plot_type = "harmony.umap")

plotDimRed(seurat_merge, col_by = "CD3D", plot_type = "harmony.umap")


saveRDS(seurat_merge, file = file.path(data_dir, "seurat_processed.rds"))

