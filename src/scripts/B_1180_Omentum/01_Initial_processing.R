library(tidyverse)
library(Seurat)

ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

# Variables
sample <- "B_1180_Omentum"
ffpe <- TRUE # Setting to TRUE will assume a probe set was used.
vars.to.regress <- NULL

# Set data dirs
base_dir <- "/Users/wellskr/Documents/Analysis/Benjamin_Bitler/visium_ovary/"

source(file.path(base_dir, "src", "scripts", "functions.R"))

results_dir <- file.path(base_dir, "results")

save_dir <- file.path(results_dir, sample, "R_analysis")

# Make output directories
ifelse(!dir.exists(file.path(save_dir)),
       dir.create(file.path(save_dir)), FALSE)

ifelse(!dir.exists(file.path(save_dir, "images")),
       dir.create(file.path(save_dir, "images")), FALSE)

ifelse(!dir.exists(file.path(save_dir, "files")),
       dir.create(file.path(save_dir, "files")), FALSE)

ifelse(!dir.exists(file.path(save_dir, "rda_obj")),
       dir.create(file.path(save_dir, "rda_obj")), FALSE)

seurat_obj <- create_spatial_seurat(results_dir = results_dir,
                                     sample_name = sample)

if (ffpe){
  plot_list <- c("nCount_Spatial", "nFeature_Spatial")
} else {
  # Add mitochondrial percent
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj,
                                                     pattern = "^MT-")
  plot_list <- c("nCount_Spatial", "nFeature_Spatial", "percent.mt")
}


plot1 <- VlnPlot(seurat_obj, features = plot_list, pt.size = 0.1) + 
  NoLegend()
plot2 <- lapply(plot_list, function(x) {
  SpatialFeaturePlot(seurat_obj, features = x) + 
  theme(legend.position = "right")
})

pdf(file.path(save_dir, "images", "initial_processing.pdf"))

print(plot1)
print(plot2)

dev.off()

# Add cell cycle scores
seurat_obj <- CellCycleScoring(seurat_obj, s.features = cc.genes$s.genes,
                               g2m.features = cc.genes$g2m.genes,
                               set.ident = FALSE)


seurat_obj <- SCTransform(seurat_obj, assay = "Spatial", verbose = FALSE,
                          vars.to.regress = vars.to.regress)

DefaultAssay(seurat_obj) <- "Spatial"

seurat_obj <- NormalizeData(seurat_obj) %>% 
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = vars.to.regress)

saveRDS(seurat_obj, file = file.path(save_dir, "rda_obj", "seurat_start.rds"))
