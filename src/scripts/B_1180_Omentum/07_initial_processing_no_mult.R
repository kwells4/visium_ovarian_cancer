library(ggplot2)
library(Seurat)
library(tidyr)
library(cowplot)
library(dplyr)
library(clustifyr)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

sample <- "B_1180_Omentum"

normalization_method <- "log" # can be SCT or log

vars.to.regress <- NULL

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

save_dir_full <- file.path(base_dir_proj, "R_analysis")

save_dir <- file.path(base_dir_proj, "R_analysis_no_multi")

# Create new directories
ifelse(!dir.exists(save_dir),
       dir.create(save_dir),
       FALSE)

ifelse(!dir.exists(file.path(save_dir, "rda_obj")),
       dir.create(file.path(save_dir, "rda_obj")),
       FALSE)

ifelse(!dir.exists(file.path(save_dir, "images")),
       dir.create(file.path(save_dir, "images")),
       FALSE)

ifelse(!dir.exists(file.path(save_dir, "files")),
       dir.create(file.path(save_dir, "files")),
       FALSE)

# Information for cell mapping
data_dir <-
  "/Users/wellskr/Documents/Analysis/Benjamin_Bitler/visium_ovary/resources"

ref_mat <- read.csv(file.path(data_dir,
                              "average_cluster_expression_my_clusters.csv"),
                    header = TRUE, row.names = 1)

# Read in the data
seurat_data <-  readRDS(file.path(save_dir_full,
                                  "rda_obj", "seurat_processed.rds"))

Idents(seurat_data) <- "high_res_cluster"

all_idents <- unique(Idents(seurat_data))
all_idents <- all_idents[!(all_idents == "multiple")]

seurat_sub <- subset(seurat_data, ident = all_idents)

seurat_obj <- SCTransform(seurat_sub, assay = "Spatial", verbose = FALSE,
                          vars.to.regress = vars.to.regress)

DefaultAssay(seurat_sub) <- "Spatial"

seurat_sub <- FindVariableFeatures(seurat_sub) %>%
  ScaleData(vars.to.regress = vars.to.regress)

saveRDS(seurat_obj, file = file.path(save_dir, "rda_obj",
                                     "seurat_no_mult_start.rds"))