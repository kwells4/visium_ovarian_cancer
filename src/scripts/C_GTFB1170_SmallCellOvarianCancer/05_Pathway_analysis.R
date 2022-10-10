library(ggplot2)
library(Seurat)
library(tidyr)
library(cowplot)
library(dplyr)
library(pathview)
library(openxlsx)
library(gprofiler2)

pval <- 0.05
logfc <- 0.5

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

sample <- "C_GTFB1170_SmallCellOvarianCancer"
cell_types <- "SCT_celltype"
clusters <- "SCT_cluster"

DE_assay <- "Spatial" # Keep spatial or RNA for now

# Set directories
base_dir <- "/Users/wellskr/Documents/Analysis/Benjamin_Bitler/visium_ovary"

gen_id <- "hsa"

path_id_list <- c(T_cell_receptor = "04660",
                  pathways_in_cancer = "05200",
                  leukocyte_transendothelial_migration = "04670",
                  transcriptional_misregulation_in_cancer = "05202",
                  p53_signaling_pathway = "04115",
                  proteoglycans_in_cancer = "05205")

normalization_method <- "SCT" # can be SCT or log


if(normalization_method == "SCT"){
  SCT <- TRUE
  seurat_assay <- "SCT"
} else {
  SCT <- FALSE
  seurat_assay <- "RNA"
}

source(file.path(base_dir, "src", "scripts", "functions.R"))

base_dir_proj <- file.path(base_dir, "results", sample)

save_dir <- file.path(base_dir_proj, "R_analysis")

out_dir <- file.path(save_dir, "images", "pathways")

out_dir_go <- file.path(save_dir, "images", "GSEA")
out_dir_text <- file.path(save_dir, "files", "GSEA")

# Create output directory
out_dir %>%
  dir.create(showWarnings = F)

out_dir_go %>%
  dir.create(showWarnings = F)

out_dir_text %>%
  dir.create(showWarnings = F)


# Read in the data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))


# Cell type DE -----------------------------------------------------------------

# Find markers
marker_genes <- find_markers(seurat_object = seurat_data,
                             seurat_assay = DE_assay,
                             test_idents = cell_types)


# Run gost
gost_output <- run_gost(seurat_de = marker_genes,
                        sources = c("GO:BP", "KEGG", "REAC", "TF"))

# Save gost output
save_gost(gost_output,
          save_dir_plots = file.path(out_dir_go, "celltype"),
          save_dir_text = file.path(out_dir_text, "celltype"))



de_to_pathview(seurat_de = marker_genes,
               path_id_list = path_id_list,
               out_dir = file.path(out_dir, "celltype"),
               seurat_assay = DE_assay,
               test_idents = cell_types,
               gen_id = gen_id)

# cluster DE -------------------------------------------------------------------

# Find markers
marker_genes <- find_markers(seurat_object = seurat_data,
                             seurat_assay = DE_assay,
                             test_idents = clusters)


# Run gost
gost_output <- run_gost(seurat_de = marker_genes,
                        sources = c("GO:BP", "KEGG", "REAC", "TF"))

# Save gost output
save_gost(gost_output,
          save_dir_plots = file.path(out_dir_go, "cluster"),
          save_dir_text = file.path(out_dir_text, "cluster"))


de_to_pathview(seurat_de = marker_genes,
               path_id_list = path_id_list,
               out_dir = file.path(out_dir, "cluster"),
               seurat_assay = DE_assay,
               test_idents = clusters,
               gen_id = gen_id)

