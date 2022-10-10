library(Seurat)
library(here)

seurat_object_list <- c("A_955_OvarianTumor" = 
                          "A_955_OvarianTumor.rds",
                        "A_GTFB1154_OvarianCancerTumor" =
                          "A_GTFB1154_OvarianCancerTumor.rds",
                        "B_1180_Omentum" =
                          "B_1180_Omentum.rds",
                        "B_GTFB1191_OvarianCancerTumor" =
                          "B_GTFB1191_OvarianCancerTumor.rds",
                        "C_GTFB1170_SmallCellOvarianCancer" =
                          "C_GTFB1170_SmallCellOvarianCancer.rds",
                        "C_GTFB_1191_OmentumTumor" =
                          "C_GTFB_1191_OmentumTumor.rds",
                        "D_GTFB1170_SmallCellOvarianCancer" =
                          "D_GTFB1170_SmallCellOvarianCancer.rds",
                        "D_GTFB_1230_OvarianTumor" =
                          "D_GTFB_1230_OvarianTumor.rds")

# Make smaller files for shiny --> remove SCT and unneeded metadata
keep_clust <- c("SCT_celltype_no_mult",
                "SCT_celltype_pc",
                "SCT_celltype_combined",
                "SCT_cluster",
                "Annotated")

`%notin%` <- Negate(`%in%`)

invisible(lapply(names(seurat_object_list), function(x){
  print(x)
  seurat_object <- readRDS(here("results", x, "R_analysis_no_multi", "rda_obj",
                                "seurat_no_mult_processed.rds"))
  
  if(x != "A_955_OvarianTumor"){
    seurat_object$Annotated <- "Not_annotated"
  }
  for(i in colnames(seurat_object[[]])){
    if(i %notin% keep_clust){
      seurat_object[[i]] <- NULL
    }
  }
  
  print(head(seurat_object))
  seurat_object[["SCT"]] <- NULL
  
  saveRDS(seurat_object, here("src/scripts/shinyapp", seurat_object_list[[x]]))
}))

# Get list of genes in all objects
all_genes <- lapply(seurat_object_list, function(x){
  seurat_obj <- readRDS(file.path("src/scripts/shinyapp", x))
  return(rownames(seurat_obj))
})

all_genes <- unique(unlist(all_genes))

saveRDS(all_genes, "src/scripts/shinyapp/all_genes.rds")

