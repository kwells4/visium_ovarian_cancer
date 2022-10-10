
image <- Images(object = seurat_data,
                assay = DefaultAssay(object = seurat_data))

coords <- GetTissueCoordinates(object = seurat_data[[image]])

# Get coords
coords <- seurat_data[["slice1"]]@coordinates

image <- seurat_data[["slice1"]]@image

spot_radius <- seurat_data[["slice1"]]@spot.radius