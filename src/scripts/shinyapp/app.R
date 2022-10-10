library(shiny)
library(tidyverse)
library(cowplot)
library(Seurat)

gene_list <- readRDS("all_genes.rds")

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

source("functions.R")

ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

# Define UI for miles per gallon app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Gene expression"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Selector for variable to plot ----
      selectInput("seurat_object", "Sample to plot",
                  names(seurat_object_list)),
      
      selectInput("meta_data", "Meta data column:",
                  c("SCT_celltype_no_mult",
                    "SCT_celltype_pc",
                    "SCT_celltype_combined",
                    "SCT_cluster",
                    "Annotated")),
      
      selectInput("gene", "Gene of interest:",
                  gene_list),
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Plot of the requested variable against mpg ----
      plotOutput("gene_plot"),
      plotOutput("spatial")
      
    )
  )
)

# Data pre-processing ----------------------------------------------------------

# Define server logic to plot various variables against mpg ----
server <- function(input, output) {
  
  objects <- reactiveValues(seurat_object = NULL)
  
  observeEvent(input$seurat_object, {
    objects$seurat_object <- readRDS(seurat_object_list[[input$seurat_object]])
  })
  
  # Generate plots
  output$gene_plot <- renderPlot({
    #seurat_object <- readRDS(seurat_object_list[[input$seurat_object]])
    
    seurat_object <- objects$seurat_object
    # if(input$seurat_object != "A_955_OvarianTumor"){
    #   seurat_object$Annotated <- "not_annotated"
    # }
    
    ncolors <- length(unique(seurat_object[[input$meta_data]][[1]]))
    colors <- grDevices::colorRampPalette(
      RColorBrewer::brewer.pal(9, "Set1"))(ncolors)
    
    names(colors) <- unique(seurat_object[[input$meta_data]][[1]])
    
    umap1 <- plotDimRed(seurat_object, col_by = input$gene,
                        plot_type = "rnasct.umap")[[1]]
    umap2 <- plotDimRed(seurat_object, col_by = input$meta_data,
                        plot_type = "rnasct.umap",
                        color = colors)[[1]]
    violin1 <- featDistPlot(seurat_object, input$gene,
                            sep_by = input$meta_data,
                            combine = FALSE,
                            color = colors)[[1]]
   
    
    cowplot::plot_grid(umap1, NULL,
                       umap2, violin1,
                       labels = c("A", "", "B", "C"),
                       nrow = 2, ncol = 2,
                       align = "hv",
                       axis = "lr")
  })
  
  output$spatial <- renderPlot({
    seurat_object <- objects$seurat_object
    #seurat_object <- readRDS(seurat_object_list[[input$seurat_object]])
    ncolors <- length(unique(seurat_object[[input$meta_data]][[1]]))
    colors <- grDevices::colorRampPalette(
      RColorBrewer::brewer.pal(9, "Set1"))(ncolors)
    
    names(colors) <- unique(seurat_object[[input$meta_data]][[1]])
    
    plot_gene <- SpatialFeaturePlot(seurat_object, features = input$gene)
    
    plot1 <- SpatialDimPlot(seurat_object, group.by = input$meta_data,
                            cols = colors)
    
    cowplot::plot_grid(plot_gene, plot1,
                       labels = c("D", "E"),
                       nrow = 1, ncol = 2,
                       align = "hv", axis = "lr")
    
  })
  
}

shinyApp(ui, server)