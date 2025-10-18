# Cell Type Annotation
library(Seurat)
library(SingleR)
library(celldex)
library(tidyverse)

# Load the Seurat object
load("seurat_obj.RData")

# Get reference dataset (Human Primary Cell Atlas)
ref <- HumanPrimaryCellAtlasData()

# Extract expression matrix
expr_matrix <- GetAssayData(nsclc.seurat.obj, slot = "data")

# Run SingleR
predictions <- SingleR(test = expr_matrix,
                      ref = ref,
                      labels = ref$label.main)

# Add cell type annotations to Seurat object
nsclc.seurat.obj$cell_type <- predictions$labels

# Visualization
DimPlot(nsclc.seurat.obj, 
        group.by = "cell_type", 
        label = TRUE, 
        label.size = 3) + 
    theme(legend.position = "right")

# Create composition plot
cell_counts <- table(Idents(nsclc.seurat.obj), nsclc.seurat.obj$cell_type)
pheatmap::pheatmap(cell_counts)

# Save annotated object
saveRDS(nsclc.seurat.obj, "annotated_seurat_obj.rds")