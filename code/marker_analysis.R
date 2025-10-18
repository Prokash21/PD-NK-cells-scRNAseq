# Marker Gene Analysis
library(Seurat)
library(tidyverse)
library(ggplot2)

# Load the Seurat object
load("seurat_obj.RData")

# Find markers for each cluster
all.markers <- FindAllMarkers(nsclc.seurat.obj,
                            only.pos = TRUE,
                            min.pct = 0.25,
                            logfc.threshold = 0.25)

# Get top markers per cluster
top10 <- all.markers %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC)

# Create heatmap of top markers
DoHeatmap(nsclc.seurat.obj, features = top10$gene) + NoLegend()

# Generate violin plots for top markers
VlnPlot(nsclc.seurat.obj, features = unique(top10$gene[1:5]), ncol = 5)

# Feature plots for spatial visualization
FeaturePlot(nsclc.seurat.obj, 
            features = unique(top10$gene[1:6]), 
            ncol = 3)

# Save results
write.csv(all.markers, "cluster_markers.csv")
saveRDS(top10, "top10_markers.rds")