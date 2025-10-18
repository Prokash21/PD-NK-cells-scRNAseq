# Trajectory Analysis
library(Seurat)
library(monocle3)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(viridis)

# Create output directory for plots
dir.create("trajectory_plots", showWarnings = FALSE)

# Load the annotated Seurat object
seurat_obj <- readRDS("annotated_seurat_obj.rds")

# Convert Seurat to CDS
DefaultAssay(seurat_obj) <- "RNA"
cds <- new_cell_data_set(GetAssayData(seurat_obj, slot = "counts"),
                        cell_metadata = seurat_obj@meta.data,
                        gene_metadata = data.frame(gene_short_name = rownames(seurat_obj),
                                                 row.names = rownames(seurat_obj)))

# Pre-process the data
cds <- preprocess_cds(cds, num_dim = 50)
cds <- align_cds(cds)

# Reduce dimensionality and cluster
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)

# Learn the trajectory
cds <- learn_graph(cds)

# Order cells
cds <- order_cells(cds)

# Generate and save all available plots

# 1. Basic trajectory plot with pseudotime
p1 <- plot_cells(cds,
                color_cells_by = "pseudotime",
                label_groups_by_cluster = FALSE,
                label_leaves = FALSE,
                label_branch_points = FALSE) +
      ggtitle("Trajectory colored by pseudotime")
ggsave("trajectory_plots/1_basic_trajectory_pseudotime.pdf", p1, width = 10, height = 8)

# 2. Trajectory with cluster labels
p2 <- plot_cells(cds,
                color_cells_by = "cluster",
                label_groups_by_cluster = TRUE,
                label_leaves = TRUE,
                label_branch_points = TRUE) +
      ggtitle("Trajectory with cluster labels")
ggsave("trajectory_plots/2_trajectory_clusters.pdf", p2, width = 10, height = 8)

# 3. Trajectory faceted by cell type
p3 <- plot_cells(cds,
                color_cells_by = "cell_type",
                label_cell_groups = FALSE,
                label_leaves = FALSE,
                label_branch_points = FALSE) +
      facet_wrap(~cell_type, nrow = 2) +
      ggtitle("Trajectory faceted by cell type")
ggsave("trajectory_plots/3_trajectory_by_celltype.pdf", p3, width = 15, height = 10)

# 4. Density plot of cells along pseudotime
cell_pseudotime <- pseudotime(cds)
pseudotime_df <- data.frame(pseudotime = cell_pseudotime,
                           cell_type = colData(cds)$cell_type)
p4 <- ggplot(pseudotime_df, aes(x = pseudotime, fill = cell_type)) +
      geom_density(alpha = 0.5) +
      theme_bw() +
      ggtitle("Cell density along pseudotime") +
      scale_fill_viridis_d()
ggsave("trajectory_plots/4_pseudotime_density.pdf", p4, width = 10, height = 6)

# 5. Gene expression dynamics along trajectory
# Select top variable genes
top_genes <- top_marker_feature_genes(cds, n_features = 6)
p5 <- plot_genes_in_pseudotime(cds[top_genes,],
                              color_cells_by = "cell_type",
                              min_expr = 0.5) +
      ggtitle("Gene expression dynamics along trajectory")
ggsave("trajectory_plots/5_gene_expression_dynamics.pdf", p5, width = 12, height = 8)

# 6. Branch point analysis
# Get branch points
branch_points <- get_earliest_principal_node(cds)
branch_cells <- which(principle_graph_test(cds, branch_points) < 0.05)
branch_genes <- row.names(subset(fData(cds),
                               gene_short_name %in% top_genes))

# Plot branch analysis
p6 <- plot_cells(cds,
                genes = branch_genes,
                show_trajectory_graph = TRUE,
                label_branch_points = TRUE,
                label_cell_groups = FALSE,
                label_leaves = TRUE) +
      ggtitle("Branch point analysis")
ggsave("trajectory_plots/6_branch_analysis.pdf", p6, width = 12, height = 8)

# 7. Combined visualization using patchwork
combined_plot <- (p1 + p2) / (p3 + p4)
ggsave("trajectory_plots/7_combined_trajectory_analysis.pdf", combined_plot, width = 20, height = 16)

# 8. 3D trajectory plot if available
if(requireNamespace("plotly", quietly = TRUE)) {
    library(plotly)
    reducedDims(cds)$UMAP_3D <- reducedDims(cds)$UMAP[,1:3]
    p8 <- plot_cells_3d(cds,
                        reduction_method = "UMAP_3D",
                        color_cells_by = "pseudotime")
    htmlwidgets::saveWidget(p8, "trajectory_plots/8_3D_trajectory.html")
}

# Save the trajectory object
saveRDS(cds, "trajectory.rds")

# Create a summary plot with key statistics
metadata <- as.data.frame(colData(cds))
p9 <- ggplot(metadata, aes(x = pseudotime, y = cell_type)) +
      geom_violin(aes(fill = cell_type), scale = "width") +
      theme_bw() +
      ggtitle("Pseudotime distribution by cell type") +
      scale_fill_viridis_d()
ggsave("trajectory_plots/9_pseudotime_summary.pdf", p9, width = 12, height = 8)

# 10. Create a heatmap of gene expression changes along pseudotime
gene_module_df <- find_gene_modules(cds[top_genes,])
p10 <- plot_gene_modules(cds[top_genes,], gene_module_df) +
       ggtitle("Gene modules along pseudotime")
ggsave("trajectory_plots/10_gene_modules_heatmap.pdf", p10, width = 12, height = 8)

# 11. Create trajectory backbone with density
p11 <- plot_cells(cds,
                 color_cells_by = "pseudotime",
                 show_trajectory_graph = TRUE,
                 graph_label_size = 4,
                 alpha = 0.5) +
       stat_density2d(aes(color = ..level..)) +
       scale_color_viridis_c() +
       ggtitle("Trajectory backbone with density contours")
ggsave("trajectory_plots/11_trajectory_density.pdf", p11, width = 10, height = 8)

# 12. Create state transitions plot
cell_states <- choose_cells(cds)
p12 <- plot_cells(cds,
                 color_cells_by = "state",
                 show_trajectory_graph = TRUE,
                 label_groups_by_cluster = TRUE,
                 group_cells_by = "state") +
       ggtitle("Cell state transitions")
ggsave("trajectory_plots/12_state_transitions.pdf", p12, width = 10, height = 8)

# 13. Generate pseudotime heatmap with cell annotations
pt_matrix <- as.matrix(exprs(cds)[top_genes,order(pseudotime(cds))])
cell_annotations <- data.frame(
    Pseudotime = pseudotime(cds)[order(pseudotime(cds))],
    Cell_Type = colData(cds)$cell_type[order(pseudotime(cds))]
)
p13 <- pheatmap::pheatmap(pt_matrix,
                         annotation_col = cell_annotations,
                         show_colnames = FALSE,
                         cluster_cols = FALSE,
                         main = "Gene expression dynamics with annotations")
ggsave("trajectory_plots/13_pseudotime_heatmap.pdf", p13, width = 12, height = 10)

# 14. Create velocity estimates if spliced/unspliced counts available
if("spliced" %in% names(assays(cds))) {
    p14 <- plot_cells(cds,
                     color_cells_by = "spliced_ratio",
                     label_cell_groups = FALSE) +
           ggtitle("RNA velocity estimation")
    ggsave("trajectory_plots/14_velocity_estimation.pdf", p14, width = 10, height = 8)
}

# 15. Generate branch point decision visualization
branch_manifold <- choose_graph_segments(cds)
p15 <- plot_cells(cds,
                 color_cells_by = "branch_manifold",
                 show_trajectory_graph = TRUE,
                 label_groups_by_cluster = FALSE) +
       ggtitle("Branch point decisions")
ggsave("trajectory_plots/15_branch_decisions.pdf", p15, width = 10, height = 8)

# 16. Create time series of cell type proportions
pseudotime_bins <- cut(pseudotime(cds), breaks = 10)
cell_type_proportions <- table(pseudotime_bins, colData(cds)$cell_type)
prop_df <- as.data.frame(prop.table(cell_type_proportions, 1))
colnames(prop_df) <- c("Pseudotime", "Cell_Type", "Proportion")

p16 <- ggplot(prop_df, aes(x = Pseudotime, y = Proportion, fill = Cell_Type)) +
       geom_bar(stat = "identity") +
       theme_bw() +
       theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
       scale_fill_viridis_d() +
       ggtitle("Cell type proportions along pseudotime")
ggsave("trajectory_plots/16_celltype_proportions.pdf", p16, width = 12, height = 8)

# 17. Create network visualization of gene modules
if(requireNamespace("igraph", quietly = TRUE)) {
    library(igraph)
    gene_module_network <- graph_from_adjacency_matrix(
        cor(t(exprs(cds)[top_genes,])) > 0.7,
        mode = "undirected", diag = FALSE
    )
    
    p17 <- plot(gene_module_network,
                vertex.size = 5,
                vertex.label.cex = 0.7,
                main = "Gene module correlation network")
    ggsave("trajectory_plots/17_gene_network.pdf", p17, width = 10, height = 10)
}

# 18. Create bifurcation analysis plot
root_cells <- get_earliest_principal_node(cds)
bifurcation_df <- data.frame(
    pseudotime = pseudotime(cds),
    cell_type = colData(cds)$cell_type,
    branch = factor(principal_graph_test(cds, root_cells))
)

p18 <- ggplot(bifurcation_df, aes(x = pseudotime, y = cell_type, color = branch)) +
       geom_jitter(alpha = 0.5) +
       scale_color_viridis_d() +
       theme_bw() +
       ggtitle("Bifurcation analysis")
ggsave("trajectory_plots/18_bifurcation_analysis.pdf", p18, width = 12, height = 8)

# Print updated analysis summary
cat("Trajectory analysis completed. Generated plots:\n",
    "1. Basic trajectory with pseudotime\n",
    "2. Trajectory with cluster labels\n",
    "3. Trajectory by cell type\n",
    "4. Pseudotime density distribution\n",
    "5. Gene expression dynamics\n",
    "6. Branch point analysis\n",
    "7. Combined visualization\n",
    "8. 3D trajectory (if plotly is installed)\n",
    "9. Pseudotime distribution summary\n",
    "10. Gene modules heatmap\n",
    "11. Trajectory backbone with density\n",
    "12. State transitions\n",
    "13. Pseudotime heatmap with annotations\n",
    "14. RNA velocity estimation (if spliced data available)\n",
    "15. Branch point decisions\n",
    "16. Cell type proportions over pseudotime\n",
    "17. Gene module correlation network\n",
    "18. Bifurcation analysis\n",
    "\nAll plots saved in 'trajectory_plots' directory\n")