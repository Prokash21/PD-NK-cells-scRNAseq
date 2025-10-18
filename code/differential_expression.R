# Differential Expression Analysis
library(Seurat)
library(tidyverse)
library(EnhancedVolcano)

# Load annotated Seurat object
seurat_obj <- readRDS("annotated_seurat_obj.rds")

# Function to perform DE analysis between groups
perform_de_analysis <- function(seurat_obj, group1, group2, ident.use = "cell_type") {
    Idents(seurat_obj) <- ident.use
    de_results <- FindMarkers(seurat_obj,
                            ident.1 = group1,
                            ident.2 = group2,
                            min.pct = 0.25,
                            logfc.threshold = 0.25)
    return(de_results)
}

# Get unique cell types
cell_types <- unique(seurat_obj$cell_type)

# Perform pairwise comparisons
for(i in 1:(length(cell_types)-1)) {
    for(j in (i+1):length(cell_types)) {
        de_res <- perform_de_analysis(seurat_obj, cell_types[i], cell_types[j])
        
        # Create volcano plot
        p <- EnhancedVolcano(de_res,
                            lab = rownames(de_res),
                            x = 'avg_log2FC',
                            y = 'p_val_adj',
                            title = paste(cell_types[i], 'vs', cell_types[j]))
        
        # Save results
        write.csv(de_res, 
                 file = paste0('DE_', 
                             gsub(' ', '_', cell_types[i]), '_vs_',
                             gsub(' ', '_', cell_types[j]), '.csv'))
        
        ggsave(paste0('volcano_', 
                     gsub(' ', '_', cell_types[i]), '_vs_',
                     gsub(' ', '_', cell_types[j]), '.pdf'), 
               plot = p)
    }
}