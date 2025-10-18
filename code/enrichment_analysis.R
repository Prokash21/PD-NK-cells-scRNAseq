# Gene Set Enrichment Analysis
library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(msigdbr)
library(tidyverse)

# Load DE results
de_files <- list.files(pattern = "DE_.*\\.csv$")

# Function to perform GSEA
perform_gsea <- function(de_results, title) {
    # Create ranked gene list
    ranked_genes <- de_results$avg_log2FC
    names(ranked_genes) <- rownames(de_results)
    ranked_genes <- sort(ranked_genes, decreasing = TRUE)
    
    # Get MSigDB gene sets
    hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
    
    # Prepare gene sets
    pathways <- split(hallmark_sets$gene_symbol, hallmark_sets$gs_name)
    
    # Run GSEA
    gsea_results <- GSEA(ranked_genes, 
                        TERM2GENE = pathways,
                        pvalueCutoff = 0.05)
    
    # Create plots
    p1 <- gseaplot2(gsea_results, 
                    geneSetID = 1:5, 
                    pvalue_table = TRUE,
                    title = title)
    
    # Save results
    write.csv(gsea_results@result, 
              file = paste0("GSEA_", title, ".csv"))
    ggsave(paste0("GSEA_plot_", title, ".pdf"), 
           plot = p1, 
           width = 12, 
           height = 8)
    
    return(gsea_results)
}

# Run GSEA for each DE comparison
for(file in de_files) {
    de_results <- read.csv(file, row.names = 1)
    title <- gsub("DE_|.csv", "", file)
    gsea_res <- perform_gsea(de_results, title)
}

# Additional GO enrichment analysis
perform_go_analysis <- function(de_results, title) {
    # Get significant genes
    sig_genes <- rownames(de_results[de_results$p_val_adj < 0.05 & 
                                   abs(de_results$avg_log2FC) > 1, ])
    
    # Run GO enrichment
    go_results <- enrichGO(gene = sig_genes,
                          OrgDb = org.Hs.eg.db,
                          keyType = "SYMBOL",
                          ont = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05)
    
    # Create plots
    p1 <- dotplot(go_results, showCategory = 20) +
        ggtitle(paste("GO enrichment -", title))
    
    # Save results
    write.csv(go_results@result,
              file = paste0("GO_", title, ".csv"))
    ggsave(paste0("GO_plot_", title, ".pdf"),
           plot = p1,
           width = 12,
           height = 8)
    
    return(go_results)
}

# Run GO analysis for each DE comparison
for(file in de_files) {
    de_results <- read.csv(file, row.names = 1)
    title <- gsub("DE_|.csv", "", file)
    go_res <- perform_go_analysis(de_results, title)
}