required_pkgs <- c("Seurat", "SeuratObject", "ggplot2", "patchwork")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(patchwork)
})

project_dir <- normalizePath(getwd())
data_dir <- file.path(project_dir, "..", "scRNA-seq_scanpy", "GSE141578_RAW_split")
results_dir <- file.path(project_dir, "results")
fig_dir <- file.path(project_dir, "figures")

if (!dir.exists(data_dir)) {
  stop("Raw 10x data not found at: ", data_dir, "\nCheck the folder path or update the data_dir variable.")
}

if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

gsm_to_condition <- c(
  "GSM4208766" = "PD",
  "GSM4208768" = "PD",
  "GSM4208769" = "PD",
  "GSM4208770" = "PD",
  "GSM4208771" = "PD",
  "GSM4208772" = "HC",
  "GSM4208773" = "HC",
  "GSM4208774" = "HC",
  "GSM4208775" = "HC",
  "GSM4208776" = "HC",
  "GSM4208777" = "HC",
  "GSM4208778" = "HC",
  "GSM4208779" = "HC",
  "GSM4208780" = "HC",
  "GSM4404053" = "PD",
  "GSM4404054" = "PD",
  "GSM4404055" = "HC",
  "GSM4404056" = "PD",
  "GSM4404057" = "PD",
  "GSM5579354" = "HC",
  "GSM5579355" = "HC",
  "GSM5579356" = "PD",
  "GSM5579357" = "HC",
  "GSM5579358" = "PD",
  "GSM5579359" = "PD",
  "GSM5579360" = "PD",
  "GSM5579361" = "PD",
  "GSM5579362" = "PD",
  "GSM5579364" = "HC",
  "GSM5579365" = "HC"
)

all_dirs <- list.dirs(data_dir, full.names = TRUE, recursive = FALSE)
sample_dirs <- sort(Filter(function(path) {
  name <- basename(path)
  if (!grepl("^GSM", name)) return(FALSE)
  gsm_id <- sub("_.*", "", name)
  gsm_id %in% names(gsm_to_condition)
}, all_dirs))

if (length(sample_dirs) == 0) {
  stop("No valid PD/HC 10x sample folders were found in: ", data_dir)
}

message("Found ", length(sample_dirs), " sample(s) to process")

sample_names <- vapply(sample_dirs, function(path) basename(path), character(1))

samples <- lapply(seq_along(sample_dirs), function(i) {
  sample_path <- sample_dirs[i]
  sample_name <- basename(sample_path)
  gsm_id <- sub("_.*", "", sample_name)
  condition <- gsm_to_condition[[gsm_id]]
  sample_label <- paste0(gsm_id, "_", condition)

  counts <- Read10X(data.dir = sample_path)
  obj <- CreateSeuratObject(
    counts = counts,
    project = sample_label,
    min.cells = 3,
    min.features = 200
  )

  obj$sample <- sample_label
  obj$gsm <- gsm_id
  obj$condition <- condition
  obj$percent.mt <- PercentageFeatureSet(obj, pattern = "^MT-")

  keep <- which(obj$nFeature_RNA > 200 & obj$nFeature_RNA < 6000 & obj$percent.mt < 20)

  if (length(keep) == 0) {
    message("Skipping empty sample after QC: ", sample_label)
    return(NULL)
  }

  obj[, keep]
})

samples <- Filter(Negate(is.null), samples)
if (length(samples) == 0) {
  stop("No samples remained after QC filtering.")
}

sample_labels <- vapply(samples, function(obj) unique(obj$sample), character(1))

# merge all non-empty samples
combined <- samples[[1]]
if (length(samples) > 1) {
  for (i in 2:length(samples)) {
    combined <- merge(
      combined,
      y = samples[[i]],
      add.cell.ids = c(sample_labels[i - 1], sample_labels[i])
    )
  }
}

combined$sample <- factor(combined$sample)
combined$condition <- factor(combined$condition)

combined <- NormalizeData(combined, normalization.method = "LogNormalize", scale.factor = 10000)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
combined <- ScaleData(combined, vars.to.regress = "percent.mt")
combined <- RunPCA(combined, npcs = 50, verbose = FALSE)
combined <- FindNeighbors(combined, dims = 1:30)
combined <- FindClusters(combined, resolution = 0.5)
combined <- RunUMAP(combined, dims = 1:30)

# Save the base analyzed object
saveRDS(combined, file.path(results_dir, "seurat_pd_gse141578.rds"))

p1 <- DimPlot(combined, reduction = "umap", group.by = "sample", pt.size = 0.5) +
  ggtitle("UMAP by sample")

p2 <- DimPlot(combined, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle("UMAP clusters")

p3 <- VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)

ggsave(file.path(fig_dir, "umap_clusters.png"), p2, width = 8, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "umap_by_sample.png"), p1, width = 8, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "qc_violin.png"), p3, width = 10, height = 6, dpi = 300)

message("Finished preprocessing and clustering.")
message("Saved object to: ", file.path(results_dir, "seurat_pd_gse141578.rds"))
message("Saved figures to: ", fig_dir)

# Return the object for downstream scripts
combined
