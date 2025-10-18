
<div align="center">
  <h1>scRNA-seq Analysis Pipeline</h1>
  <img src="https://img.shields.io/badge/R-Seurat-blue?logo=R" alt="Seurat">
  <img src="https://img.shields.io/badge/License-MIT-green" alt="MIT License">
</div>

---

<details>
<summary><strong>Table of Contents</strong></summary>

- [Overview](#overview)
- [Folder Structure](#folder-structure)
- [Main Scripts](#main-scripts)
- [Usage](#usage)
- [Output](#output)
- [License](#license)
- [Author](#author)
</details>

---

## 🧬 Overview

This repository contains a comprehensive pipeline for single-cell RNA sequencing (**scRNA-seq**) analysis using **R** and **Seurat**. The workflow covers:

- Quality control
- Clustering
- Cell type annotation
- Trajectory analysis
- Differential expression
- Enrichment analysis

---

## 📁 Folder Structure

| Folder      | Description                       |
| ----------- | --------------------------------- |
| `code/`     | R scripts for each analysis step   |
| `figure/`   | Output figures and plots           |

---

## 📝 Main Scripts

| Script                    | Purpose                                   |
|---------------------------|-------------------------------------------|
| `scRNA_run.R`             | Main pipeline for preprocessing & clustering |
| `marker_analysis.R`       | Marker gene identification & visualization |
| `cell_annotation.R`       | Automated cell type annotation             |
| `trajectory_analysis.R`   | Trajectory & pseudotime analysis           |
| `differential_expression.R` | Differential expression analysis         |
| `enrichment_analysis.R`   | Gene set enrichment & GO analysis          |

---

## 🚀 Usage


1. **Install required R packages:**
  - All required packages are listed in [`requirements.txt`](./requirements.txt).
  - You can install them manually or use the following R commands:
    ```r
    pkgs <- readLines("requirements.txt")
    install.packages(pkgs[!pkgs %in% c("SingleR", "celldex", "monocle3", "clusterProfiler", "org.Hs.eg.db", "enrichplot", "msigdbr", "BiocManager")])
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(pkgs[pkgs %in% c("SingleR", "celldex", "monocle3", "clusterProfiler", "org.Hs.eg.db", "enrichplot", "msigdbr")])
    ```
2. **Run scripts in order** as described above.
3. **Output files** (`.csv`, `.pdf`, plots) are saved in the respective folders.

---

## 📊 Output

- Cluster markers
- Cell type annotations
- Trajectory plots
- Differential expression results
- Enrichment analysis tables and figures

---

## 📄 License

This project is licensed under the MIT License.

---

## 👤 Author

**Prokash21**
