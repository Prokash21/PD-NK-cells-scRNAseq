# GSE141578 analysis

Instructions to reproduce data processing and analysis for GSE141578.

Setup (macOS / Linux):

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Run the analysis script:

```bash
python scripts/analyze_gse141578.py
```

Outputs:
- data/processed/gse141578_adata.h5ad
- figures/umap_sample_leiden.png
- results/ranked_genes_groups.csv
