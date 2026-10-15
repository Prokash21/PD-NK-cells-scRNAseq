import warnings; warnings.filterwarnings("ignore")
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]   # PD_CSF_NK/
DATA = ROOT / "data"                         # h5ad files from 01_qc_integration_annotation.ipynb
OUT = ROOT / "results"; OUT.mkdir(exist_ok=True)
pd.set_option("display.width", 200, "display.max_columns", 20)

adata = sc.read_h5ad(DATA / "integrated.h5ad")
full = adata.raw.to_adata()            # log-normalized, all 15,553 genes
full.obs = adata.obs.copy()
full.obsm["X_umap"] = adata.obsm["X_umap"]

# ---------- 1. NK markers vs all other cells ----------
sc.tl.rank_genes_groups(full, "cell type", groups=["NK cell"], reference="rest",
                        method="wilcoxon", pts=True)
nk = sc.get.rank_genes_groups_df(full, "NK cell")
nk = nk.rename(columns={"pct_nz_group": "pct_NK", "pct_nz_reference": "pct_rest"})
sig = nk[(nk.pvals_adj < 0.05) & (nk.logfoldchanges > 1) & (nk.pct_NK > 0.25)]
sig = sig.sort_values("scores", ascending=False)

# NK vs CD8+ T (closest neighbour; removes shared cytotoxic genes)
sc.tl.rank_genes_groups(full, "cell type", groups=["NK cell"], reference="CD8+ T-cell",
                        method="wilcoxon", pts=True, key_added="nk_vs_cd8")
vcd8 = sc.get.rank_genes_groups_df(full, "NK cell", key="nk_vs_cd8")
vcd8 = vcd8[(vcd8.pvals_adj < 0.05) & (vcd8.logfoldchanges > 1) & (vcd8.pct_nz_group > 0.25)]
specific = sig[sig.names.isin(vcd8.names)]

# scVI DE for cluster 6 already stored in the object
scvi = adata.uns["scvi_markers"]
scvi6 = scvi[scvi.group1.astype(str) == "6"].sort_values("lfc_mean", ascending=False)

sig.to_csv(OUT / "NK_markers_vs_all.csv", index=False)
specific.to_csv(OUT / "NK_markers_vs_all_and_CD8T.csv", index=False)
scvi6.to_csv(OUT / "NK_markers_scVI_cluster6.csv")

print("NK cells:", (full.obs["cell type"] == "NK cell").sum())
print("\n== NK vs all: n sig =", len(sig)); print(sig.head(40).round(3).to_string(index=False))
print("\n== also up vs CD8+ T: n =", len(specific)); print(specific.head(30).names.tolist())
print("\n== scVI cluster 6 top 30:", scvi6.index[:30].tolist(), " n =", len(scvi6))

# ---------- 2. PD vs HC within NK ----------
pd_ids = {"GSM4208766","GSM4208768","GSM4208769","GSM4208770","GSM4208771","GSM4404053"}
hc_ids = {f"GSM42087{i}" for i in range(72, 81)}
gsm = full.obs.Sample.str.split("_").str[0]
full.obs["condition"] = np.where(gsm.isin(pd_ids), "PD", np.where(gsm.isin(hc_ids), "HC", "excluded"))
nkc = full[(full.obs["cell type"] == "NK cell") & (full.obs.condition != "excluded")].copy()
print("\nNK cells per sample/condition:\n", nkc.obs.groupby(["condition", "Sample"], observed=True).size())

sc.pp.filter_genes(nkc, min_cells=10)
sc.tl.rank_genes_groups(nkc, "condition", groups=["PD"], reference="HC", method="wilcoxon", pts=True)
de = sc.get.rank_genes_groups_df(nkc, "PD")

# pseudobulk sanity check: per-sample mean expression, PD vs HC (Mann-Whitney)
# samples with < 5 NK cells are dropped
counts = nkc.obs.Sample.value_counts()
keep = counts[counts >= 5].index
X = nkc[nkc.obs.Sample.isin(keep)]
pb = pd.DataFrame(X.X.toarray() if hasattr(X.X, "toarray") else X.X,
                  columns=X.var_names).groupby(X.obs.Sample.values).mean()
cond = X.obs.groupby("Sample", observed=True).condition.first().loc[pb.index]
p_pb = {}
for g in de.names:
    a, b = pb.loc[cond == "PD", g], pb.loc[cond == "HC", g]
    p_pb[g] = stats.mannwhitneyu(a, b).pvalue
de["pseudobulk_p"] = de.names.map(p_pb)
de["n_PD_samples"] = (cond == "PD").sum(); de["n_HC_samples"] = (cond == "HC").sum()
de.to_csv(OUT / "NK_PD_vs_HC_DE.csv", index=False)

sde = de[(de.pvals_adj < 0.05) & (de.logfoldchanges.abs() > 0.5)]
print("\n== NK PD vs HC: sig genes (cell-level) =", len(sde),
      "| samples used in pseudobulk PD/HC:", (cond == "PD").sum(), (cond == "HC").sum())
cols = ["names", "logfoldchanges", "pvals_adj", "pct_nz_group", "pseudobulk_p"]
print(sde[cols].round(4).to_string(index=False))
print("\nTop 15 each direction by rank (nominal):")
print("UP in PD:"); print(de[de.logfoldchanges > 0].sort_values("scores", ascending=False)[cols].head(15).round(4).to_string(index=False))
print("DOWN in PD:"); print(de[de.logfoldchanges < 0].sort_values("scores")[cols].head(15).round(4).to_string(index=False))

# ---------- figure ----------
sc.settings.figdir = str(ROOT / "figures")
top = specific.names.head(20).tolist()
sc.pl.dotplot(full, top, groupby="cell type", standard_scale="var", show=False,
              save="_NK_markers.pdf")
