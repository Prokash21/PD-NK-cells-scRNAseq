"""NK-cell characterisation in CSF (GSE141578): publication figures + stats.

Outputs -> figures/NK_Fig*.pdf|png and results/*.csv
"""
import warnings; warnings.filterwarnings("ignore")
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy import stats
from statsmodels.stats.multitest import multipletests
from adjustText import adjust_text

ROOT = Path(__file__).resolve().parents[1]   # PD_CSF_NK/
DATA = ROOT / "data"                         # h5ad files from 01_qc_integration_annotation.ipynb
FIG, RES = ROOT / "figures", ROOT / "results"
FIG.mkdir(exist_ok=True); RES.mkdir(exist_ok=True)

# ---------------- style ----------------
INK, INK2, MUTED, GRID, BG_CELL = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#d9d8d3"
C_HC, C_PD = "#2a78d6", "#eb6834"
C_NK = "#2a78d6"
C_STATE = {"NK SELL+IL7R+": "#1baf7a", "NK CCL5+ effector": "#4a3aa7"}
CT_PAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "axes.titlesize": 8, "axes.labelsize": 7, "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5, "legend.fontsize": 6.5, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5,
    "ytick.major.size": 2.5, "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 600, "figure.dpi": 150,
    "legend.frameon": False,
})
MM = 1 / 25.4

def panel(ax, letter, x=-0.18, y=1.13):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", color=INK)

def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"{name}.{ext}", bbox_inches="tight", facecolor="white")
    plt.close(fig)

def pstr(p):
    return "p < 0.001" if p < 0.001 else f"p = {p:.3f}" if p < 0.01 else f"p = {p:.2f}"

def cond_box(ax, df, y, title, ylabel, p=None):
    """Per-sample box + jittered points, HC vs PD."""
    rng = np.random.default_rng(1)
    for i, (c, col) in enumerate([("HC", C_HC), ("PD", C_PD)]):
        v = df.loc[df.condition == c, y].values
        ax.boxplot(v, positions=[i], widths=0.55, showfliers=False, patch_artist=True,
                   boxprops=dict(facecolor="white", edgecolor=col, lw=1),
                   medianprops=dict(color=col, lw=1.6), whiskerprops=dict(color=col, lw=0.8),
                   capprops=dict(color=col, lw=0.8))
        ax.scatter(i + rng.uniform(-0.14, 0.14, len(v)), v, s=14, color=col,
                   edgecolor="white", lw=0.5, zorder=3)
    ax.set_xticks([0, 1], [f"HC\n(n={(df.condition=='HC').sum()})", f"PD\n(n={(df.condition=='PD').sum()})"])
    ax.set_xlim(-0.6, 1.6); ax.set_title(title, color=INK); ax.set_ylabel(ylabel)
    if p is None:
        p = stats.mannwhitneyu(df.loc[df.condition == "HC", y], df.loc[df.condition == "PD", y]).pvalue
    lo, hi = ax.get_ylim(); h = hi + (hi - lo) * 0.04
    ax.plot([0, 0, 1, 1], [h, h + (hi-lo)*0.03, h + (hi-lo)*0.03, h], color=INK2, lw=0.6)
    ax.text(0.5, h + (hi - lo) * 0.05, pstr(p), ha="center", va="bottom", fontsize=6.5, color=INK)
    ax.set_ylim(lo, h + (hi - lo) * 0.18)
    return p

# ---------------- data ----------------
adata = sc.read_h5ad(DATA / "integrated.h5ad")
full = adata.raw.to_adata(); full.obs = adata.obs.copy()
for k in ("X_umap", "X_scVI"): full.obsm[k] = adata.obsm[k]
pd_ids = {"GSM4208766", "GSM4208768", "GSM4208769", "GSM4208770", "GSM4208771", "GSM4404053"}
hc_ids = {f"GSM42087{i}" for i in range(72, 81)}
gsm = full.obs.Sample.astype(str).str.split("_").str[0]
full.obs["condition"] = np.select([gsm.isin(pd_ids), gsm.isin(hc_ids)], ["PD", "HC"], "excluded")
full.obs["is_NK"] = full.obs["cell type"] == "NK cell"
ct_order = ["CD4+ T-cell", "CD8+ T-cell", "NK cell", "B-cell", "Monocyte", "Macrophage", "DC", "pDC"]
full.obs["cell type"] = pd.Categorical(full.obs["cell type"], categories=ct_order)
ct_col = dict(zip(ct_order, CT_PAL))
LABEL_OFF = {"Macrophage": (-0.6, 0.9, "right"), "DC": (0.6, 0.6, "left"), "Monocyte": (0, -1.0, "center"),
             "pDC": (0.4, 0.5, "left"), "NK cell": (-0.3, 0, "right"), "B-cell": (0.5, 0, "left")}
MIN_NK = 10   # min NK cells per sample for per-sample NK statistics

# =========================================================================
# FIGURE 1 — NK identity in CSF
# =========================================================================
sc.tl.rank_genes_groups(full, "cell type", groups=["NK cell"], reference="rest", method="wilcoxon", pts=True)
vol = sc.get.rank_genes_groups_df(full, "NK cell")
vol["nlp"] = -np.log10(vol.pvals_adj.clip(lower=1e-300))

fig = plt.figure(figsize=(180 * MM, 150 * MM))
gs = GridSpec(2, 3, figure=fig, width_ratios=[1, 1, 1.05], height_ratios=[1, 1.05], wspace=0.35, hspace=0.45)
um = full.obsm["X_umap"]

ax = fig.add_subplot(gs[0, 0])
order = np.random.default_rng(0).permutation(full.n_obs)
sizes = np.where(full.obs["cell type"].isin(["CD4+ T-cell", "CD8+ T-cell"]).values, 0.4, 1.6)
ax.scatter(um[order, 0], um[order, 1], s=sizes[order], c=full.obs["cell type"].map(ct_col).astype(str).values[order],
           lw=0, rasterized=True)
for ct in ct_order:
    m = (full.obs["cell type"] == ct).values
    x, y = np.median(um[m], axis=0)
    dx, dy, ha = LABEL_OFF.get(ct, (0, 0, "center"))
    ax.text(x + dx, y + dy, ct.replace(" T-cell", " T").replace("-cell", ""), fontsize=6, ha=ha, va="center",
            color=INK, fontweight="bold", path_effects=[mpl.patheffects.withStroke(linewidth=2, foreground="white")])
ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
ax.set_title(f"CSF immune cells (n = {full.n_obs:,})"); panel(ax, "a", -0.12)

ax = fig.add_subplot(gs[0, 1])
m = full.obs.is_NK.values
ax.scatter(um[~m, 0], um[~m, 1], s=0.4, c=BG_CELL, lw=0, rasterized=True)
ax.scatter(um[m, 0], um[m, 1], s=1.5, c=C_NK, lw=0, rasterized=True)
ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
ax.set_title(f"NK cells (n = {m.sum()}, {100*m.mean():.1f}%)"); panel(ax, "b", -0.12)

ax = fig.add_subplot(gs[0, 2])
up = (vol.pvals_adj < 0.05) & (vol.logfoldchanges > 1) & (vol.pct_nz_group > 0.25)
dn = (vol.pvals_adj < 0.05) & (vol.logfoldchanges < -1) & (vol.pct_nz_group > 0.25)
ns = ~(up | dn)
ax.scatter(vol.logfoldchanges[ns], vol.nlp[ns], s=2, c=BG_CELL, lw=0, rasterized=True)
ax.scatter(vol.logfoldchanges[dn], vol.nlp[dn], s=4, c=MUTED, lw=0, rasterized=True)
ax.scatter(vol.logfoldchanges[up], vol.nlp[up], s=5, c=C_NK, lw=0, rasterized=True)
lab = vol[up].sort_values("nlp", ascending=False).head(14)
ymax = vol.nlp.max()
for k, r in enumerate(lab.itertuples()):
    ty = ymax * (1.0 - k * 0.062)
    ax.plot([r.logfoldchanges, 12.3], [r.nlp, ty], color=MUTED, lw=0.35)
    ax.text(12.5, ty, r.names, fontsize=5.5, color=INK, style="italic", va="center")
ax.axvline(1, color=GRID, lw=0.6, ls="--"); ax.axvline(-1, color=GRID, lw=0.6, ls="--")
ax.set_xlim(-8, 16); ax.set_xticks([-5, 0, 5, 10])
ax.set_xlabel("log$_2$ fold change (NK vs other cells)"); ax.set_ylabel("−log$_{10}$ adjusted p")
ax.set_title(f"{up.sum()} NK-enriched genes"); panel(ax, "c")

# d: curated dotplot
groups = {
    "Cytotoxic": ["NKG7", "GNLY", "PRF1", "GZMA", "GZMB", "CTSW", "CST7"],
    "NK receptors": ["KLRD1", "KLRF1", "KLRC1", "NCR1", "NCR3", "SH2D1B", "KIR2DL4"],
    "Adaptors": ["TYROBP", "FCER1G", "IL2RB"],
    "Chemokines": ["XCL1", "XCL2", "CCL5"],
    "TF / identity": ["ID2", "HOPX", "NCAM1", "FCGR3A"],
    "T-cell": ["CD3D", "CD3E", "TRDC", "CD8A"],
}
genes = [g for v in groups.values() for g in v if g in full.var_names]
ax = fig.add_subplot(gs[1, :])
X = full[:, genes].X.toarray()
ct = full.obs["cell type"].values
frac = pd.DataFrame(X > 0, columns=genes).groupby(ct, observed=True).mean().loc[ct_order]
mean = pd.DataFrame(X, columns=genes).groupby(ct, observed=True).mean().loc[ct_order]
scaled = (mean - mean.min()) / (mean.max() - mean.min())
cmap = mpl.colors.LinearSegmentedColormap.from_list("seq", ["#f0efec", "#9ec5f4", "#3987e5", "#1c5cab", "#0d366b"])
yy, xx = np.meshgrid(range(len(ct_order)), range(len(genes)), indexing="ij")
ax.scatter(xx.ravel(), yy.ravel(), s=(frac.values.ravel() * 60) + 0.5, c=scaled.values.ravel(),
           cmap=cmap, vmin=0, vmax=1, edgecolor=INK2, lw=0.2)
ax.set_xticks(range(len(genes)), genes, rotation=90, style="italic")
ax.set_yticks(range(len(ct_order)), ct_order); ax.invert_yaxis()
ax.set_xlim(-0.7, len(genes) - 0.3); ax.set_ylim(len(ct_order) - 0.4, -1.9)
ax.spines[["left", "bottom"]].set_visible(False); ax.tick_params(length=0)
ax.axhspan(1.5, 2.5, color="#eaf2fc", zorder=0)
i = 0
for gname, gl in groups.items():
    n = len([g for g in gl if g in full.var_names])
    ax.plot([i - 0.3, i + n - 0.7], [-1.1, -1.1], color=INK2, lw=0.8)
    ax.text(i + (n - 1) / 2, -1.3, gname, ha="center", va="bottom", fontsize=6.5, color=INK)
    i += n
cax = ax.inset_axes([1.01, 0.55, 0.012, 0.4])
cb = fig.colorbar(mpl.cm.ScalarMappable(cmap=cmap), cax=cax); cb.set_label("Scaled mean\nexpression", fontsize=6)
cb.outline.set_linewidth(0.4); cb.ax.tick_params(labelsize=5.5)
for j, f in enumerate([0.25, 0.5, 0.75, 1.0]):
    ax.scatter([], [], s=f * 60 + 0.5, c="white", edgecolor=INK2, lw=0.4, label=f"{int(f*100)}%")
ax.legend(title="% expressing", loc="lower left", bbox_to_anchor=(1.0, -0.05), fontsize=5.5, title_fontsize=6,
          labelspacing=0.9, borderaxespad=0)
panel(ax, "d", -0.07, 1.02)
save(fig, "NK_Fig1_identity")
vol[up].sort_values("scores", ascending=False).to_csv(RES / "NK_markers_vs_all.csv", index=False)

# =========================================================================
# FIGURE 2 — NK heterogeneity
# =========================================================================
nk = full[full.obs.is_NK].copy()
sc.pp.neighbors(nk, use_rep="X_scVI", n_neighbors=15, random_state=0)
sc.tl.umap(nk, random_state=0)
sc.tl.leiden(nk, resolution=0.5, key_added="sub", flavor="igraph", n_iterations=2, random_state=0)
# name clusters by marker content (SELL high = CD56bright-like)
sell = pd.Series(nk[:, "SELL"].X.toarray().ravel(), index=nk.obs_names).groupby(nk.obs["sub"]).mean()
bright = sell.idxmax()
nk.obs["NK state"] = np.where(nk.obs["sub"] == bright, "NK SELL+IL7R+", "NK CCL5+ effector")
states = list(C_STATE)
nk.obs["NK state"] = pd.Categorical(nk.obs["NK state"], categories=states)
sc.tl.rank_genes_groups(nk, "NK state", method="wilcoxon", pts=True)
state_de = pd.concat([sc.get.rank_genes_groups_df(nk, s).assign(state=s) for s in states])
state_de.to_csv(RES / "NK_state_markers.csv", index=False)

fig = plt.figure(figsize=(180 * MM, 110 * MM))
gs = GridSpec(2, 4, figure=fig, width_ratios=[1.1, 1.3, 0.8, 0.8], wspace=0.6, hspace=0.6)
ax = fig.add_subplot(gs[:, 0])
u = nk.obsm["X_umap"]
for s in states:
    mm = (nk.obs["NK state"] == s).values
    ax.scatter(u[mm, 0], u[mm, 1], s=5, c=C_STATE[s], lw=0, label=f"{s} (n={mm.sum()})", rasterized=True)
ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.06), markerscale=2, handletextpad=0.2)
ax.set_title("NK sub-clustering (scVI latent)"); panel(ax, "a", -0.12)

sgenes = ["SELL", "IL7R", "TCF7", "GPR183", "XCL1", "KLRC1", "CCL5", "CCL4", "GZMH", "GZMA", "GZMM",
          "CD160", "KRT81", "KRT86", "PRF1", "FCGR3A"]
sgenes = [g for g in sgenes if g in nk.var_names]
ax = fig.add_subplot(gs[:, 1])
X = nk[:, sgenes].X.toarray()
fr = pd.DataFrame(X > 0, columns=sgenes).groupby(nk.obs["NK state"].values, observed=True).mean().loc[states]
mn = pd.DataFrame(X, columns=sgenes).groupby(nk.obs["NK state"].values, observed=True).mean().loc[states]
sm = (mn - mn.min()) / (mn.max() - mn.min() + 1e-9)
yy, xx = np.meshgrid(range(len(sgenes)), range(len(states)), indexing="ij")
ax.scatter(xx.ravel(), yy.ravel(), s=fr.T.values.ravel() * 70 + 0.5, c=sm.T.values.ravel(), cmap=cmap,
           vmin=0, vmax=1, edgecolor=INK2, lw=0.2)
ax.set_yticks(range(len(sgenes)), sgenes, style="italic"); ax.invert_yaxis()
ax.set_xticks(range(len(states)), [s.replace("NK ", "") for s in states], rotation=30, ha="right")
ax.set_xlim(-0.6, len(states) - 0.4); ax.tick_params(length=0)
ax.spines[["left", "bottom"]].set_visible(False)
ax.set_title("State markers"); panel(ax, "b", -0.6)

# c: CD3 complex / TCR transcripts — NK vs T
cg = ["CD3E", "CD3D", "CD3G", "TRDC", "TRAC", "TRBC2"]
cg = [g for g in cg if g in full.var_names]
ax = fig.add_subplot(gs[0, 2:])
sub = full[full.obs["cell type"].isin(["CD4+ T-cell", "CD8+ T-cell", "NK cell"])]
pct = pd.DataFrame(sub[:, cg].X.toarray() > 0, columns=cg).groupby(sub.obs["cell type"].values, observed=True).mean().loc[["CD4+ T-cell", "CD8+ T-cell", "NK cell"]] * 100
w = 0.26
for k, (c, col) in enumerate(zip(pct.index, ["#86b6ef", "#256abf", C_PD])):
    ax.bar(np.arange(len(cg)) + (k - 1) * w, pct.loc[c], width=w - 0.03, color=col, label=c.replace("-cell", ""))
ax.set_xticks(range(len(cg)), cg, style="italic"); ax.set_ylabel("% cells detected"); ax.set_ylim(0, 105)
ax.legend(ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.22), handlelength=1, columnspacing=0.8)
ax.yaxis.grid(True, color=GRID, lw=0.4); ax.set_axisbelow(True)
panel(ax, "c", -0.13, 1.1)

# d: state composition by condition
nkc = nk[nk.obs.condition != "excluded"]
comp = pd.crosstab(nkc.obs.Sample, nkc.obs["NK state"])
comp = comp[comp.sum(1) >= MIN_NK]
compdf = pd.DataFrame({"frac_eff": comp["NK CCL5+ effector"] / comp.sum(1) * 100})
compdf["condition"] = nkc.obs.groupby("Sample", observed=True).condition.first().loc[compdf.index].values
ax = fig.add_subplot(gs[1, 2])
p_state = cond_box(ax, compdf, "frac_eff", "CCL5+ effector NK", "% of NK cells")
panel(ax, "d", -0.45)

# e: NK frequency per sample
allc = full.obs[full.obs.condition != "excluded"]
freq = allc.groupby("Sample", observed=True).agg(nk=("is_NK", "sum"), tot=("is_NK", "size"), condition=("condition", "first"))
freq["pct"] = freq.nk / freq.tot * 100
ax = fig.add_subplot(gs[1, 3])
p_freq = cond_box(ax, freq, "pct", "NK frequency", "% of CSF immune cells")
panel(ax, "e", -0.45)
save(fig, "NK_Fig2_heterogeneity")
freq.to_csv(RES / "NK_frequency_per_sample.csv"); compdf.to_csv(RES / "NK_state_fraction_per_sample.csv")

# =========================================================================
# FIGURE 3 — PD vs HC functional programmes + pseudobulk DE
# =========================================================================
modules = {
    "Cytotoxicity": ["PRF1", "GZMA", "GZMB", "GZMH", "GZMM", "GNLY", "NKG7", "CTSW", "CST7"],
    "Activating receptors": ["NCR1", "NCR3", "KLRK1", "KLRC2", "CD226", "FCGR3A", "KLRF1", "SH2D1B", "CD160", "CD244"],
    "Inhibitory receptors": ["KLRC1", "KIR2DL1", "KIR2DL3", "KIR3DL1", "KIR3DL2", "LILRB1", "TIGIT", "HAVCR2", "PDCD1", "KLRG1"],
    "IFN-I response": ["ISG15", "IFI6", "IFIT1", "IFIT3", "MX1", "OAS1", "STAT1", "IFI44L", "IRF7", "XAF1"],
    "Chemokines / cytokines": ["CCL3", "CCL4", "CCL5", "XCL1", "XCL2", "IFNG", "TNF", "CCL4L2"],
    "Tissue residency": ["CD69", "ITGA1", "ITGAE", "CXCR6", "RGS1", "CD9"],
}
for k, g in modules.items():
    sc.tl.score_genes(nkc, [x for x in g if x in nkc.var_names], score_name=k, random_state=0)
keep = nkc.obs.Sample.value_counts(); keep = keep[keep >= MIN_NK].index
ms = nkc.obs[nkc.obs.Sample.isin(keep)].groupby("Sample", observed=True)[list(modules) + ["condition"]] \
    .agg({**{k: "mean" for k in modules}, "condition": "first"})
praw = [stats.mannwhitneyu(ms.loc[ms.condition == "HC", k], ms.loc[ms.condition == "PD", k]).pvalue for k in modules]
padj = multipletests(praw, method="fdr_bh")[1]
modstats = pd.DataFrame({"module": list(modules), "p": praw, "p_BH": padj,
                         "median_HC": [ms.loc[ms.condition == "HC", k].median() for k in modules],
                         "median_PD": [ms.loc[ms.condition == "PD", k].median() for k in modules]})
modstats.to_csv(RES / "NK_module_scores_PD_vs_HC.csv", index=False); ms.to_csv(RES / "NK_module_scores_per_sample.csv")

# pseudobulk DESeq2 on raw counts
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
comb = sc.read_h5ad(DATA / "combined.h5ad")
comb.var_names_make_unique()
cells = nkc.obs_names[nkc.obs.Sample.isin(keep)]
assert set(cells) <= set(comb.obs_names), "barcode mismatch with combined.h5ad"
cnt = comb[cells]
samp = nkc.obs.loc[cells, "Sample"].astype(str).values
M = cnt.X.tocsr() if hasattr(cnt.X, "tocsr") else cnt.X
pb = pd.DataFrame({s: np.asarray(M[samp == s].sum(0)).ravel() for s in np.unique(samp)}, index=cnt.var_names).T
pb = pb.loc[:, (pb >= 10).sum(0) >= 3].astype(int)
meta = pd.DataFrame({"condition": nkc.obs.groupby("Sample", observed=True).condition.first().loc[pb.index].values}, index=pb.index)
dds = DeseqDataSet(counts=pb, metadata=meta, design="~condition", quiet=True)
dds.deseq2()
st = DeseqStats(dds, contrast=["condition", "PD", "HC"], quiet=True); st.summary()
pbde = st.results_df.dropna(subset=["padj"]).sort_values("pvalue")
pbde.to_csv(RES / "NK_pseudobulk_DESeq2_PD_vs_HC.csv")

fig = plt.figure(figsize=(180 * MM, 105 * MM))
gs = GridSpec(2, 5, figure=fig, width_ratios=[1, 1, 1, 0.25, 2.1], wspace=0.75, hspace=0.75)
for i, (k, p, q) in enumerate(zip(modules, praw, padj)):
    ax = fig.add_subplot(gs[i // 3, i % 3])
    cond_box(ax, ms.reset_index(), k, k, "Module score" if i % 3 == 0 else "", p=p)
    if i == 0: panel(ax, "a", -0.6)
ax = fig.add_subplot(gs[:, 4])
v = pbde.copy(); v["nlp"] = -np.log10(v.pvalue)
sig = v.padj < 0.1
ax.scatter(v.log2FoldChange[~sig], v.nlp[~sig], s=3, c=BG_CELL, lw=0, rasterized=True)
ax.scatter(v.log2FoldChange[sig & (v.log2FoldChange > 0)], v.nlp[sig & (v.log2FoldChange > 0)], s=7, c=C_PD, lw=0)
ax.scatter(v.log2FoldChange[sig & (v.log2FoldChange < 0)], v.nlp[sig & (v.log2FoldChange < 0)], s=7, c=C_HC, lw=0)
top = v.head(14)
texts = [ax.text(r.log2FoldChange, r.nlp, r.Index, fontsize=5.5, style="italic", color=INK) for r in top.itertuples()]
adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.4))
ax.axvline(0, color=GRID, lw=0.6)
lim = np.nanpercentile(v.log2FoldChange.abs(), 99.8); ax.set_xlim(-lim * 1.1, lim * 1.1)
ax.set_xlabel("log$_2$ fold change (PD / HC)"); ax.set_ylabel("−log$_{10}$ p")
ax.set_title(f"NK pseudobulk DESeq2 ({(meta.condition=='PD').sum()} PD vs {(meta.condition=='HC').sum()} HC)\n"
             f"{sig.sum()} genes FDR < 0.1")
ax.text(0.02, 0.02, "↓ in PD", transform=ax.transAxes, color=C_HC, va="bottom", fontsize=6.5, fontweight="bold")
ax.text(0.98, 0.02, "↑ in PD", transform=ax.transAxes, color=C_PD, va="bottom", ha="right", fontsize=6.5, fontweight="bold")
panel(ax, "b", -0.2, 1.1)
save(fig, "NK_Fig3_PD_vs_HC")

# =========================================================================
# FIGURE 4 — pathway enrichment
# =========================================================================
import gseapy as gp
enr_ok = True
try:
    bg = full.var_names.tolist()
    nk_up = vol[up].names.tolist()
    e1 = gp.enrichr(gene_list=nk_up, gene_sets=["GO_Biological_Process_2023"], organism="human",
                    background=bg, outdir=None).results
    e1.to_csv(RES / "NK_markers_GO_BP_enrichr.csv", index=False)
    rnk = pbde["stat"].dropna().sort_values(ascending=False)
    pre = gp.prerank(rnk=rnk, gene_sets="MSigDB_Hallmark_2020", min_size=10, max_size=500,
                     permutation_num=2000, seed=0, outdir=None, threads=4).res2d
    pre.to_csv(RES / "NK_PD_vs_HC_GSEA_hallmark.csv", index=False)
except Exception as ex:
    enr_ok = False; print("Enrichment skipped:", ex)

if enr_ok:
    fig, axes = plt.subplots(1, 2, figsize=(180 * MM, 70 * MM), gridspec_kw=dict(width_ratios=[1, 1], wspace=1.25))
    ax = axes[0]
    t = e1.sort_values("Adjusted P-value").head(12).iloc[::-1]
    t["term"] = t.Term.str.replace(r"\s*\(GO:\d+\)", "", regex=True)
    t["term"] = t.term.apply(lambda s: s if len(s) < 48 else s[:46] + "…")
    ax.barh(t.term, -np.log10(t["Adjusted P-value"]), color=C_NK, height=0.7)
    ax.set_xlabel("−log$_{10}$ adjusted p"); ax.set_title(f"GO BP: NK-enriched genes (n = {len(nk_up)})")
    ax.tick_params(axis="y", length=0); ax.xaxis.grid(True, color=GRID, lw=0.4); ax.set_axisbelow(True)
    panel(ax, "a", -1.1)
    ax = axes[1]
    pre["NES"] = pre.NES.astype(float); pre["FDR q-val"] = pre["FDR q-val"].astype(float)
    t = pre.reindex(pre.NES.abs().sort_values(ascending=False).index).head(12).sort_values("NES")
    t["Term"] = t.Term.str.replace("HALLMARK_", "")
    cols = [C_PD if n > 0 else C_HC for n in t.NES]
    ax.barh(t.Term, t.NES, color=cols, height=0.7)
    for y, (n, q) in enumerate(zip(t.NES, t["FDR q-val"])):
        ax.text(n + (0.05 if n > 0 else -0.05), y, f"q={q:.2f}", va="center", ha="left" if n > 0 else "right",
                fontsize=5.5, color=INK2)
    ax.axvline(0, color=INK2, lw=0.6); ax.set_xlabel("NES (PD vs HC)")
    lim = t.NES.abs().max() * 1.45; ax.set_xlim(-lim, lim)
    ax.set_title("GSEA Hallmark: NK cells, PD vs HC"); ax.tick_params(axis="y", length=0)
    panel(ax, "b", -1.1)
    save(fig, "NK_Fig4_pathways")

# ---------------- console summary ----------------
print("NK states:", nk.obs["NK state"].value_counts().to_dict())
print("Top state markers:"); [print(" ", s, state_de[state_de.state == s].names.head(12).tolist()) for s in states]
print("CD3/TCR detection (%):\n", pct.round(1).to_string())
print(f"NK freq: HC median {freq[freq.condition=='HC'].pct.median():.2f}% vs PD {freq[freq.condition=='PD'].pct.median():.2f}%, {pstr(p_freq)}")
print(f"Effector state: HC median {compdf[compdf.condition=='HC'].frac_eff.median():.1f}% vs PD {compdf[compdf.condition=='PD'].frac_eff.median():.1f}%, {pstr(p_state)}")
print("Samples per condition (>= MIN_NK):", ms.condition.value_counts().to_dict())
print(modstats.round(4).to_string(index=False))
print("DESeq2 FDR<0.1:", (pbde.padj < 0.1).sum(), "| FDR<0.05:", (pbde.padj < 0.05).sum())
print(pbde.head(15)[["baseMean", "log2FoldChange", "pvalue", "padj"]].round(4).to_string())
if enr_ok:
    print(e1.sort_values("Adjusted P-value").head(10)[["Term", "Adjusted P-value", "Genes"]].to_string(index=False))
    print(pre.sort_values("FDR q-val").head(10)[["Term", "NES", "NOM p-val", "FDR q-val"]].to_string(index=False))
