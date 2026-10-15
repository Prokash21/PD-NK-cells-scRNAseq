"""Overview figure: CSF cell-type composition per donor + NK specificity of candidate genes.
Reads only obs from the h5ad (h5py), so scanpy is not needed."""
import warnings; warnings.filterwarnings("ignore")
from pathlib import Path
import h5py, numpy as np, pandas as pd
import matplotlib as mpl, matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]   # PD_CSF_NK/
DATA = ROOT / "data"                         # h5ad files from 01_qc_integration_annotation.ipynb
FIG, RES = ROOT / "figures", ROOT / "results"

# ---------------- style (same as 03_nk_figures.py) ----------------
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
C_HC, C_PD, C_NK = "#2a78d6", "#eb6834", "#2a78d6"
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

def panel(ax, letter, x=-0.18, y=1.06):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", color=INK)

# ---------------- data ----------------
with h5py.File(DATA / "integrated.h5ad", "r") as f:
    def cat(k):
        g = f["obs"][k]
        cats = [c.decode() if isinstance(c, bytes) else c for c in g["categories"][:]]
        return pd.Categorical.from_codes(g["codes"][:], cats)
    obs = pd.DataFrame({"cell_type": cat("cell type"), "Sample": cat("Sample")})

pd_ids = {"GSM4208766", "GSM4208768", "GSM4208769", "GSM4208770", "GSM4208771", "GSM4404053"}
hc_ids = {f"GSM42087{i}" for i in range(72, 81)}
gsm = obs.Sample.astype(str).str.split("_").str[0]
obs["condition"] = np.select([gsm.isin(pd_ids), gsm.isin(hc_ids)], ["PD", "HC"], "excluded")
obs = obs[obs.condition != "excluded"]

ct_order = ["CD4+ T-cell", "CD8+ T-cell", "Macrophage", "NK cell", "DC", "Monocyte", "pDC", "B-cell"]
frac = pd.crosstab(obs.Sample, obs.cell_type).reindex(columns=ct_order, fill_value=0)
frac = frac.div(frac.sum(1), axis=0) * 100
cond = obs.groupby("Sample", observed=True).condition.first()
donor = frac.index.astype(str).str.split("_").str[1]
order = sorted(frac.index, key=lambda s: (cond[s], -frac.loc[s, "NK cell"]))
frac.to_csv(RES / "celltype_fraction_per_donor.csv")

genes = pd.read_csv(RES / "NK_markers_vs_all_and_CD8T.csv").head(20)
cand = set(pd.read_csv(RES / "candidate_genes.csv").gene)
for g in cand:                                   # make sure every candidate is shown
    if g not in set(genes.names):
        genes = pd.concat([genes, pd.read_csv(RES / "NK_markers_vs_all_and_CD8T.csv").query("names == @g")])
genes = genes.sort_values("pct_NK")

# ---------------- figure ----------------
fig = plt.figure(figsize=(180 * MM, 92 * MM))
gs = fig.add_gridspec(1, 2, width_ratios=[1.35, 1], wspace=0.55)

# a: composition per donor
ax = fig.add_subplot(gs[0])
is_pd = np.array([cond[s] == "PD" for s in order])
y = np.arange(len(order)) + is_pd * 0.8           # small gap between HC and PD blocks
left = np.zeros(len(order))
for ct, col in zip(ct_order, CT_PAL):
    v = frac.loc[order, ct].values
    ax.barh(y, v, left=left, color=col, height=0.78, edgecolor="white", linewidth=0.6, label=ct)
    left += v
ax.set_yticks(y, [s.split("_")[1] for s in order])
ax.invert_yaxis(); ax.set_xlim(0, 100); ax.set_xlabel("% of CSF immune cells")
for grp, m in (("HC", ~is_pd), ("PD", is_pd)):
    ax.plot([102, 102], [y[m].min() - 0.35, y[m].max() + 0.35], color=INK2, lw=0.8, clip_on=False)
    ax.text(104, y[m].mean(), f"{grp}\n(n = {m.sum()})", va="center", ha="left", color=INK, fontsize=6.5)
ax.set_title("Cell-type composition per donor")
ax.legend(ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.13), handlelength=1, columnspacing=1)
ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0)
panel(ax, "a", -0.14)

# b: dumbbell, % NK vs % other cells
ax = fig.add_subplot(gs[1])
yy = np.arange(len(genes))
pn, pr = genes.pct_NK.values * 100, genes.pct_rest.values * 100
ax.hlines(yy, pr, pn, color=GRID, linewidth=2, zorder=1)
ax.scatter(pr, yy, s=16, color=MUTED, edgecolor="white", linewidth=0.8, zorder=2, label="other CSF cells")
ax.scatter(pn, yy, s=16, color=C_NK, edgecolor="white", linewidth=0.8, zorder=3, label="NK cells")
ax.set_yticks(yy, genes.names)
for t in ax.get_yticklabels():
    t.set_fontstyle("italic")
    if t.get_text() in cand:
        t.set_fontweight("bold"); t.set_color(INK)
    else:
        t.set_color(INK2)
ax.set_xlim(0, 100); ax.set_xlabel("% of cells expressing")
ax.grid(axis="x", color=GRID, linewidth=0.5); ax.set_axisbelow(True)
ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0)
ax.set_title("Top NK-specific genes\n(candidates in bold)")
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2, handletextpad=0.3)
panel(ax, "b", -0.3)

for ext in ("pdf", "png"):
    fig.savefig(FIG / f"Fig4_composition_candidates.{ext}", bbox_inches="tight", facecolor="white")
print(frac.loc[order].round(1).to_string())
