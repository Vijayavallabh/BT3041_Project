"""Step 2: quality control and exploratory analysis (no labels are used to build any embedding).

Outputs: figures/fig2-6*, figures/figS1-S5*, results/explore/*
Run from the project root:  python scripts/02_explore.py
"""
import json, sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.cluster import DBSCAN
from sklearn.decomposition import PCA
from sklearn.metrics import cohen_kappa_score, roc_auc_score
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.common import (PROC, RES, SEED, GROUP_ORDER, GROUP_COL, PAM50_ORDER, PAM50_COL, MARKERS,
                        INK2, MUTED, style, save, write_params, drop_excluded)

N_TOP, N_PC, UMAP_PCS, UMAP_NN, UMAP_MIN_DIST = 2000, 50, 30, [15, 30], 0.3
N_PERM, N_BOOT = 1000, 1000
PARAMS = {
    "samples": "primary tumours not excluded (see 01 exclusions), HiSeqV2 log2(RSEM upper-quartile normalised + 1)",
    "expressed_gene_filter": "mean log2 expression >= 1 across tumours",
    "gene_exclusions": "sex-specific genes (chrY + XIST/TSIX) removed before variance ranking (12 male patients would otherwise drive a component)",
    "n_top_variance_genes": N_TOP,
    "scaling": "each gene centred and scaled to unit variance (z-score) across tumours",
    "pca_components": N_PC,
    "pc_orientation": "PC1 sign fixed so that ESR1 has a positive loading (ER-high tumours on the right in every PCA figure)",
    "umap": {"input": f"first {UMAP_PCS} PCs", "n_neighbors": UMAP_NN, "min_dist": UMAP_MIN_DIST, "metric": "euclidean", "random_state": SEED},
    "umap_island": "DBSCAN(eps=0.7, min_samples=10) on each UMAP; the island = the DBSCAN cluster with the most PAM50-Basal tumours",
    "marker_tests": "two-sided Mann-Whitney U on float64 values, TNBC vs HR+/HER2- and HER2+ vs HR+/HER2-, BH-FDR over all 16 tests; effect = rank-biserial r = 2U/(n1 n2) - 1 (positive = higher in the first group); supervised set only",
    "mrna_vs_status": f"AUC of ESR1/PGR/ERBB2 mRNA for IHC ER/PR and HER2 (v1 rule), {N_BOOT} bootstrap 95% CI; Spearman reported with its prevalence ceiling sqrt(3p(1-p))",
    "batch_metric": f"eta-squared (one-way ANOVA R^2) of PC1-10 on tissue source site (sites with <10 tumours pooled as 'other'), clinical group, PAM50, histology, stage; chance level = 99th percentile of {N_PERM} label permutations",
    "pca_outlier_rule": "robust z (median, 1.4826*MAD) beyond 5 on any of PC1-5; flagged, not removed",
    "heatmap_pam50": "z-score per gene, clipped to +/-3; rows ordered by average linkage on correlation; columns ordered by PAM50 then clinical group",
    "rppa_spearman": "total-protein antibodies only (phospho/cleaved excluded), primary tumours with both RNA and protein; pan-isoform antibodies kept in the table but left out of the median",
    "star_crosscheck": "same PCA pipeline on STAR log2(CPM+1) for the same tumours; per-PC |Pearson| between score vectors",
    "tumour_normal_umap": "not made: the tumour+normal PCA (figS3) already separates normals; recorded here as a deliberate omission",
    "figS1b_sites": "the six sites (pooled factor, >=10 tumours) with the largest share of the site eta-squared on the most site-associated PC",
    "erbb2_auc_by_source": "ERBB2 mRNA AUC for HER2 status computed separately for each HER2_source (FISH, IHC score, IHC status text)",
    "literature_like_for_like": "HER2E share among clinical HER2+ recomputed with the TCGA 2012 HER2 labels and PAM50 calls (BRCA_clinicalMatrix)",
}
OUT = RES / "explore"
OUT.mkdir(parents=True, exist_ok=True)
style()
rep = {}
rng = np.random.default_rng(SEED)


def log(k, v):
    rep[k] = v
    print(f"{k}: {v}", flush=True)


lab = drop_excluded(pd.read_csv(PROC / "tcga_labels.tsv", sep="\t", index_col=0))
X = pd.read_parquet(PROC / "tcga_hiseq_tumour_log2.parquet").loc[lab.index].astype("float64")
N = pd.read_parquet(PROC / "tcga_hiseq_normal_log2.parquet").astype("float64")
sexg = set(pd.read_csv(PROC / "sex_specific_genes.txt", header=None)[0])
sup = lab.supervised_set
log("n_tumours", len(X))
log("clinical_group_all", lab.clinical_group.value_counts().to_dict())
log("clinical_group_supervised", lab[sup].clinical_group.value_counts().to_dict())


def top_var(M, n):
    M = M.loc[:, (M.mean() >= 1) & ~M.columns.isin(sexg)]
    return M[M.var().sort_values(ascending=False).index[:n]]


def z(M):
    return (M - M.mean()) / M.std(ddof=0)


Xv = top_var(X, N_TOP)
Z = z(Xv)
pca = PCA(n_components=N_PC, random_state=SEED).fit(Z)
sign = np.sign(pca.components_[0, list(Z.columns).index("ESR1")]) if "ESR1" in Z.columns else 1
S = pd.DataFrame(pca.transform(Z), index=Z.index, columns=[f"PC{i}" for i in range(1, N_PC + 1)])
S["PC1"] *= sign
ev = pca.explained_variance_ratio_
log("pca_variance_explained_PC1_10", [round(float(v), 4) for v in ev[:10]])
S.to_csv(OUT / "pca_scores_tumours.tsv", sep="\t")
load = pd.DataFrame(pca.components_[:5].T, index=Z.columns, columns=[f"PC{i}" for i in range(1, 6)])
load["PC1"] *= sign
load.to_csv(OUT / "pca_loadings_PC1_5.tsv", sep="\t")
log("pca_top_loadings", {pc: {"positive": load[pc].nlargest(10).index.tolist(), "negative": load[pc].nsmallest(10).index.tolist()} for pc in load})

rz = (S.iloc[:, :5] - S.iloc[:, :5].median()) / (1.4826 * (S.iloc[:, :5] - S.iloc[:, :5].median()).abs().median())
log("pca_outliers_flagged", sorted(S.index[rz.abs().max(axis=1) > 5]))


# ---- batch / biology variance decomposition with a permutation chance level
def eta2(y, g):
    g = pd.Series(g)
    m = g.notna().values
    y, gv = np.asarray(y)[m], g.values[m]
    codes, inv = np.unique(gv.astype(str), return_inverse=True)
    mu = y.mean()
    return float(sum((inv == k).sum() * (y[inv == k].mean() - mu) ** 2 for k in range(len(codes))) / ((y - mu) ** 2).sum())


fac = {"Tissue source site": lab.tss.where(lab.tss.map(lab.tss.value_counts()) >= 10, "other"),
       "Clinical group": lab.clinical_group.where(lab.clinical_group != "Incomplete"),
       "PAM50 (PanCan)": lab.PAM50_pancan, "Histology": lab.histology, "Stage": lab.stage}
E, null99 = {}, {}
for k, v in fac.items():
    v = v.reindex(S.index)
    E[k] = [eta2(S[f"PC{i}"].values, v) for i in range(1, 11)]
    null99[k] = float(np.percentile([eta2(S["PC1"].values, rng.permutation(v.values)) for _ in range(N_PERM)], 99))
E = pd.DataFrame(E, index=[f"PC{i}\n({ev[i-1]*100:.1f}%)" for i in range(1, 11)]).T
E.to_csv(OUT / "pc_variance_by_factor_eta2.tsv", sep="\t")
log("eta2_null99", {k: round(v, 3) for k, v in null99.items()})
log("eta2_PC1_PC5", E.iloc[:, :5].round(3).to_dict())
tss_pc = int(np.argmax(E.loc["Tissue source site"].values)) + 1
tss_max = float(E.loc["Tissue source site"].max())

fig, ax = plt.subplots(figsize=(10.5, 3.4))
Eshow = E.copy()
Eshow.index = [f"{k}  (chance ≤ {null99[k]:.2f})" for k in E.index]
sns.heatmap(Eshow, annot=True, fmt=".2f", cmap="Blues", vmin=0, vmax=0.8, cbar_kws={"label": "η² (share of PC variance)"},
            linewidths=2, linecolor="white", ax=ax, annot_kws={"size": 8})
ax.set_title(f"PC1-PC4 follow tumour biology; tissue source site peaks on PC{tss_pc} (η² {tss_max:.2f}, chance ≤ {null99['Tissue source site']:.2f})", fontsize=10.5)
ax.set_xlabel(""); ax.tick_params(axis="x", labelsize=8)
save(fig, "figS1_pc_variance_by_factor")


def scatter_groups(ax, xy, groups, order, cols, s=9, title="", bg_label="Incomplete"):
    groups = pd.Series(np.asarray(groups, dtype=object))
    bg = (groups.isna() | ~groups.isin(order)).values
    if bg.any():
        ax.scatter(xy[bg, 0], xy[bg, 1], s=s * 0.6, c="#d9d8d2", lw=0, label=f"{bg_label} (n={bg.sum()})", rasterized=True)
    for g in order:
        m = (groups == g).values
        ax.scatter(xy[m, 0], xy[m, 1], s=s, c=cols[g], lw=0.3, edgecolors="white", label=f"{g} (n={m.sum()})",
                   rasterized=True, alpha=0.9)
    ax.set_title(title)
    ax.set_xticks([]); ax.set_yticks([])


# ---- PCA tumours: clinical group + ESR1
xy = S[["PC1", "PC2"]].values
er_auc_pc1 = roc_auc_score((lab.ER == "Positive")[lab.ER.notna()], S.PC1[lab.ER.notna()])
log("PC1_AUC_for_IHC_ER", round(float(er_auc_pc1), 3))
fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), gridspec_kw={"width_ratios": [1.25, 1]})
scatter_groups(axes[0], xy, lab.clinical_group, GROUP_ORDER, GROUP_COL, title="Clinical receptor group (IHC/FISH)")
axes[0].legend(loc="upper left", bbox_to_anchor=(1.0, 1.0), fontsize=8, markerscale=1.6, handletextpad=0.2)
axes[0].set_xlabel(f"PC1 ({ev[0]*100:.1f}% of variance)"); axes[0].set_ylabel(f"PC2 ({ev[1]*100:.1f}%)")
sc = axes[1].scatter(xy[:, 0], xy[:, 1], c=X.loc[S.index, "ESR1"], s=9, cmap="viridis", lw=0, rasterized=True)
axes[1].set_title("ESR1 mRNA (log2)"); axes[1].set_xticks([]); axes[1].set_yticks([])
axes[1].set_xlabel("PC1"); axes[1].set_ylabel("PC2")
plt.colorbar(sc, ax=axes[1], shrink=0.7, label="log2 expression")
fig.suptitle(f"PCA of {len(S):,} TCGA breast tumours (top {N_TOP:,} variable genes): PC1 separates IHC ER status (AUC {er_auc_pc1:.2f})",
             fontsize=11, fontweight="bold", x=0.01, ha="left")
fig.tight_layout()
save(fig, "fig2_pca_clinical_group")

# PAM50 as small multiples
fig, axes = plt.subplots(1, 5, figsize=(14, 3.3), sharex=True, sharey=True)
for ax, g in zip(axes, PAM50_ORDER):
    m = (lab.PAM50_pancan.reindex(S.index) == g).values
    ax.scatter(xy[~m, 0], xy[~m, 1], s=4, c="#e1e0d9", lw=0, rasterized=True)
    ax.scatter(xy[m, 0], xy[m, 1], s=7, c=PAM50_COL[g], lw=0, rasterized=True)
    ax.set_title(f"{g} (n={m.sum()})"); ax.set_xticks([]); ax.set_yticks([])
axes[0].set_ylabel("PC2"); axes[2].set_xlabel("PC1")
fig.suptitle("Same PCA, PAM50 subtypes (comparator labels, PanCanAtlas)", fontsize=11, fontweight="bold", x=0.01, ha="left", y=1.08)
save(fig, "figS2_pca_pam50_facets")

# histology and tissue source site on the PCA (plan 2.2 step 6)
# histology uses hues that are not clinical-group colours; sites are the six that contribute most to the site eta2 on that PC
HCOL = {"IDC": "#7a5c3e", "ILC": "#e87ba4", "Other": "#00a3b4"}
site_f = fac["Tissue source site"].reindex(S.index)
pcv = S[f"PC{tss_pc}"]
contrib = pcv.groupby(site_f).agg(lambda v: len(v) * (v.mean() - pcv.mean()) ** 2) / ((pcv - pcv.mean()) ** 2).sum()
top_sites = contrib.drop("other", errors="ignore").nlargest(6).index
log("site_contribution_to_eta2_on_site_PC", {"PC": tss_pc, **contrib.sort_values(ascending=False).head(8).round(3).to_dict()})
fig = plt.figure(figsize=(14, 6.2))
gs_ = fig.add_gridspec(2, 4, width_ratios=[1.6, 1, 1, 1])
ax = fig.add_subplot(gs_[:, 0])
scatter_groups(ax, xy, lab.histology, ["IDC", "ILC", "Other"], HCOL, title="Histology", bg_label="—")
ax.legend(fontsize=8, markerscale=1.5, loc="upper left", bbox_to_anchor=(0, -0.02), ncol=2)
ax.set_xlabel("PC1"); ax.set_ylabel("PC2")
for i, site in enumerate(top_sites):
    a = fig.add_subplot(gs_[i // 3, 1 + i % 3])
    m = (site_f == site).values
    a.scatter(S.PC1[~m], S[f"PC{tss_pc}"][~m], s=3, c="#e1e0d9", lw=0, rasterized=True)
    a.scatter(S.PC1[m], S[f"PC{tss_pc}"][m], s=6, c="#6250d6", lw=0, rasterized=True)
    a.axhline(S[f"PC{tss_pc}"][m].mean(), color="#6250d6", lw=1, ls="--")
    a.set_title(f"site {site} (n={m.sum()}): {contrib[site]:.3f} of η²", fontsize=9); a.set_xticks([]); a.set_yticks([])
    if i % 3 == 0:
        a.set_ylabel(f"PC{tss_pc}")
    if i // 3 == 1:
        a.set_xlabel("PC1")
fig.suptitle(f"Histology on PC1/PC2; the six sites that contribute most to the site effect on PC{tss_pc} (η² {tss_max:.2f}; dashed = site mean)",
             fontsize=11, fontweight="bold", x=0.01, ha="left")
save(fig, "figS1b_pca_histology_site")

# ---- tumour + normal PCA (tumours only in the background/group layers)
TN = pd.concat([X, N[X.columns]])
TNv = z(TN[Xv.columns])
pt = PCA(n_components=2, random_state=SEED).fit(TNv)
st = pt.transform(TNv)
st[:, 0] *= np.sign(pt.components_[0, list(TNv.columns).index("ESR1")])
istum = (TN.index.str[13:15] == "01")
fig, ax = plt.subplots(figsize=(6.4, 4.8))
scatter_groups(ax, st[istum], lab.clinical_group.reindex(TN.index[istum]), GROUP_ORDER, GROUP_COL, s=8)
ax.scatter(st[~istum, 0], st[~istum, 1], s=16, c=GROUP_COL["Normal"], marker="D", lw=0.3, edgecolors="white",
           label=f"Adjacent normal (n={(~istum).sum()})")
ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.0, 1.0), markerscale=1.4)
ax.set_xlabel(f"PC1 ({pt.explained_variance_ratio_[0]*100:.1f}%)"); ax.set_ylabel(f"PC2 ({pt.explained_variance_ratio_[1]*100:.1f}%)")
ax.set_title("Tumours and adjacent normal tissue")
save(fig, "figS3_pca_tumour_normal")

# ---- UMAP
import umap
island = {}
for nn in UMAP_NN:
    U = umap.UMAP(n_neighbors=nn, min_dist=UMAP_MIN_DIST, random_state=SEED).fit_transform(S.iloc[:, :UMAP_PCS].values)
    Ud = pd.DataFrame(U, index=S.index, columns=["UMAP1", "UMAP2"])
    db = DBSCAN(eps=0.7, min_samples=10).fit_predict(U)
    basal = (lab.PAM50_pancan.reindex(S.index) == "Basal").values
    cand = pd.Series(db[basal & (db >= 0)])
    isl = (db == cand.value_counts().index[0]) if len(cand) else np.zeros(len(db), bool)
    Ud["island"] = isl
    Ud.to_csv(OUT / f"umap_nn{nn}.tsv", sep="\t")
    cg = lab.clinical_group.reindex(S.index).values
    ntn = int((cg == "TNBC").sum())
    a_ = int(((cg == "TNBC") & isl).sum())
    island[nn] = {"n_island": int(isl.sum()), "TNBC_in_island": f"{a_}/{ntn}", "HR+/HER2- in island": int(((cg == "HR+/HER2-") & isl).sum()),
                  "HER2+ in island": int(((cg == "HER2+") & isl).sum()), "Incomplete in island": int(((cg == "Incomplete") & isl).sum()),
                  "PAM50 Basal in island": f"{int((basal & isl).sum())}/{int(basal.sum())}"}
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8))
    scatter_groups(axes[0], U, cg, GROUP_ORDER, GROUP_COL, title="Clinical receptor group (IHC/FISH)")
    scatter_groups(axes[1], U, lab.PAM50_pancan.reindex(S.index), ["LumA", "LumB", "Her2", "Basal"], PAM50_COL,
                   title="PAM50 subtype (comparator)", bg_label="Normal-like")
    for a in axes:
        a.set_xlabel("UMAP1"); a.set_ylabel("UMAP2")
        a.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0), fontsize=8, markerscale=1.6, handletextpad=0.2)
    fig.suptitle(f"UMAP (n_neighbors={nn}): a separate island of {isl.sum()} tumours holds {a_}/{ntn} TNBC ({a_/ntn:.0%}), "
                 f"plus {island[nn]['HR+/HER2- in island']} HR+/HER2−, {island[nn]['HER2+ in island']} HER2+ "
                 f"and {island[nn]['Incomplete in island']} with incomplete labels",
                 fontsize=11, fontweight="bold", x=0.01, ha="left")
    fig.tight_layout()
    save(fig, f"fig3_umap_nn{nn}" if nn == 15 else f"figS4_umap_nn{nn}")
log("umap_island", island)

# ---- marker violins + tests (supervised set)
L = lab[sup]
rows = []
for g in MARKERS:
    for grp_name in ["TNBC", "HER2+"]:
        a = X.loc[L.index[L.clinical_group == grp_name], g].values
        b = X.loc[L.index[L.clinical_group == "HR+/HER2-"], g].values
        u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        rows.append((g, f"{grp_name} vs HR+/HER2-", len(a), len(b), float(np.median(a)), float(np.median(b)),
                     2 * u / (len(a) * len(b)) - 1, float(p)))
mt = pd.DataFrame(rows, columns=["gene", "contrast", "n1", "n2", "median1", "median2", "rank_biserial_r", "p"])
mt["log2_median_diff"] = mt.median1 - mt.median2
mt["q_BH"] = multipletests(mt.p, method="fdr_bh")[1]
assert (mt.p > 0).all(), "p-value underflow"
mt.to_csv(OUT / "marker_tests.tsv", sep="\t", index=False)
nG = L.clinical_group.value_counts()

d = X.loc[L.index, MARKERS].join(L.clinical_group)
d = d[d.clinical_group.isin(GROUP_ORDER)].melt(id_vars="clinical_group", var_name="gene", value_name="log2")
fig, axes = plt.subplots(2, 4, figsize=(13, 6.4))
for ax, g in zip(axes.flat, MARKERS):
    sns.violinplot(data=d[d.gene == g], x="clinical_group", y="log2", order=GROUP_ORDER, hue="clinical_group", palette=GROUP_COL,
                   legend=False, inner="quartile", cut=0, linewidth=0.8, ax=ax, density_norm="width", saturation=1)
    ax.set_title(g, fontstyle="italic"); ax.set_xlabel(""); ax.set_ylabel("log2 expression" if g in ("ESR1", "FOXA1") else "")
    ax.set_xticks(range(3)); ax.set_xticklabels([f"HR+/HER2−\nn={nG['HR+/HER2-']}", f"HER2+\nn={nG['HER2+']}", f"TNBC\nn={nG['TNBC']}"], fontsize=8)
    con = "HER2+ vs HR+/HER2-" if g == "ERBB2" else "TNBC vs HR+/HER2-"
    r = mt[(mt.gene == g) & (mt.contrast == con)].iloc[0]
    lo, hi = ax.get_ylim(); ax.set_ylim(lo, hi + (hi - lo) * 0.22)
    ax.text(0.98, 0.98, f"{con.split(' vs')[0]} vs HR+: r = {r.rank_biserial_r:+.2f}, q = {r.q_BH:.0e}", transform=ax.transAxes,
            ha="right", va="top", fontsize=7.5, color=INK2)
fig.suptitle("Marker genes by clinical group (Mann-Whitney; r = rank-biserial, + = higher than HR+/HER2−; BH-FDR over 16 tests)",
             fontsize=11, fontweight="bold", x=0.01, ha="left")
fig.tight_layout()
save(fig, "fig4_marker_violins")

# ---- mRNA vs clinical status: AUC with bootstrap CI (primary), Spearman with its ceiling
conc = {}
for gene, col in [("ESR1", "ER"), ("PGR", "PR"), ("ERBB2", "HER2")]:
    m = L[col].notna()
    y = (L.loc[m, col] == "Positive").astype(int).values
    x = X.loc[L.index[m], gene].values
    boots = []
    for _ in range(N_BOOT):
        i = rng.integers(0, len(y), len(y))
        if y[i].min() != y[i].max():
            boots.append(roc_auc_score(y[i], x[i]))
    p_ = y.mean()
    conc[f"{gene}_vs_{col}"] = {"n": int(m.sum()), "auc": round(float(roc_auc_score(y, x)), 3),
                                "auc_95ci": [round(float(v), 3) for v in np.percentile(boots, [2.5, 97.5])],
                                "spearman": round(float(stats.spearmanr(x, y)[0]), 3),
                                "spearman_ceiling_at_this_prevalence": round(float(np.sqrt(3 * p_ * (1 - p_))), 3)}
log("mrna_vs_clinical_status", conc)
m = L.HER2_v0.notna()
log("ERBB2_auc_with_HER2_v0_rule", round(float(roc_auc_score((L.HER2_v0[m] == "Positive").astype(int), X.loc[L.index[m], "ERBB2"])), 3))
by_src = {}
for src_, d_ in L[L.HER2.notna()].groupby("HER2_source"):
    y_ = (d_.HER2 == "Positive").astype(int)
    if 0 < y_.sum() < len(y_):
        by_src[src_] = {"n": len(d_), "positives": int(y_.sum()), "auc": round(float(roc_auc_score(y_, X.loc[d_.index, "ERBB2"])), 3)}
log("ERBB2_auc_by_HER2_source", by_src)

# ---- PAM50 heatmap ordered by PAM50, annotated by clinical group
p50 = pd.read_csv(PROC / "pam50_centroids_hgnc.tsv", sep="\t", index_col=0).index
order = lab.assign(p=pd.Categorical(lab.PAM50_pancan, PAM50_ORDER),
                   g=pd.Categorical(lab.clinical_group, GROUP_ORDER + ["Incomplete"])).sort_values(["p", "g"]).index
H = z(X.loc[order, p50]).clip(-3, 3).T
col_colors = pd.DataFrame({"PAM50": lab.loc[order, "PAM50_pancan"].map(PAM50_COL),
                           "Clinical": lab.loc[order, "clinical_group"].map(GROUP_COL)})
m = lab.clinical_group.isin(GROUP_ORDER) & lab.PAM50_pancan.notna()
tn, ba = (lab.clinical_group[m] == "TNBC"), (lab.PAM50_pancan[m] == "Basal")
a, b, c_, dd_ = int((tn & ba).sum()), int((tn & ~ba).sum()), int((~tn & ba).sum()), int((~tn & ~ba).sum())
cg = sns.clustermap(H, col_cluster=False, row_cluster=True, method="average", metric="correlation", cmap="RdBu_r",
                    center=0, vmin=-3, vmax=3, col_colors=col_colors, xticklabels=False, yticklabels=True,
                    figsize=(11, 8.5), dendrogram_ratio=(0.08, 0.02), colors_ratio=0.025,
                    cbar_pos=(1.02, 0.35, 0.012, 0.15), cbar_kws={"label": "z-score"})
cg.ax_heatmap.tick_params(axis="y", labelsize=6.5); cg.ax_heatmap.set_ylabel("")
cg.ax_heatmap.set_xlabel(f"{len(order):,} tumours, ordered by PAM50 subtype then clinical group")
from matplotlib.patches import Patch
hand = [Patch(color=PAM50_COL[k], label=k) for k in PAM50_ORDER] + [Patch(color="none", label="")] + \
       [Patch(color=GROUP_COL[k], label=k) for k in GROUP_ORDER + ["Incomplete"]]
cg.ax_heatmap.legend(handles=hand, ncol=10, loc="lower center", bbox_to_anchor=(0.5, 1.06), fontsize=8, frameon=False)
cg.figure.suptitle(f"PAM50 genes: {a}/{a+b} ({a/(a+b):.0%}) of clinical TNBC fall in the PAM50 Basal block", y=1.04,
                   fontsize=11, fontweight="bold")
save(cg.figure, "fig5_pam50_heatmap")

# ---- clinical group vs PAM50 concordance
ct = pd.crosstab(lab.clinical_group, lab.PAM50_pancan).reindex(index=GROUP_ORDER + ["Incomplete"], columns=PAM50_ORDER).fillna(0).astype(int)
ct.to_csv(OUT / "clinical_group_vs_pam50.tsv", sep="\t")
OR, _ = stats.fisher_exact([[a, b], [c_, dd_]])
se = np.sqrt(1 / a + 1 / b + 1 / c_ + 1 / dd_)
log("TNBC_vs_Basal", {"TNBC_called_Basal": f"{a}/{a+b} ({a/(a+b):.1%})", "Basal_that_are_TNBC": f"{a}/{a+c_} ({a/(a+c_):.1%})",
                      "odds_ratio": round(float(OR), 1), "OR_95CI_woolf": [round(float(np.exp(np.log(OR) - 1.96 * se)), 1), round(float(np.exp(np.log(OR) + 1.96 * se)), 1)],
                      "kappa": round(float(cohen_kappa_score(tn, ba)), 3),
                      "TCGA2012_reference": "65/76 (86%) TNBC basal-like; 65/81 (80%) basal-like TNBC (Nature 2012, doi 10.1038/nature11412; "
                                            "TCGA 2012 analysed a subset of these same tumours, so this is not an independent cohort)"})
h2 = lab.clinical_group[m] == "HER2+"; he = lab.PAM50_pancan[m] == "Her2"
her2e_share = (h2 & he).sum() / h2.sum()


def share(hmask, pcol):
    k_ = hmask & pcol.notna()
    return f"{int((k_ & (pcol == 'Her2')).sum())}/{int(k_.sum())}"


h12 = lab.HER2_Final_Status_nature2012 == "Positive"
hv1 = lab.clinical_group == "HER2+"
like = {"TCGA2012_HER2_x_TCGA2012_PAM50": share(h12, lab.PAM50_nature2012),
        "TCGA2012_HER2_x_PanCan_PAM50": share(h12, lab.PAM50_pancan),
        "v1_HER2pos_with_a_2012_PAM50_call_x_PanCan": share(hv1 & lab.PAM50_nature2012.notna(), lab.PAM50_pancan),
        "v1_HER2pos_with_a_2012_PAM50_call_x_TCGA2012_PAM50": share(hv1, lab.PAM50_nature2012),
        "v1_HER2pos_by_HER2_source": {k_: share(hv1 & (lab.HER2_source == k_), lab.PAM50_pancan)
                                      for k_ in lab.HER2_source[hv1].dropna().unique()}}
t12 = like["TCGA2012_HER2_x_TCGA2012_PAM50"].split("/")
like12 = int(t12[0]) / int(t12[1])
log("HER2pos_vs_HER2E", {"HER2+_called_Her2E": f"{int((h2 & he).sum())}/{int(h2.sum())} ({her2e_share:.1%})",
                         "Her2E_that_are_HER2+": f"{int((h2 & he).sum())}/{int(he.sum())}",
                         "HER2+_PAM50_distribution": lab.PAM50_pancan[m][h2].value_counts().to_dict(),
                         "like_for_like": like,
                         "TCGA2012_reference": "'only ~50% of clinically HER2-positive tumors fall into this HER2E-mRNA-subtype' (doi 10.1038/nature11412)"})
erp = (lab.ER == "Positive") & lab.PAM50_pancan.notna()
nl1 = erp & lab.PAM50_pancan.isin(["Basal", "Her2"])
nl2 = erp & lab.PAM50_pancan.isin(["Basal", "Her2", "Normal"])
log("ERpos_called_nonluminal", {"Basal_or_Her2E": f"{int(nl1.sum())}/{int(erp.sum())} ({nl1.sum()/erp.sum():.1%})",
                                "Basal_Her2E_or_Normal": f"{int(nl2.sum())}/{int(erp.sum())} ({nl2.sum()/erp.sum():.1%})",
                                "Nielsen2010_reference": "786 ER+ tamoxifen-treated tumours, qRT-PCR PAM50: HER2E 8.1% + Basal 0.6% = 8.7%; "
                                                         "with Normal-like 2.0% = 10.7% (doi 10.1158/1078-0432.CCR-10-1282)"})

fig, ax = plt.subplots(figsize=(6.8, 3.6))
sns.heatmap(ct.div(ct.sum(1), axis=0), annot=ct.values, fmt="d", cmap="Blues", vmin=0, vmax=1, linewidths=2,
            linecolor="white", cbar_kws={"label": "share of row"}, ax=ax)
ax.set_xlabel("PAM50 subtype (PanCanAtlas)"); ax.set_ylabel(""); ax.tick_params(axis="y", rotation=0)
ax.set_title(f"{a/(a+b):.0%} of TNBC are PAM50 Basal; {her2e_share:.0%} of HER2+ are Her2-enriched "
             f"({like12:.0%} with the TCGA 2012 HER2 labels and PAM50 calls)", fontsize=9.5)
save(fig, "fig6_clinical_vs_pam50")

# ---- mRNA-protein preview (RPPA)
R = pd.read_parquet(PROC / "tcga_rppa.parquet")
abm = pd.read_csv(PROC / "rppa_antibody_map.tsv", sep="\t")
abm = abm[abm.total_protein & abm.gene.isin(X.columns)]
common = R.index.intersection(X.index)
mp = pd.DataFrame([(r.antibody, r.gene, r.isoform_ambiguous, stats.spearmanr(X.loc[common, r.gene], R.loc[common, r.antibody])[0])
                   for r in abm.itertuples()], columns=["antibody", "gene", "isoform_ambiguous", "spearman"]).sort_values("spearman", ascending=False)
mp.to_csv(OUT / "mrna_protein_spearman.tsv", sep="\t", index=False)
mpu = mp[~mp.isoform_ambiguous]
log("mrna_protein_spearman", {"n_tumours": len(common), "n_pairs_unambiguous": len(mpu), "median": round(float(mpu.spearman.median()), 3),
                              "selected": mp.set_index("antibody").spearman.reindex(["ERALPHA", "PR", "HER2", "GATA3", "AR", "BCL2", "CYCLINB1", "ECADHERIN", "EGFR"]).round(2).to_dict()})
fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4), gridspec_kw={"width_ratios": [1.25, 1]})
top = mpu.head(25)
axes[0].barh(top.antibody[::-1] + " (" + top.gene[::-1] + ")", top.spearman[::-1], color="#2a78d6", height=0.7)
axes[0].axvline(mpu.spearman.median(), color=MUTED, lw=1, ls="--")
axes[0].text(mpu.spearman.median() + 0.01, -0.9, f"median of {len(mpu)} = {mpu.spearman.median():.2f}", fontsize=7, color=INK2)
axes[0].set_xlabel("Spearman ρ, mRNA vs protein (same tumour)"); axes[0].tick_params(axis="y", labelsize=7)
axes[0].set_title(f"Top 25 of {len(mpu)} total-protein antibodies (n={len(common)} tumours)")
cgp = lab.loc[common, "clinical_group"]
for g in ["Incomplete"] + GROUP_ORDER:
    mm = cgp == g
    axes[1].scatter(X.loc[common, "ESR1"][mm], R.loc[common, "ERALPHA"][mm], s=8, c=GROUP_COL[g], lw=0, label=g, alpha=0.9)
axes[1].set_xlabel("ESR1 mRNA (log2)"); axes[1].set_ylabel("ERα protein (RPPA)")
axes[1].set_title(f"ESR1 mRNA vs ERα protein, ρ = {mp.set_index('antibody').spearman['ERALPHA']:.2f}")
axes[1].legend(fontsize=8, markerscale=1.5, loc="upper left")
fig.tight_layout()
save(fig, "figS5_mrna_protein")

# ---- STAR cross-check of the PCA
C = pd.read_parquet(PROC / "tcga_star_counts_tumour_filtered.parquet")
ann = pd.read_csv(PROC / "gencode_v36_annotation.tsv", sep="\t", index_col=0)
C = C.loc[C.index.intersection(X.index)]
Lg = np.log2(C.div(C.sum(1), axis=0) * 1e6 + 1)
Lg.columns = ann.hgnc.reindex(Lg.columns).values
Lg = Lg.T.groupby(level=0).max().T
gg = [g for g in Xv.columns if g in Lg.columns]
ps = PCA(n_components=5, random_state=SEED).fit_transform(z(Lg[gg]))
cr = [round(abs(float(np.corrcoef(ps[:, i], S.loc[Lg.index, f"PC{i+1}"])[0, 1])), 3) for i in range(5)]
log("star_vs_hiseq_pc_score_abs_correlation_PC1_5", {"n_tumours": len(Lg), "genes": len(gg), "r": cr})

write_params(OUT, PARAMS)
json.dump(rep, open(OUT / "explore_report.json", "w"), indent=2, default=str)
print("done")
