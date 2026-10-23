"""Step 3b: clustering figures, drawn from the tables written by 03_cluster.py (runs in seconds).

Outputs: figures/fig7-9*, figures/figS7-S8*, results/cluster/netanely_comparison.tsv, ari_ceiling_k*.json,
         cluster_purity_site_k*.tsv, cluster_confounder_check.json
Run from the project root:  python scripts/03b_cluster_figures.py
"""
import json, sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from sklearn.metrics import adjusted_rand_score as ARI
from src.common import PROC, RES, GROUP_ORDER, PAM50_ORDER, CLUSTER_COL, METHOD_COL, INK, INK2, MUTED, style, save, drop_excluded, write_params

style()
OUT = RES / "cluster"
rep = json.load(open(OUT / "cluster_report.json"))
k_star = rep["k_star_lowest_PAC"]
K_REPORT = sorted({k_star, 4, 5})
scan = pd.read_csv(OUT / "k_selection.tsv", sep="\t", index_col=0)
ag = pd.read_csv(OUT / "cluster_agreement.tsv", sep="\t")
clus = pd.read_csv(PROC / "tcga_clusters.tsv", sep="\t", index_col=0)
lab = drop_excluded(pd.read_csv(PROC / "tcga_labels.tsv", sep="\t", index_col=0)).loc[clus.index]
cdf = pd.read_csv(OUT / "consensus_cdf_quantiles.tsv", sep="\t", index_col=0)
sens = rep["sensitivity_ARI_vs_main"]
magree = rep["method_agreement_ARI"]

# ---- k selection
fig, axes = plt.subplots(1, 4, figsize=(15, 3.3))
for ax, col, t in zip(axes, ["PAC", "silhouette", "inertia"], ["PAC (lower = more stable)", "mean silhouette", "within-cluster SS (elbow)"]):
    ax.plot(scan.index, scan[col], "-o", color="#2a78d6", ms=6, lw=2)
    ax.set_xlabel("k"); ax.set_title(t); ax.set_xticks(scan.index)
    ax.axvline(k_star, color=MUTED, ls="--", lw=1)
axes[0].text(k_star + 0.15, scan.PAC.max() * 0.95, f"k* = {k_star}", fontsize=8, color=INK2, va="top")
for i, k in enumerate(scan.index):
    axes[3].plot(cdf[str(k)].values, cdf.index.values, lw=1.6, color=CLUSTER_COL[i], label=f"k={k}")
axes[3].axvspan(0.1, 0.9, color="#f0efec", zorder=0)
axes[3].set_xlabel("consensus index"); axes[3].set_ylabel("CDF"); axes[3].set_title("consensus CDF (shaded = PAC window)")
axes[3].legend(fontsize=7, ncol=2, loc="lower right", frameon=True, facecolor="white", edgecolor="none")
sil_txt = "silhouette agrees" if int(scan.silhouette.idxmax()) == k_star else f"silhouette peaks at k = {int(scan.silhouette.idxmax())}"
fig.suptitle(f"Choosing k: consensus clustering (100 × 80% resamples) is most stable at k = {k_star}; {sil_txt}",
             fontsize=11, fontweight="bold", x=0.01, ha="left")
fig.tight_layout()
save(fig, "fig7_k_selection")


# ---- contingency tables
def contingency(ax, col, ref, order, title):
    t = pd.crosstab(clus[col].map(lambda c: f"C{c}"), ref).reindex(columns=order).fillna(0).astype(int)
    sns.heatmap(t.div(t.sum(0), axis=1), annot=t.values, fmt="d", cmap="Blues", vmin=0, vmax=1, cbar=False,
                linewidths=2, linecolor="white", ax=ax, annot_kws={"size": 9})
    ax.set_title(title, fontsize=10); ax.set_xlabel(""); ax.set_ylabel(""); ax.tick_params(axis="y", rotation=0)


for k in K_REPORT:
    a1 = ag[(ag.clustering == f"kmeans_k{k}") & (ag.reference == "Clinical group")].iloc[0]
    a2 = ag[(ag.clustering == f"kmeans_k{k}") & (ag.reference == "PAM50 (no Normal-like)")].iloc[0]
    fig, axes = plt.subplots(1, 2, figsize=(10, 0.55 * k + 2.4), gridspec_kw={"width_ratios": [3, 4]})
    contingency(axes[0], f"kmeans_k{k}", lab.clinical_group, GROUP_ORDER,
                f"vs clinical group (n={a1.n}): ARI {a1.ARI:.2f} [{a1.ARI_lo:.2f}–{a1.ARI_hi:.2f}]")
    contingency(axes[1], f"kmeans_k{k}", lab.PAM50_pancan, ["LumA", "LumB", "Her2", "Basal"],
                f"vs PAM50 excl. Normal-like (n={a2.n}): ARI {a2.ARI:.2f} [{a2.ARI_lo:.2f}–{a2.ARI_hi:.2f}]")
    stab = (f"PAC {scan.PAC[k]:.3f}; vs Ward ARI {magree[f'k{k}']['kmeans_vs_ward']:.2f}; "
            f"vs 5,000-gene run ARI {sens[f'top5000_excl_PAM50_k{k}']:.2f}")
    fig.suptitle(f"k-means clusters (k={k}) against labels; cell = count, shade = share of column", fontsize=11,
                 fontweight="bold", x=0.01, ha="left")
    fig.text(0.01, -0.02, f"Stability: {stab}", fontsize=8, color=INK2)
    fig.tight_layout()
    save(fig, f"fig8_contingency_k{k}" if k == k_star else f"figS7_contingency_k{k}")

# ---- ARI forest: methods at k*, k-means at k=5, and label-vs-label reference points
refs = ["Clinical group", "PAM50 (no Normal-like)", "TCGA 2012 SigClust", "ER status (IHC)"]
meth = [f"kmeans_k{k_star}", f"ward_k{k_star}", f"avgcorr_k{k_star}"]
if magree[f"k{k_star}"]["kmeans_vs_consensus"] < 0.999:
    meth.insert(1, f"consensus_k{k_star}")
meth += ["kmeans_k5"]
names = {f"kmeans_k{k_star}": f"k-means, k={k_star}", f"consensus_k{k_star}": f"consensus, k={k_star}", f"ward_k{k_star}": f"Ward, k={k_star}",
         f"avgcorr_k{k_star}": f"average-linkage 1−r, k={k_star}", "kmeans_k5": "k-means, k=5"}
nref = ag[ag.clustering == f"kmeans_k{k_star}"].set_index("reference").n
fig, ax = plt.subplots(figsize=(10, 4.6))
for j, mth in enumerate(meth):
    s = ag[ag.clustering == mth].set_index("reference").reindex(refs)
    y = np.arange(len(refs)) + (j - (len(meth) - 1) / 2) * 0.14
    ax.errorbar(s.ARI, y, xerr=[s.ARI - s.ARI_lo, s.ARI_hi - s.ARI], fmt="o", ms=6, color=METHOD_COL.get(mth, CLUSTER_COL[j]), lw=1.5, label=names[mth])
lv = ag[ag.clustering == "PAM50 (no Normal-like)"].set_index("reference")
p4 = lab.PAM50_pancan.where(lab.PAM50_pancan.isin(["LumA", "LumB", "Her2", "Basal"]))
p2 = p4.map({"LumA": "luminal", "LumB": "luminal", "Her2": "Basal+Her2", "Basal": "Basal+Her2"})
for i, r in enumerate(refs):
    if r in lv.index:
        ax.plot(lv.loc[r, "ARI"], i, marker="D", ms=8, mfc="white", mec=INK, mew=1.5, ls="none",
                label="PAM50, 4 classes (not class-matched)" if i == 0 else None)
# ceiling: best ARI any k*-cluster partition built by merging the reference classes can reach
REFCOL = {"Clinical group": lab.clinical_group.where(lab.clinical_group.isin(GROUP_ORDER)),
          "PAM50 (no Normal-like)": lab.PAM50_pancan.where(lab.PAM50_pancan.isin(["LumA", "LumB", "Her2", "Basal"])),
          "TCGA 2012 SigClust": lab.SigClust_2012, "ER status (IHC)": lab.ER}


def merges(classes, k):
    """All ways to split a list of classes into exactly k non-empty blocks (set partitions)."""
    if k == 1:
        yield [classes]; return
    if len(classes) == k:
        yield [[c] for c in classes]; return
    if len(classes) < k:
        return
    first, rest = classes[0], classes[1:]
    for part in merges(rest, k - 1):
        yield [[first]] + part
    for part in merges(rest, k):
        for i in range(len(part)):
            yield part[:i] + [[first] + part[i]] + part[i + 1:]


ceiling = {}
for r, v in REFCOL.items():
    v = v.dropna().astype(str)
    cls = sorted(v.unique())
    best = np.nan
    if len(cls) >= k_star and (len(cls) <= 8 or k_star <= 3):   # keeps the number of set partitions small
        best = max(ARI(v, v.map({c: i for i, b in enumerate(part) for c in b})) for part in merges(cls, k_star))
    ceiling[r] = {"n_classes": len(cls), f"best_ARI_with_{k_star}_merged_groups": None if np.isnan(best) else round(float(best), 3)}
json.dump(ceiling, open(OUT / f"ari_ceiling_k{k_star}.json", "w"), indent=2)
# class-matched comparator: PAM50 collapsed to 2 groups (luminal vs Basal+Her2), scored on the same tumours
pam2 = {}
for i, r in enumerate(refs):
    if r == "PAM50 (no Normal-like)":
        continue
    v = REFCOL[r]
    mm = v.notna() & p2.notna()
    pam2[r] = round(float(ARI(v[mm].astype(str), p2[mm])), 3)
    ax.plot(pam2[r], i, marker="D", ms=8, mfc=MUTED, mec=INK, mew=1.0, ls="none",
            label="PAM50 as 2 groups (luminal vs Basal+Her2)" if len(pam2) == 1 else None)
json.dump(pam2, open(OUT / "ari_pam50_two_groups_vs_refs.json", "w"), indent=2)
ystar = [(j - (len(meth) - 1) / 2) * 0.14 for j, mth in enumerate(meth) if mth.endswith(f"_k{k_star}")]
for i, r in enumerate(refs):
    b = ceiling[r][f"best_ARI_with_{k_star}_merged_groups"]
    if b is not None:
        ax.plot([b, b], [i + min(ystar) - 0.06, i + max(ystar) + 0.06], color=INK2, lw=2.2, solid_capstyle="butt",
                label=f"best ARI possible with {k_star} clusters (k = {k_star} rows only)" if i == 0 else None)
ax.set_yticks(range(len(refs)))
ax.set_yticklabels([f"{r}\n(n={nref[r]:,}; {ceiling[r]['n_classes']} classes)" for r in refs])
ax.invert_yaxis()
ax.set_xlabel("Adjusted Rand index (95% bootstrap CI)"); ax.set_xlim(0, 1)
ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1.0))
ax.set_title("Agreement of label-free clusters with each reference")
save(fig, "fig9_ari_forest")

# ---- clusters on the UMAP (direct labels at cluster medians)
U = pd.read_csv(RES / "explore/umap_nn15.tsv", sep="\t", index_col=0).loc[clus.index]
for k in sorted({k_star, 5}):
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    for c in range(k):
        m = (clus[f"kmeans_k{k}"] == c).values
        ax.scatter(U.UMAP1[m], U.UMAP2[m], s=8, c=CLUSTER_COL[c], lw=0, label=f"C{c} (n={m.sum()})", rasterized=True)
        ax.text(np.median(U.UMAP1[m]), np.median(U.UMAP2[m]), f"C{c}", fontsize=10, fontweight="bold", ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.85))
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
    ax.legend(fontsize=8, markerscale=1.5, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    ax.set_title(f"k-means clusters (k={k}) on the UMAP")
    save(fig, f"figS8_umap_clusters_k{k}")

# ---- are clusters driven by tumour purity or by tissue source site? (simple marker-gene proxies)
from scipy.stats import chi2_contingency
PROXY = {"adipocyte": ["ADIPOQ", "PLIN1", "LEP", "CIDEA", "FABP4", "LPL", "CD36"],
         "stromal": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FAP", "PDGFRB"],
         "immune": ["PTPRC", "CD3E", "CD8A", "MS4A1", "CD68", "CD2", "LCK"],
         "proliferation": ["MKI67", "TOP2A", "CCNB1", "BIRC5", "MYBL2", "AURKA"]}
X = pd.read_parquet(PROC / "tcga_hiseq_tumour_log2.parquet").loc[clus.index].astype("float64")
Zx = (X - X.mean()) / X.std(ddof=0)
score = pd.DataFrame({k_: Zx[[g for g in gl if g in Zx.columns]].mean(1) for k_, gl in PROXY.items()})
site = lab.tss.where(lab.tss.map(lab.tss.value_counts()) >= 10, "other")


def cramers_v(a, b):
    t = pd.crosstab(a, b)
    chi2, p, _, _ = chi2_contingency(t)
    return float(np.sqrt(chi2 / (t.values.sum() * (min(t.shape) - 1)))), float(p)


conf = {}
for k in K_REPORT:
    c = clus[f"kmeans_k{k}"]
    tab = score.groupby(c).mean().round(2)
    tab.index = [f"C{i}" for i in tab.index]
    top = pd.crosstab(c, site).drop(columns="other", errors="ignore")
    tab["top_site"] = [f"{top.loc[i].idxmax()} ({top.loc[i].max() / (c == i).sum():.0%})" for i in top.index]
    tab.to_csv(OUT / f"cluster_purity_site_k{k}.tsv", sep="\t")
    v_site, p_site = cramers_v(c, site)
    v_p50, _ = cramers_v(c[lab.PAM50_pancan.notna()], lab.PAM50_pancan.dropna())
    conf[f"k{k}"] = {"cramers_V_vs_site": round(v_site, 3), "p_site": float(f"{p_site:.2e}"), "cramers_V_vs_PAM50": round(v_p50, 3),
                     "proxy_means_z": tab.drop(columns="top_site").to_dict(orient="index")}
json.dump(conf, open(OUT / "cluster_confounder_check.json", "w"), indent=2)
print(json.dumps(conf, indent=1))

# ---- Netanely 2016 comparison table (k = 5)
prof = pd.read_csv(OUT / "cluster_profile_kmeans_k5.tsv", sep="\t")


def netanely(r):
    """Heuristic by majority PAM50 class, not Netanely's published assignment (they also clustered 113 normals)."""
    if r.top_PAM50 == "Basal" and r.top_PAM50_share >= 0.8:
        return "basal"
    if r.top_PAM50 == "Her2":
        return "HER2-enriched"
    if r.top_PAM50 == "LumA" and r.top_PAM50_share >= 0.75:
        return "homogeneous LumA"
    return "mixed LumA/LumB"


prox5 = pd.DataFrame(conf["k5"]["proxy_means_z"]).T if "k5" in conf else None
prof["closest_Netanely2016_cluster"] = prof.apply(netanely, axis=1)
if prox5 is not None:
    prof = prof.join(prox5.add_prefix("proxy_z_"), on="cluster")
    prof["note"] = np.where((prof.proxy_z_adipocyte > 0.5) & (prof.proxy_z_proliferation < -0.5),
                            "high adipocyte/stromal/immune proxies, low proliferation: lower tumour content is possible", "")
cols = ["cluster", "n", "PAM50", "HR+/HER2-", "HER2+", "TNBC", "ILC_share", "ER_neg_share", "closest_Netanely2016_cluster"]
cols += [c for c in prof.columns if c.startswith("proxy_z_")] + (["note"] if "note" in prof else []) + ["top_genes_up"]
prof[cols].to_csv(OUT / "netanely_comparison.tsv", sep="\t", index=False)
print(prof[["cluster", "n", "PAM50", "ILC_share", "closest_Netanely2016_cluster"]].to_string(index=False))

# ---- per-cluster stability: best-match Jaccard of each k-means cluster in every alternative run at the same k
for k in K_REPORT:
    main = clus[f"kmeans_k{k}"]
    alts = [c for c in clus.columns if c != f"kmeans_k{k}" and c.endswith(f"_k{k}")]
    jac = pd.DataFrame(index=[f"C{c}" for c in range(k)], columns=alts, dtype=float)
    for alt in alts:
        for c in range(k):
            A = main == c
            jac.loc[f"C{c}", alt] = max(((A & (clus[alt] == d)).sum() / (A | (clus[alt] == d)).sum()) for d in clus[alt].unique())
    jac.columns = [c.replace(f"_k{k}", "").replace("sens_", "").replace("_kmeans", "") for c in jac.columns]
    jac["min"], jac["median"] = jac.min(1), jac.median(1)
    if k == 5:
        jac.insert(0, "closest_Netanely2016_cluster", prof.set_index("cluster").closest_Netanely2016_cluster)
    jac.round(2).to_csv(OUT / f"cluster_stability_jaccard_k{k}.tsv", sep="\t")
    print(f"k={k} best-match Jaccard per cluster:\n", jac.round(2).to_string())
write_params(OUT, {
    "purity_proxies": {k_: v for k_, v in PROXY.items()},
    "purity_proxy_score": "mean z-score (across the analysed tumours) of the proxy genes present",
    "site_pooling": "tissue source sites with < 10 tumours pooled as 'other'; Cramer's V from the chi-squared test of cluster x site",
    "ari_ceiling": f"maximum ARI over all partitions of the reference classes into k* = {k_star} blocks (exhaustive set-partition enumeration)",
    "pam50_two_groups": "PAM50 collapsed to luminal (LumA, LumB) vs Basal+Her2, Normal-like left out; scored on the tumours each reference labels",
    "cluster_stability": "for each k-means cluster, the best-match Jaccard index (|A and B| / |A or B|) with any cluster of each "
                         "alternative run at the same k (Ward, average-linkage, consensus, and the gene-set sensitivity runs); "
                         "Jaccard < 0.5 is commonly read as a cluster that dissolves (Hennig 2007)",
    "netanely_mapping": "heuristic by majority PAM50 class: basal if top class Basal with share >= 0.8; HER2-enriched if top class Her2; "
                        "homogeneous LumA if top class LumA with share >= 0.75; else mixed LumA/LumB",
}, name="params_figures.json")
print("done")
