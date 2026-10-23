"""Step 4: METABRIC quality control and a first look at whether the TCGA structure replicates.

Uses METABRIC labels from other assays only (ER by IHC, HER2 by SNP6 copy number); the array-derived
PAM50 + claudin-low call is shown for comparison only.
Outputs: figures/fig10_*, figures/figS9_*, results/metabric/*, data/processed/metabric_clusters.tsv
Run from the project root:  python scripts/04_metabric_explore.py   (after 03_cluster.py)
"""
import json, sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score as ARI, roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.common import PROC, RES, SEED, INK2, style, save, write_params, drop_excluded

PARAMS = {
    "samples": "all 1,980 METABRIC tumours with expression (Illumina HT-12, log2 intensity as distributed by cBioPortal)",
    "genes": "the TCGA main clustering genes (top 2,000 non-PAM50, sex genes removed) that are present in METABRIC",
    "scaling": "z-score per gene within METABRIC (each cohort standardised separately). This removes platform offsets but "
               "assumes the two cohorts have a broadly similar subtype mix, because a cohort with more ER-negative tumours "
               "would have its gene means shifted; the ER-IHC negative share is logged per cohort to check this "
               "(ER_negative_share_by_cohort)",
    "quantile_mapping_sensitivity": "each METABRIC gene mapped onto the TCGA distribution of the same gene (rank -> TCGA "
                                    "quantile), then z-scored with the TCGA mean and SD and assigned by the same nearest-"
                                    "centroid rule; reported per k as cluster shares, ER-negative shares and ARI vs the main "
                                    "assignment. Not used for the main result",
    "pc1_orientation": "PC1 sign chosen so that PC1 correlates positively with ESR1 expression",
    "batch_check": "eta-squared of PC1-10 on METABRIC COHORT (1-5, known processing batches), ER_IHC, HER2 gain, "
                   "CLAUDIN_SUBTYPE (NC dropped)",
    "centroid_projection": "each METABRIC tumour assigned to the nearest TCGA k-means centroid by Euclidean distance "
                           "in z-units over the shared genes (the k-means rule itself); checked by re-applying the "
                           "rule to TCGA, where it must reproduce the k-means labels",
    "centroid_projection_sensitivity": "Pearson-correlation nearest centroid, reported but not used",
    "ER_negative_share": "share of ER-IHC negative among tumours with a known ER-IHC result",
}
OUT = RES / "metabric"
OUT.mkdir(parents=True, exist_ok=True)
style()
rep = {}


def log(k, v):
    rep[k] = v
    print(f"{k}: {v}", flush=True)


M = pd.read_parquet(PROC / "metabric_log2.parquet")
ml = pd.read_csv(PROC / "metabric_labels.tsv", sep="\t", index_col=0).loc[M.index]
G = pd.read_csv(RES / "cluster/cluster_genes_main.txt", header=None)[0].tolist()
Gm = [g for g in G if g in M.columns]
log("main_cluster_genes_present_in_metabric", f"{len(Gm)}/{len(G)}")


def z(A):
    return (A - A.mean()) / A.std(ddof=0)


Z = z(M[Gm])

pca = PCA(n_components=10, random_state=SEED).fit(Z)
S = pd.DataFrame(pca.transform(Z), index=Z.index, columns=[f"PC{i}" for i in range(1, 11)])
if np.corrcoef(S.PC1, M["ESR1"])[0, 1] < 0:
    S["PC1"] *= -1
ev = pca.explained_variance_ratio_


def eta2(y, g):
    m = g.notna(); y, g = y[m], g[m]
    return float(sum(len(v) * (v.mean() - y.mean()) ** 2 for _, v in y.groupby(g)) / ((y - y.mean()) ** 2).sum())


fac = {"Cohort (batch)": ml.COHORT.astype("Int64").astype(str).where(ml.COHORT.notna()), "ER (IHC)": ml.ER_IHC,
       "HER2 gain (SNP6)": ml.HER2_gain, "PAM50+claudin-low": ml.CLAUDIN_SUBTYPE.where(ml.CLAUDIN_SUBTYPE != "NC")}
E = pd.DataFrame({k: [eta2(S[f"PC{i}"], v) for i in range(1, 11)] for k, v in fac.items()},
                 index=[f"PC{i}\n({ev[i-1]*100:.1f}%)" for i in range(1, 11)]).T
E.to_csv(OUT / "pc_variance_by_factor_eta2.tsv", sep="\t")
log("eta2_PC1_PC3", E.iloc[:, :3].round(3).to_dict())
coh = E.loc["Cohort (batch)"]
coh_max_pc = int(np.argmax(coh.values)) + 1
log("cohort_eta2_max", {"PC": coh_max_pc, "eta2": round(float(coh.max()), 3), "PC1_3_max": round(float(coh.iloc[:3].max()), 3)})
er_known = ml.ER_IHC.notna()
pc1_er_auc = float(roc_auc_score((ml.ER_IHC == "Positive")[er_known], S.PC1[er_known]))
log("PC1_AUC_for_ER_IHC", round(pc1_er_auc, 3))

# PCA coloured by ER-IHC / HER2-gain group
grp = ml.contrast_group.map({"ER+/HER2-notgained": "ER+/HER2 not gained", "HER2-gained": "HER2 gained",
                             "ER-/HER2-notgained": "ER−/HER2 not gained"})
col = {"ER+/HER2 not gained": "#2a78d6", "HER2 gained": "#eb6834", "ER−/HER2 not gained": "#1baf7a"}
fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4), gridspec_kw={"width_ratios": [1.1, 1]})
ax = axes[0]
bg = grp.isna().values
ax.scatter(S.PC1[bg], S.PC2[bg], s=5, c="#d9d8d2", lw=0, label=f"label missing (n={bg.sum()})")
for g, c in col.items():
    m = (grp == g).values
    ax.scatter(S.PC1[m], S.PC2[m], s=7, c=c, lw=0.2, edgecolors="white", label=f"{g} (n={m.sum()})", rasterized=True)
ax.set_xlabel(f"PC1 ({ev[0]*100:.1f}%)"); ax.set_ylabel(f"PC2 ({ev[1]*100:.1f}%)"); ax.set_xticks([]); ax.set_yticks([])
ax.legend(fontsize=8, markerscale=1.6, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.1))
ax.set_title(f"METABRIC PCA on {len(Gm):,} of the {len(G):,} TCGA clustering genes", fontsize=10)
sns.heatmap(E, annot=True, fmt=".2f", cmap="Blues", vmin=0, vmax=0.8, cbar=False, linewidths=2, linecolor="white",
            ax=axes[1], annot_kws={"size": 7})
axes[1].tick_params(axis="x", labelsize=7); axes[1].tick_params(axis="y", rotation=0)
axes[1].set_title("η² of each PC by factor", fontsize=10)
coh_txt = (f"processing cohort explains ≤{coh.iloc[:3].max():.2f} of PC1–3 "
           f"(most on PC{coh_max_pc}, η² {coh.max():.2f})")
fig.suptitle(f"Independent cohort (n={len(M):,}): PC1 is again the ER axis (AUC {pc1_er_auc:.2f} for IHC ER); {coh_txt}",
             fontsize=10.5, fontweight="bold", x=0.01, ha="left")
fig.tight_layout()
save(fig, "fig10_metabric_pca")

# ------------------------------------------------ project TCGA k-means centroids (Euclidean nearest centroid)
tl = drop_excluded(pd.read_csv(PROC / "tcga_labels.tsv", sep="\t", index_col=0))
tc = pd.read_csv(PROC / "tcga_clusters.tsv", sep="\t", index_col=0)
TX = pd.read_parquet(PROC / "tcga_hiseq_tumour_log2.parquet").loc[tc.index, G].astype("float64")
TZ = z(TX)


def nearest(Zs, cen):
    d = ((Zs.values[:, None, :] - cen.T.values[None, :, :]) ** 2).sum(-1)
    return pd.Series(d.argmin(1), index=Zs.index)


def quantile_map(Ms, Ts):
    """Map each METABRIC gene onto TCGA's distribution of that gene (rank -> TCGA quantile), then z-score with TCGA's mean and SD."""
    q = (Ms.rank(method="average") - 0.5) / len(Ms)
    Q = pd.DataFrame({g: np.quantile(Ts[g].values, q[g].values) for g in Ms.columns}, index=Ms.index)
    return (Q - Ts.mean()) / Ts.std(ddof=0)


def pearson(Zs, cen):
    r = np.corrcoef(np.vstack([Zs.values, cen.T.values]))[: len(Zs), len(Zs):]
    return pd.Series(r.argmax(1), index=Zs.index)


tcga_er_known = tl.loc[tc.index, "ER"].dropna()
log("ER_negative_share_by_cohort", {"TCGA_clustered_tumours_ER_IHC_known": round(float((tcga_er_known == "Negative").mean()), 3),
                                    "METABRIC_ER_IHC_known": round(float((ml.ER_IHC[er_known] == "Negative").mean()), 3)})
QZ = quantile_map(M[Gm], TX[Gm])
clus, check = {}, {}
for f in sorted(PROC.glob("tcga_kmeans_k*_centroids_z.tsv")):
    k = int(f.stem.split("kmeans_k")[1].split("_")[0])
    cen_all = pd.read_csv(f, sep="\t", index_col=0)
    # self-check on TCGA: the rule must give back the k-means labels on all genes, and nearly so on the shared genes
    ari_all = ARI(tc[f"kmeans_k{k}"], nearest(TZ[G], cen_all.loc[G]))
    ari_shared = ARI(tc[f"kmeans_k{k}"], nearest(TZ[Gm], cen_all.loc[Gm]))
    assert ari_all > 0.999, f"nearest-centroid rule does not reproduce TCGA k-means k={k} (ARI {ari_all:.3f})"
    check[f"k{k}"] = {"ARI_all_genes": round(ari_all, 4), "ARI_shared_genes": round(ari_shared, 4)}
    cen = cen_all.loc[Gm]
    a = nearest(Z, cen)
    ap = pearson(Z, cen)
    aq = nearest(QZ, cen)
    clus[f"tcga_centroid_k{k}"] = a
    clus[f"tcga_centroid_k{k}_pearson_sensitivity"] = ap
    order = ["LumA", "LumB", "Her2", "Basal", "claudin-low", "Normal"]
    ct = pd.crosstab(a.map(lambda c: f"C{c}"), ml.CLAUDIN_SUBTYPE).reindex(columns=order).fillna(0).astype(int)
    ct.to_csv(OUT / f"centroid_k{k}_vs_claudin_subtype.tsv", sep="\t")
    ct2 = pd.crosstab(a.map(lambda c: f"C{c}"), grp).reindex(columns=list(col)).fillna(0).astype(int)
    ct2.to_csv(OUT / f"centroid_k{k}_vs_ER_HER2_group.tsv", sep="\t")
    m4 = ml.CLAUDIN_SUBTYPE.isin(["LumA", "LumB", "Her2", "Basal"])
    tcga_er = tl.loc[tc.index, "ER"]
    tk = tc[f"kmeans_k{k}"]
    stats = {
        "sizes_METABRIC": a.value_counts().reindex(range(k), fill_value=0).tolist(),
        "share_METABRIC": (a.value_counts(normalize=True).reindex(range(k), fill_value=0)).round(3).tolist(),
        "share_TCGA": (tk.value_counts(normalize=True).reindex(range(k), fill_value=0)).round(3).tolist(),
        "ARI_vs_PAM50_4class": round(ARI(a[m4], ml.CLAUDIN_SUBTYPE[m4]), 3),
        "ARI_vs_ER_IHC": round(ARI(a[er_known], ml.ER_IHC[er_known]), 3),
        "ER_neg_share_METABRIC": (ml.ER_IHC[er_known] == "Negative").groupby(a[er_known]).mean().round(3).to_dict(),
        "ER_neg_share_TCGA": (tcga_er[tcga_er.notna()] == "Negative").groupby(tk[tcga_er.notna()]).mean().round(3).to_dict(),
        "pearson_sensitivity": {"sizes": ap.value_counts().reindex(range(k), fill_value=0).tolist(),
                                "ARI_vs_euclidean": round(ARI(a, ap), 3)},
        "quantile_mapping_sensitivity": {
            "share_METABRIC": aq.value_counts(normalize=True).reindex(range(k), fill_value=0).round(3).tolist(),
            "ER_neg_share_METABRIC": (ml.ER_IHC[er_known] == "Negative").groupby(aq[er_known]).mean().round(3).to_dict(),
            "ARI_vs_main": round(ARI(a, aq), 3)},
    }
    log(f"centroid_k{k}", stats)
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 0.5 * k + 2.2), gridspec_kw={"width_ratios": [3, 5]})
    for ax_, t, ttl in [(axes[0], ct2, f"vs ER (IHC) / HER2 (SNP6), n={int(ct2.values.sum()):,}: ARI vs ER {stats['ARI_vs_ER_IHC']:.2f}"),
                        (axes[1], ct, f"vs PAM50 + claudin-low (array call), n={int(ct.values.sum()):,}: ARI vs 4-class PAM50 "
                                      f"{stats['ARI_vs_PAM50_4class']:.2f}")]:
        sns.heatmap(t.div(t.sum(0).replace(0, np.nan), axis=1), annot=t.values, fmt="d", cmap="Blues", vmin=0, vmax=1,
                    cbar=False, linewidths=2, linecolor="white", ax=ax_, annot_kws={"size": 8})
        ax_.set_xlabel(""); ax_.set_ylabel(""); ax_.tick_params(axis="y", rotation=0); ax_.tick_params(axis="x", labelsize=8)
        ax_.set_xticklabels([l.get_text().replace(" not gained", "\nnot gained").replace("HER2 gained", "HER2\ngained")
                             for l in ax_.get_xticklabels()], rotation=0)
        ax_.set_title(ttl, fontsize=9.5)
    axes[0].set_ylabel("TCGA cluster (nearest centroid)")
    fig.suptitle(f"TCGA k-means k={k} centroids applied to all {len(a):,} METABRIC tumours; panels show those with each label; "
                 "cell = count, shade = share of column",
                 fontsize=10.5, fontweight="bold", x=0.01, ha="left")
    sm_ = a.value_counts(normalize=True).reindex(range(k), fill_value=0)   # unrounded, so the footer rounds once
    st_ = tk.value_counts(normalize=True).reindex(range(k), fill_value=0)
    fig.text(0.01, -0.03, "Shares of METABRIC / TCGA tumours per cluster: "
             + ", ".join(f"C{c} {sm:.0%}/{st:.0%}" for c, (sm, st) in enumerate(zip(sm_, st_))),
             fontsize=8, color=INK2)
    fig.tight_layout()
    save(fig, f"figS9_metabric_centroids_k{k}")
log("projection_rule_self_check_on_TCGA", check)
pd.DataFrame(clus).to_csv(PROC / "metabric_clusters.tsv", sep="\t")

write_params(OUT, PARAMS)
json.dump(rep, open(OUT / "metabric_report.json", "w"), indent=2, default=str)
print("done")
