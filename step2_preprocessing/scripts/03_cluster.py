"""Step 3: unsupervised clustering of tumours and agreement with clinical groups and PAM50 (SQ1).

No label is used to build any clustering; labels are only compared afterwards.
This script does the computation (about 10-15 min, mostly consensus clustering);
scripts/03b_cluster_figures.py draws the figures from the saved tables.
Outputs: results/cluster/*, data/processed/tcga_clusters.tsv, data/processed/tcga_kmeans_k*_centroids_z.tsv
Run from the project root:  python scripts/03_cluster.py
"""
import json, sys, time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from sklearn.cluster import KMeans
from sklearn.metrics import (adjusted_rand_score as ARI, normalized_mutual_info_score as NMI,
                             silhouette_score, silhouette_samples, calinski_harabasz_score)

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.common import PROC, RES, SEED, GROUP_ORDER, write_params, drop_excluded

N_TOP, KS, N_INIT, N_RESAMPLE, SUB_FRAC, N_BOOT = 2000, list(range(2, 9)), 50, 100, 0.8, 1000
PARAMS = {
    "samples": "primary tumours not excluded (see 01 exclusions); males included (unsupervised)",
    "expressed_gene_filter": "mean log2 >= 1",
    "main_gene_set": f"top {N_TOP:,} variance genes after removing the 50 PAM50 genes and sex-specific genes (chrY, XIST, TSIX)",
    "sensitivity_gene_sets": ["top 2,000 incl. PAM50 genes (variance-ranked)", "main 2,000 plus all 50 PAM50 genes forced in",
                              "top 5,000 excl. PAM50 genes", "top 1,000 excl. PAM50 genes"],
    "scaling": "z-score per gene across tumours",
    "kmeans": {"k": KS, "n_init": N_INIT, "random_state": SEED},
    "hierarchical": ["Ward linkage, Euclidean distance on z-scores", "average linkage, 1 - Pearson correlation"],
    "consensus": {"base": "k-means (n_init=5, random_state = resample index)", "resamples": N_RESAMPLE, "sample_fraction": SUB_FRAC,
                  "same subsample used for every k in a resample": True,
                  "final_partition": "average linkage on 1 - consensus, cut at k", "PAC_interval": [0.1, 0.9]},
    "k_selection_rule": "pre-specified: primary k = the k in 2..8 with the lowest PAC (Senbabaoglu 2014); k = 4 (next PAC minimum) and k = 5 (Netanely 2016 used K-means K = 5 with correlation distance on the top 2,000 genes, with adjacent normals included) also reported",
    "agreement": f"ARI and NMI against clinical group (3 classes), PanCan PAM50 without Normal-like (4 classes), TCGA 2012 SigClust clusters (13 levels) and IHC ER status; 95% CI from {N_BOOT} bootstrap resamples of tumours (seed 42). The CIs are conditional on the fitted partition (clusters are not re-fitted in each resample); clustering stability is reported separately (PAC, method and gene-set ARIs). ARI depends on the number of classes, so label-vs-label agreement is reported as a reference",
    "cluster_ids": "relabelled by size: C0 = largest cluster",
}
OUT = RES / "cluster"
OUT.mkdir(parents=True, exist_ok=True)
rep = {}
rng = np.random.default_rng(SEED)


def log(k, v):
    rep[k] = v
    print(f"{k}: {v}", flush=True)


lab = drop_excluded(pd.read_csv(PROC / "tcga_labels.tsv", sep="\t", index_col=0))
X = pd.read_parquet(PROC / "tcga_hiseq_tumour_log2.parquet").loc[lab.index].astype("float64")
p50 = list(pd.read_csv(PROC / "pam50_centroids_hgnc.tsv", sep="\t", index_col=0).index)
sexg = set(pd.read_csv(PROC / "sex_specific_genes.txt", header=None)[0])
Xe = X.loc[:, (X.mean() >= 1) & ~X.columns.isin(sexg)]
var = Xe.var().sort_values(ascending=False)


def gene_set(n, drop_pam50=True):
    v = var[[g for g in var.index if not (drop_pam50 and g in p50)]]
    return list(v.index[:n])


def z(M):
    return (M - M.mean()) / M.std(ddof=0)


G = gene_set(N_TOP)
Z = z(Xe[G])
Zv = Z.values
n = len(Z)
log("n_tumours", n)
log("n_genes_main", len(G))
log("pam50_genes_that_pass_top2000_variance_cut", int(len(set(gene_set(N_TOP, False)) & set(p50))))
pd.Series(G).to_csv(OUT / "cluster_genes_main.txt", index=False, header=False)

# ------------------------------------------------ k-means scan
km, scan = {}, []
for k in KS:
    m = KMeans(n_clusters=k, n_init=N_INIT, random_state=SEED).fit(Zv)
    km[k] = m
    scan.append((k, m.inertia_, silhouette_score(Zv, m.labels_), calinski_harabasz_score(Zv, m.labels_)))
scan = pd.DataFrame(scan, columns=["k", "inertia", "silhouette", "calinski_harabasz"]).set_index("k")

# ------------------------------------------------ consensus clustering
t0 = time.time()
M = {k: np.zeros((n, n), np.float32) for k in KS}
I = np.zeros((n, n), np.float32)
m_sub = int(round(SUB_FRAC * n))
for r in range(N_RESAMPLE):
    idx = np.sort(rng.choice(n, m_sub, replace=False))
    I[np.ix_(idx, idx)] += 1
    for k in KS:
        l = KMeans(n_clusters=k, n_init=5, random_state=r).fit_predict(Zv[idx])
        oh = np.eye(k, dtype=np.float32)[l]
        M[k][np.ix_(idx, idx)] += oh @ oh.T
log("consensus_runtime_s", round(time.time() - t0))
iu = np.triu_indices(n, 1)
pac, cons_lab, cdf = {}, {}, {}
for k in KS:
    C = M[k] / np.maximum(I, 1)
    np.fill_diagonal(C, 1)
    c = C[iu][I[iu] > 0]
    pac[k] = float(((c > 0.1) & (c < 0.9)).mean())
    cdf[k] = np.quantile(c, np.linspace(0, 1, 501))
    cons_lab[k] = fcluster(linkage(squareform(1 - C, checks=False), "average"), k, "maxclust") - 1
pd.DataFrame(cdf, index=np.linspace(0, 1, 501)).to_csv(OUT / "consensus_cdf_quantiles.tsv", sep="\t")
scan["PAC"] = pd.Series(pac)
scan.to_csv(OUT / "k_selection.tsv", sep="\t")
log("k_scan", scan.round(4).to_dict(orient="index"))
k_star = int(scan.PAC.idxmin())
log("k_star_lowest_PAC", k_star)
K_REPORT = sorted({k_star, 4, 5})

# ------------------------------------------------ hierarchical
Lw = linkage(Zv, "ward")
La = linkage(squareform(1 - np.corrcoef(Zv), checks=False), "average")


def by_size(s):
    s = pd.Series(s, index=Z.index)
    return s.map({c: i for i, c in enumerate(s.value_counts().index)})


clus = pd.DataFrame(index=Z.index)
for k in sorted(set(K_REPORT) | {2, 3}):
    clus[f"kmeans_k{k}"] = by_size(km[k].labels_)
    clus[f"ward_k{k}"] = by_size(fcluster(Lw, k, "maxclust") - 1)
    clus[f"avgcorr_k{k}"] = by_size(fcluster(La, k, "maxclust") - 1)
    clus[f"consensus_k{k}"] = by_size(cons_lab[k])
log("cluster_sizes", {c: clus[c].value_counts().sort_index().tolist() for c in clus})
log("method_agreement_ARI", {f"k{k}": {"kmeans_vs_ward": round(ARI(clus[f'kmeans_k{k}'], clus[f'ward_k{k}']), 3),
                                       "kmeans_vs_consensus": round(ARI(clus[f'kmeans_k{k}'], clus[f'consensus_k{k}']), 3),
                                       "kmeans_vs_avgcorr": round(ARI(clus[f'kmeans_k{k}'], clus[f'avgcorr_k{k}']), 3)}
                             for k in K_REPORT})

# ------------------------------------------------ agreement with labels (bootstrap CIs)
refs = {
    "Clinical group": lab.clinical_group.where(lab.clinical_group.isin(GROUP_ORDER)),
    "PAM50 (no Normal-like)": lab.PAM50_pancan.where(lab.PAM50_pancan.isin(["LumA", "LumB", "Her2", "Basal"])),
    "TCGA 2012 SigClust": lab.SigClust_2012,
    "ER status (IHC)": lab.ER,
}


def boot(a, b):
    m = a.notna() & b.notna()
    a, b = a[m].astype(str).values, b[m].astype(str).values
    r = np.random.default_rng(SEED)
    ar, nm = [], []
    for _ in range(N_BOOT):
        i = r.integers(0, len(a), len(a))
        ar.append(ARI(a[i], b[i])); nm.append(NMI(a[i], b[i]))
    return (len(a), ARI(a, b), *np.percentile(ar, [2.5, 97.5]), NMI(a, b), *np.percentile(nm, [2.5, 97.5]))


def agreement(cl, tag):
    return [(tag, rn, *boot(cl, ref.reindex(cl.index))) for rn, ref in refs.items()]


rows = []
for c in clus:
    rows += agreement(clus[c], c)
for rn, ref in refs.items():   # label-vs-label reference points
    if rn != "PAM50 (no Normal-like)":
        rows.append(("PAM50 (no Normal-like)", rn, *boot(refs["PAM50 (no Normal-like)"], ref)))
cols_ = ["clustering", "reference", "n", "ARI", "ARI_lo", "ARI_hi", "NMI", "NMI_lo", "NMI_hi"]
ag = pd.DataFrame(rows, columns=cols_)

# ------------------------------------------------ sensitivity: other gene sets
sens = []
for name, gs in [("top2000_incl_PAM50", gene_set(N_TOP, False)), ("main_plus_all_PAM50", G + p50),
                 ("top5000_excl_PAM50", gene_set(5000)), ("top1000_excl_PAM50", gene_set(1000))]:
    gs = [g for g in dict.fromkeys(gs) if g in X.columns]
    Zs = z(X[gs]).values
    for k in K_REPORT:
        l = by_size(KMeans(n_clusters=k, n_init=N_INIT, random_state=SEED).fit_predict(Zs))
        clus[f"sens_{name}_kmeans_k{k}"] = l
        sens += agreement(l, f"sens_{name}_kmeans_k{k}")
        rep.setdefault("sensitivity_ARI_vs_main", {})[f"{name}_k{k}"] = round(ARI(l, clus[f"kmeans_k{k}"]), 3)
ag = pd.concat([ag, pd.DataFrame(sens, columns=cols_)])
ag.to_csv(OUT / "cluster_agreement.tsv", sep="\t", index=False)
log("sensitivity_ARI_vs_main", rep["sensitivity_ARI_vs_main"])
show = ag[ag.clustering.isin([f"kmeans_k{k}" for k in K_REPORT] + [f"consensus_k{k_star}", f"ward_k{k_star}"])]
log("agreement_main", show.round(3).to_dict(orient="records"))

clus.to_csv(PROC / "tcga_clusters.tsv", sep="\t")

# ------------------------------------------------ centroids for projecting other cohorts (z-units, main gene set)
for k in K_REPORT:
    cen = Z.groupby(clus[f"kmeans_k{k}"]).mean().T
    cen.columns = [f"C{c}" for c in cen.columns]
    cen.to_csv(PROC / f"tcga_kmeans_k{k}_centroids_z.tsv", sep="\t")


# ------------------------------------------------ cluster profiles (shares over tumours with a complete clinical group)
def describe(col):
    d = []
    for c, idx in clus.groupby(col).groups.items():
        L = lab.loc[idx]
        Lc = L[L.clinical_group.isin(GROUP_ORDER)]
        cg = Lc.clinical_group.value_counts(normalize=True)
        pm = L.PAM50_pancan.value_counts(normalize=True)
        top_up = (Z.loc[idx].mean() - Z.drop(idx).mean()).sort_values(ascending=False)
        d.append({"cluster": f"C{c}", "n": len(idx), "n_with_clinical_group": len(Lc),
                  "HR+/HER2-": round(cg.get("HR+/HER2-", 0), 3), "HER2+": round(cg.get("HER2+", 0), 3), "TNBC": round(cg.get("TNBC", 0), 3),
                  "top_PAM50": pm.index[0], "top_PAM50_share": round(pm.iloc[0], 3),
                  "PAM50": "; ".join(f"{k} {v:.0%}" for k, v in pm.items()),
                  "ILC_share": round((L.histology == "ILC").mean(), 3), "ER_neg_share": round((L.ER == "Negative").sum() / max(L.ER.notna().sum(), 1), 3),
                  "top_genes_up": ", ".join(top_up.index[:10])})
    return pd.DataFrame(d)


for k in K_REPORT:
    dd = describe(f"kmeans_k{k}")
    dd.to_csv(OUT / f"cluster_profile_kmeans_k{k}.tsv", sep="\t", index=False)
    log(f"profile_kmeans_k{k}", dd.drop(columns=["PAM50"]).to_dict(orient="records"))

sil = silhouette_samples(Zv, clus[f"kmeans_k{k_star}"].values)
log("silhouette_negative_fraction_kstar", round(float((sil < 0).mean()), 3))
log(f"clinical_group_by_kmeans_k{k_star}", pd.crosstab(lab.clinical_group, clus[f"kmeans_k{k_star}"]).to_dict(orient="index"))

write_params(OUT, PARAMS)
json.dump(rep, open(OUT / "cluster_report.json", "w"), indent=2, default=str)
print("done")
