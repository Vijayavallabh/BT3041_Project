"""Step 3.4b: external validation in METABRIC (plan section 4): classifier transfer and DE replication.

Transfer: the models frozen by script 13 (CHANGELOG `metabric-freeze`) predict METABRIC after z-scoring each gene
within METABRIC. Labels from other assays only: ER_IHC (IHC) and HER2_SNP6 GAIN (copy number); the ER-/HER2-
group is ER_IHC negative and HER2 not gained (PR is not used; not called TNBC).
Primary metric AUC with a 2,000-sample bootstrap CI; balanced accuracy and MCC at the frozen threshold (0.5) and at a
prevalence-matched threshold (the share called positive equals the METABRIC prevalence).
Robustness: quantile-map each METABRIC gene onto the TCGA distribution, then z-score with TCGA's mean and SD.
Fixed feature set: the TCGA DE genes of the matched contrast (ER-/HER2- vs ER+/HER2-) as features of an L1 logistic
fitted on TCGA (legitimate here because METABRIC is a true external test).
Baselines: ESR1 alone, ERBB2 alone. Agreement with METABRIC's array-derived ER_STATUS / HER2_STATUS is reported
only as agreement with an expression-threshold call.
Secondary (circular comparator): 4-class PAM50 model from TCGA vs METABRIC CLAUDIN_SUBTYPE (macro-F1, kappa).

DE replication: OLS per gene, log2 intensity ~ group + age, ER-IHC-/HER2-not-gained vs ER-IHC+/HER2-not-gained,
BH over genes; compared with TCGA `matched_erher2` (share replicated with the same sign at q < 0.05, Spearman of
effects, Hallmark NES agreement via prerank on the METABRIC t statistic).
"""
import sys, json
from pathlib import Path
import numpy as np
import pandas as pd
import joblib
from scipy import stats
from sklearn.metrics import roc_auc_score, balanced_accuracy_score, matthews_corrcoef, f1_score, cohen_kappa_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
import gseapy as gp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.common3 import PROC, RES3, SEED, load_labels, read_gmt, bh, write_params3, dump  # noqa: E402

OUT = RES3 / "transfer"
FR = OUT / "frozen"
frozen = json.load(open(FR / "frozen_models.json"))
rng = np.random.default_rng(SEED)
N_BOOT = 2000

lab = load_labels()
Xt = pd.read_parquet(PROC / "tcga_hiseq_tumour_log2.parquet").loc[lab.index]
Xm = pd.read_parquet(PROC / "metabric_log2.parquet")
ml = pd.read_csv(PROC / "metabric_labels.tsv", sep="\t", index_col=0).loc[Xm.index]
shared = [g for g in Xt.columns if g in set(Xm.columns)]
Xt, Xm = Xt[shared], Xm[shared]
Zm = (Xm - Xm.mean()) / Xm.std(ddof=0)
Zm = Zm.fillna(0.0)

# quantile mapping of METABRIC onto TCGA (per gene), then TCGA z-scores
def qmap(Xsrc, Xref):
    out = np.empty_like(Xsrc.to_numpy(dtype=float))
    R = np.sort(Xref.to_numpy(dtype=float), axis=0)
    q = (Xsrc.rank(pct=True).to_numpy() - 0.5 / len(Xsrc))
    nref = R.shape[0]
    for j in range(Xsrc.shape[1]):
        out[:, j] = np.interp(q[:, j], (np.arange(nref) + 0.5) / nref, R[:, j])
    return pd.DataFrame(out, index=Xsrc.index, columns=Xsrc.columns)


Qm = qmap(Xm, Xt)
Zq = ((Qm - Xt.mean()) / Xt.std(ddof=0)).fillna(0.0)

y_targets = {
    "ER": ml.ER_IHC.map({"Positive": 1, "Negative": 0}),
    "HER2": ml.HER2_gain,
    "ERneg_HER2neg": ml.contrast_group.map({"ER-/HER2-notgained": 1, "ER+/HER2-notgained": 0, "HER2-gained": 0}),
}
y_targets["ER_ablated"], y_targets["HER2_ablated"] = y_targets["ER"], y_targets["HER2"]


def boot_auc(y, s):
    y, s = np.asarray(y), np.asarray(s)
    v = []
    for _ in range(N_BOOT):
        i = rng.integers(0, len(y), len(y))
        if y[i].min() != y[i].max():
            v.append(roc_auc_score(y[i], s[i]))
    return [float(x) for x in np.percentile(v, [2.5, 97.5])]


def thr_metrics(y, s, thr):
    p = (s > thr).astype(int)
    return {"threshold": float(thr), "called_pos": int(p.sum()), "balanced_acc": balanced_accuracy_score(y, p),
            "MCC": matthews_corrcoef(y, p), "F1": f1_score(y, p)}


cv = pd.read_csv(RES3 / "classify" / "summary.tsv", sep="\t")
rows, preds = [], pd.DataFrame(index=Xm.index)
for name, info in frozen["models"].items():
    m = joblib.load(FR / f"{name}.joblib")
    y = y_targets[name].dropna().astype(int)
    for scaling, Z in [("z_within_cohort", Zm), ("quantile_map_to_TCGA", Zq)]:
        s = pd.Series(m["pipeline"].predict_proba(Z.loc[y.index, m["genes"]].to_numpy())[:, 1], index=y.index)
        if scaling == "z_within_cohort":
            preds[f"p_{name}"] = pd.Series(m["pipeline"].predict_proba(Z.loc[:, m["genes"]].to_numpy())[:, 1], index=Xm.index)
        auc = roc_auc_score(y, s)
        prev_thr = np.quantile(s, 1 - y.mean())
        tag = "transfer_target" if name == "ERneg_HER2neg" else "main"
        cvrow = cv[(cv.tag == tag) & (cv.target == info["target"]) & (cv.ablation == info["ablation"]) & (cv.model == "L1-logistic")]
        rows.append({"model": name, "scaling": scaling, "n": len(y), "n_pos": int(y.sum()), "prevalence": round(y.mean(), 3),
                     "AUC": auc, "AUC_CI": boot_auc(y, s), "TCGA_CV_AUC": float(cvrow.AUC_pooled.iloc[0]) if len(cvrow) else None,
                     "frozen_0.5": thr_metrics(y, s, 0.5), "prevalence_matched": thr_metrics(y, s, prev_thr)})
        print(name, scaling, f"AUC {auc:.3f}", flush=True)

# single-gene baselines in METABRIC (no fitting; AUC is threshold-free)
base = {}
for tgt, gene, sign in [("ER", "ESR1", 1), ("HER2", "ERBB2", 1), ("ERneg_HER2neg", "ESR1", -1)]:
    y = y_targets[tgt].dropna().astype(int)
    s = sign * Xm.loc[y.index, gene]
    base[tgt] = {"gene": gene, "AUC": roc_auc_score(y, s), "AUC_CI": boot_auc(y, s)}

# DE genes of the matched TCGA contrast as a fixed feature set (fit on TCGA, test on METABRIC)
de = pd.read_csv(RES3 / "de" / "matched_erher2_results.tsv", sep="\t", index_col=0)
de_genes = sorted(set(de.symbol[de.DE]) & set(shared))
lab["tn"] = np.where(lab.erher2.eq("ERneg_HER2neg"), 1, np.where(lab.ER.notna() & lab.her2.notna(), 0, np.nan))
yt = lab.tn.dropna().astype(int)
Zt = (Xt - Xt.mean()) / Xt.std(ddof=0)
C_tn = frozen["models"]["ERneg_HER2neg"]["C"]
lr = LogisticRegression(l1_ratio=1.0, solver="liblinear", C=C_tn, class_weight="balanced", max_iter=5000, random_state=SEED)
lr.fit(Zt.loc[yt.index, de_genes].fillna(0).to_numpy(), yt)
y = y_targets["ERneg_HER2neg"].dropna().astype(int)
s = lr.predict_proba(Zm.loc[y.index, de_genes].to_numpy())[:, 1]
fixed_de = {"n_genes": len(de_genes), "C": C_tn, "nonzero": int((lr.coef_ != 0).sum()), "AUC": roc_auc_score(y, s),
            "AUC_CI": boot_auc(y, s), "frozen_0.5": thr_metrics(y, s, 0.5)}

# agreement with METABRIC's array-derived calls (not independent)
agree = {}
for name, col in [("ER", "ER_STATUS"), ("HER2", "HER2_STATUS")]:
    ok = ml[col].isin(["Positive", "Negative"])
    yy = (ml.loc[ok, col] == "Positive").astype(int)
    pp = (preds.loc[ok, f"p_{name}"] > 0.5).astype(int)
    agree[f"{name}_vs_{col}"] = {"n": int(ok.sum()), "agreement": float((yy == pp).mean()), "kappa": cohen_kappa_score(yy, pp),
                                 "AUC": roc_auc_score(yy, preds.loc[ok, f"p_{name}"])}

# secondary: 4-class PAM50 trained on TCGA (PanCan calls) vs METABRIC CLAUDIN_SUBTYPE
p50 = lab.PAM50_pancan
idx = p50.index[p50.isin(["LumA", "LumB", "Her2", "Basal"])]
mu, var = Xt.loc[idx].mean(), Xt.loc[idx].var()
keep = var[mu >= 1].sort_values(ascending=False).index[:5000]
pam_pipe = Pipeline([("select", SelectKBest(f_classif, k=1000)),
                     ("model", LogisticRegression(l1_ratio=1.0, solver="saga", C=0.1, max_iter=3000, class_weight="balanced",
                                                  random_state=SEED))])
pam_pipe.fit(Zt.loc[idx, keep].fillna(0).to_numpy(), p50.loc[idx])
cm = ml.CLAUDIN_SUBTYPE.isin(["LumA", "LumB", "Her2", "Basal"])
pp = pam_pipe.predict(Zm.loc[cm, keep].to_numpy())
pam_sec = {"n": int(cm.sum()), "accuracy": float((pp == ml.CLAUDIN_SUBTYPE[cm]).mean()),
           "macro_F1": f1_score(ml.CLAUDIN_SUBTYPE[cm], pp, average="macro"), "kappa": cohen_kappa_score(ml.CLAUDIN_SUBTYPE[cm], pp),
           "confusion": pd.crosstab(ml.CLAUDIN_SUBTYPE[cm], pp).to_dict()}
preds["pred_PAM50_from_TCGA"] = pam_pipe.predict(Zm.loc[:, keep].to_numpy())
preds["pred_ERneg_HER2neg"] = (preds.p_ERneg_HER2neg > 0.5).astype(int)
preds.to_csv(OUT / "metabric_predictions.tsv", sep="\t")

tr = pd.DataFrame(rows)
tr.to_csv(OUT / "transfer_metrics.tsv", sep="\t", index=False)

# ------------------------------------------------ DE replication
grp = ml.contrast_group.map({"ER-/HER2-notgained": 1, "ER+/HER2-notgained": 0})
d = ml.assign(g=grp).dropna(subset=["g", "AGE_AT_DIAGNOSIS"])
Y = Xm.loc[d.index].to_numpy(dtype=float)
D = np.column_stack([np.ones(len(d)), d.g.to_numpy(), (d.AGE_AT_DIAGNOSIS - d.AGE_AT_DIAGNOSIS.mean()).to_numpy()])
beta, *_ = np.linalg.lstsq(D, Y, rcond=None)
resid = Y - D @ beta
dof = len(d) - D.shape[1]
s2 = (resid ** 2).sum(0) / dof
se = np.sqrt(s2 * np.linalg.inv(D.T @ D)[1, 1])
t = beta[1] / se
p = 2 * stats.t.sf(np.abs(t), dof)
mrep = pd.DataFrame({"beta_ERneg": beta[1], "t": t, "p": p, "q_BH": bh(p)}, index=shared)
mrep.index.name = "symbol"
mrep.to_csv(OUT / "metabric_de_ols.tsv", sep="\t")
tde = de.dropna(subset=["symbol"]).drop_duplicates("symbol").set_index("symbol")
com = tde.index.intersection(mrep.index)
tdeg = tde.loc[com][tde.loc[com].DE]
rep_same = (mrep.loc[tdeg.index].q_BH < 0.05) & (np.sign(mrep.loc[tdeg.index].beta_ERneg) == np.sign(tdeg.log2FC))
repl = {"n_ERneg": int(d.g.sum()), "n_ERpos": int((d.g == 0).sum()), "genes_shared": len(com),
        "TCGA_DEGs_on_array": int(len(tdeg)), "replicated_same_sign_q05": int(rep_same.sum()), "fraction": float(rep_same.mean()),
        "sign_agreement": float((np.sign(mrep.loc[tdeg.index].beta_ERneg) == np.sign(tdeg.log2FC)).mean()),
        "spearman_effect_all_shared": stats.spearmanr(tde.loc[com, "log2FC_shrunk"], mrep.loc[com, "beta_ERneg"]).statistic,
        "spearman_effect_TCGA_DEGs": stats.spearmanr(tdeg.log2FC_shrunk, mrep.loc[tdeg.index, "beta_ERneg"]).statistic,
        "metabric_genes_q05": int((mrep.q_BH < 0.05).sum()),
        "metabric_genes_q05_abs_beta_gt1": int(((mrep.q_BH < 0.05) & (mrep.beta_ERneg.abs() > 1)).sum())}
up_rep = rep_same[tdeg.log2FC > 0].mean(); dn_rep = rep_same[tdeg.log2FC < 0].mean()
repl["fraction_up"], repl["fraction_down"] = float(up_rep), float(dn_rep)
# Hallmark NES in METABRIC vs TCGA matched contrast
rnk = mrep.t.sort_values(ascending=False)
pre = gp.prerank(rnk=rnk, gene_sets=read_gmt(), permutation_num=1000, seed=SEED, threads=4, outdir=None, min_size=15,
                 max_size=500, verbose=False)
nm = pre.res2d.assign(Term=lambda x: x.Term.str.replace("HALLMARK_", "", regex=False)).set_index("Term")
nm[["NES", "FDR q-val", "Lead_genes"]].to_csv(OUT / "metabric_gsea_prerank.tsv", sep="\t")
G = pd.read_csv(RES3 / "enrichment" / "gsea_prerank_all_runs.tsv", sep="\t")
gt = G[G.run == "matched_erher2"].set_index("Term")
j = gt[["NES", "FDR_q"]].join(nm[["NES", "FDR q-val"]].astype(float), rsuffix="_MB")
repl["hallmark_NES_spearman"] = stats.spearmanr(j.NES, j.NES_MB).statistic
sig = j[(j.FDR_q < 0.05)]
repl["hallmark_TCGA_FDR05_same_sign_in_METABRIC"] = f"{int((np.sign(sig.NES) == np.sign(sig.NES_MB)).sum())}/{len(sig)}"
repl["hallmark_TCGA_FDR05_same_sign_and_FDR05_in_METABRIC"] = int(((np.sign(sig.NES) == np.sign(sig.NES_MB)) & (sig["FDR q-val"] < 0.05)).sum())
j.to_csv(OUT / "hallmark_nes_tcga_vs_metabric.tsv", sep="\t")

report = {"shared_genes": len(shared), "transfer": rows, "single_gene_baselines": base, "fixed_DE_feature_set": fixed_de,
          "agreement_with_array_derived_calls": agree, "pam50_secondary": pam_sec, "de_replication": repl}
dump(report, OUT / "metabric_report.json")
write_params3(OUT, {"bootstrap": N_BOOT, "scaling_main": "z-score within METABRIC", "robustness": "quantile map to TCGA then TCGA z",
                    "prevalence_matched_threshold": "quantile of scores so share called positive = METABRIC prevalence",
                    "de_replication": "OLS log2 ~ group + age (centred), BH", "pam50_secondary": "TCGA PanCan 4-class, L1 multinomial k=1000 C=0.1"},
              name="params_metabric.json")
print(json.dumps({"base": base, "fixed_de": fixed_de, "repl": repl, "pam": {k: pam_sec[k] for k in ["accuracy", "macro_F1", "kappa"]}},
                 indent=1, default=float))
