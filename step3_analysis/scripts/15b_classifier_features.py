"""Step 3.5b: what the classifiers use once the receptor genes, the 17q12 amplicon and PAM50 are removed (plan 3D/3E).

Stable features: genes with a non-zero L1-logistic weight in >= 80% of the 15 outer folds (12/15).
ORA (Fisher, one-sided) of the stable features against the 50 Hallmarks; universe = genes that passed the in-fold
expression/variance filter in >= 80% of folds (the genes the model could have chosen). BH over the 50 sets.
For ER the remaining features are read as the ER regulon (GATA3, AGR3, XBP1, ...), not as new markers.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.common3 import RES3, PROC, SEED, read_gmt, bh, dump, load_labels  # noqa: E402
from sklearn.model_selection import RepeatedStratifiedKFold  # noqa: E402

lab = load_labels()
Xall = pd.read_parquet(PROC / "tcga_hiseq_tumour_log2.parquet").loc[lab.index]
amp = pd.read_csv(PROC / "erbb2_amplicon_genes.tsv", sep="\t")
pam50 = set(pd.read_csv(PROC / "pam50_centroids_hgnc.tsv", sep="\t", index_col=0).index)
ABL_B = {"ESR1", "PGR", "ERBB2"} | set(amp.loc[amp.within_1500kb, "hgnc"].dropna()) | pam50


def target_ids(target):
    col = {"ER": "ER", "HER2": "her2"}[target]
    ids = lab.index[lab[col].isin(["Positive", "Negative"])]
    return ids, (lab.loc[ids, col] == "Positive").astype(int).to_numpy()

fe = pd.read_csv(RES3 / "classify" / "feature_selection_frequency.tsv", sep="\t")
H = read_gmt()
N_FOLDS, MIN = 15, 12
report = {}
rows = []
for target in ["ER", "HER2"]:
    for abl in ["a_all_genes", "b_no_receptor_amplicon1.5Mb_PAM50"]:
        d = fe[(fe.tag == "main") & (fe.target == target) & (fe.ablation == abl) & (fe.model == "L1-logistic")]
        stable = set(d.gene[d.folds_selected >= MIN])
        # universe: rerun the (label-free) in-fold expression/variance filter on the same 15 outer training folds
        ids, y = target_ids(target)
        X = Xall.loc[ids, [g for g in Xall.columns if not (abl.startswith("b_") and g in ABL_B)]]
        cnt = pd.Series(0, index=X.columns)
        for tr, _ in RepeatedStratifiedKFold(n_splits=5, n_repeats=3, random_state=SEED).split(np.zeros(len(y)), y):
            Xt = X.iloc[tr]
            mu, var = Xt.mean(), Xt.var(ddof=0)
            cnt[var[mu >= 1.0].sort_values(ascending=False).index[:5000]] += 1
        universe = set(cnt.index[cnt >= MIN]) | stable
        res = []
        for term, gs in H.items():
            gs = set(gs) & universe
            a = len(stable & gs); b = len(stable - gs); c = len(gs - stable); dd = len(universe) - a - b - c
            orr, p = stats.fisher_exact([[a, b], [c, dd]], alternative="greater")
            res.append({"target": target, "ablation": abl, "Term": term.replace("HALLMARK_", ""), "overlap": a, "set_in_universe": len(gs),
                        "odds_ratio": orr, "p": p, "genes": ";".join(sorted(stable & gs))})
        r = pd.DataFrame(res)
        r["q_BH"] = bh(r.p)
        rows.append(r)
        top = d.sort_values(["folds_selected"], ascending=False)
        report[f"{target}_{abl}"] = {"stable_features": len(stable), "universe": len(universe),
                                     "top_stable": top.gene[top.folds_selected >= MIN].head(25).tolist(),
                                     "hallmarks_q_lt_0.05": r[r.q_BH < 0.05].sort_values("p").Term.tolist()}
        print(target, abl, report[f"{target}_{abl}"], flush=True)
pd.concat(rows).sort_values(["target", "ablation", "p"]).to_csv(RES3 / "classify" / "stable_feature_ora.tsv", sep="\t", index=False)
dump(report, RES3 / "classify" / "stable_feature_report.json")
