"""Step 3.4a: freeze the METABRIC transfer models on TCGA, before any METABRIC label is read (plan section 4.1).

Reads only TCGA labels, TCGA expression and the METABRIC *gene list* (column names of the expression matrix).
For each target, the hyperparameters are the ones the inner CV of script 12 picked most often across the 15 outer
folds (L1-logistic). The model is refitted once on all TCGA tumours:
  shared genes -> expression/variance filter on TCGA log2 (same rule as in CV) -> z-score within TCGA ->
  SelectKBest(f_classif, k) -> L1 logistic (class_weight balanced).
METABRIC will be z-scored within its own cohort (script 14); the scaler is not carried across platforms.
Frozen threshold: probability 0.5 (balanced class weights make 0.5 the natural cut).
Writes results/transfer/frozen/*.joblib, frozen_models.json (with a SHA-256 of each model) and appends a
`metabric-freeze` entry to CHANGELOG.md.
"""
import sys, json, hashlib, datetime
from pathlib import Path
import numpy as np
import pandas as pd
import joblib
import pyarrow.parquet as pq
from sklearn.pipeline import Pipeline
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.common3 import PROC, RES3, ROOT3, SEED, load_labels, write_params3, dump  # noqa: E402

OUT = RES3 / "transfer" / "frozen"
OUT.mkdir(parents=True, exist_ok=True)
MIN_MEAN, N_VAR = 1.0, 5000
lab = load_labels()
X = pd.read_parquet(PROC / "tcga_hiseq_tumour_log2.parquet").loc[lab.index]
mb_genes = pq.read_schema(PROC / "metabric_log2.parquet").names        # gene names only, no labels, no values
shared = [g for g in X.columns if g in set(mb_genes)]
X = X[shared]
amp = pd.read_csv(PROC / "erbb2_amplicon_genes.tsv", sep="\t")
pam50 = set(pd.read_csv(PROC / "pam50_centroids_hgnc.tsv", sep="\t", index_col=0).index)
ABL_B = {"ESR1", "PGR", "ERBB2"} | set(amp.loc[amp.within_1500kb, "hgnc"].dropna()) | pam50
bp = pd.read_csv(RES3 / "classify" / "best_params_per_fold.tsv", sep="\t")

lab["tn"] = np.where(lab.erher2.eq("ERneg_HER2neg"), "Positive", np.where(lab.ER.notna() & lab.her2.notna(), "Negative", None))
TARGETS = {"ER": ("ER", "main", "a_all_genes"), "ER_ablated": ("ER", "main", "b_no_receptor_amplicon1.5Mb_PAM50"),
           "HER2": ("her2", "main", "a_all_genes"), "HER2_ablated": ("her2", "main", "b_no_receptor_amplicon1.5Mb_PAM50"),
           "ERneg_HER2neg": ("tn", "transfer_target", "a_all_genes")}
frozen = {"created_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
          "n_shared_genes": len(shared), "threshold": 0.5, "models": {}}
for name, (col, tag, abl) in TARGETS.items():
    tgt = {"ER": "ER", "her2": "HER2", "tn": "ERneg_HER2neg"}[col]
    p = bp[(bp.tag == tag) & (bp.target == tgt) & (bp.ablation == abl) & (bp.model == "L1-logistic")]
    k = int(p["select__k"].mode().iloc[0]); C = float(p["model__C"].mode().iloc[0])
    ids = lab.index[lab[col].isin(["Positive", "Negative"])]
    y = (lab.loc[ids, col] == "Positive").astype(int).to_numpy()
    genes = [g for g in shared if not (abl.startswith("b_") and g in ABL_B)]
    Xa = X.loc[ids, genes]
    mu, var = Xa.mean(), Xa.var(ddof=0)
    ok = mu.index[mu >= MIN_MEAN]
    keep = var[ok].sort_values(ascending=False).index[:N_VAR].tolist()
    Z = (Xa[keep] - Xa[keep].mean()) / Xa[keep].std(ddof=0)
    pipe = Pipeline([("select", SelectKBest(f_classif, k=min(k, len(keep)))),
                     ("model", LogisticRegression(l1_ratio=1.0, solver="liblinear", C=C, class_weight="balanced",
                                                  max_iter=5000, random_state=SEED))]).fit(Z.to_numpy(), y)
    coef = pd.Series(pipe.named_steps["model"].coef_.ravel(), index=np.array(keep)[pipe.named_steps["select"].get_support()])
    coef = coef[coef != 0].sort_values()
    path = OUT / f"{name}.joblib"
    joblib.dump({"pipeline": pipe, "genes": keep, "target_column": col}, path)
    frozen["models"][name] = {"target": tgt, "ablation": abl, "k": k, "C": C, "n_train": len(y), "n_pos": int(y.sum()),
                              "genes_after_filter": len(keep), "nonzero_weights": int(len(coef)),
                              "top_negative": coef.head(10).round(3).to_dict(), "top_positive": coef.tail(10).round(3).to_dict(),
                              "fold_choice_frequency": p.groupby(["select__k", "model__C"]).size().rename("n").reset_index().to_dict("records"),
                              "sha256": hashlib.sha256(open(path, "rb").read()).hexdigest()}
    coef.to_csv(OUT / f"{name}_weights.tsv", sep="\t", header=["weight"])
    print(name, {k_: frozen["models"][name][k_] for k_ in ["k", "C", "n_train", "n_pos", "nonzero_weights"]}, flush=True)

dump(frozen, OUT / "frozen_models.json")
write_params3(OUT, {"rule": "modal inner-CV choice of (k, C) across 15 outer folds, L1-logistic; refit on all TCGA",
                    "shared_genes": len(shared), "filter": f"mean log2 >= {MIN_MEAN}, top {N_VAR} variance on TCGA",
                    "scaling": "z-score within cohort", "threshold": 0.5})
cl = ROOT3 / "CHANGELOG.md"
entry = (f"\n## metabric-freeze ({frozen['created_utc']})\n\nTransfer models frozen on TCGA before any METABRIC label was "
         f"read (script 13). Shared genes: {len(shared)}. Threshold 0.5. Models and SHA-256:\n\n" +
         "\n".join(f"- `{n}`: k = {m['k']}, C = {m['C']}, {m['nonzero_weights']} non-zero weights, sha256 `{m['sha256'][:16]}…`"
                   for n, m in frozen["models"].items()) + "\n")
txt = cl.read_text() if cl.exists() else "# Step 3 changelog\n"
cl.write_text(txt + entry)
