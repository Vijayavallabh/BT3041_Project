"""Step 3.3: predicting clinical ER and HER2 status from tumour expression (plan section 3E).

Targets (female supervised set, step 2 exclusions, HiSeqV2 log2):  IHC ER (Positive = 1) and IHC/FISH HER2
(Positive = 1); plus ER-/HER2- vs the rest (PR ignored) for METABRIC transfer.
Models: L1 logistic regression, linear SVM, Random Forest (500 trees), PLS-DA.
Ablations: (a) all genes; (b) without ESR1, PGR, ERBB2, the 17q12 amplicon (every gene starting within 1.5 Mb of
ERBB2) and the 50 PAM50 genes; (c, HER2 only) as (b) with a 2 Mb window.
Nested CV: outer stratified 5-fold x 3 repeats; inner stratified 5-fold grid search on ROC-AUC.
Everything that learns from data sits inside the Pipeline, so it is refitted on each training fold only:
expression filter -> top-variance filter -> StandardScaler -> SelectKBest(f_classif, k tuned) -> model.
Baseline: the receptor gene alone (ESR1 / ERBB2), threshold by Youden's J on the training fold.
Tests: Nadeau-Bengio corrected resampled t-test on the 15 paired fold AUCs (model vs baseline, and ablation b vs a);
paired bootstrap of the pooled out-of-fold AUC as a second check. Controls: shuffled labels, majority class.
Secondary (labelled circular): PanCan PAM50, 4- and 5-class, L1 logistic, all genes vs non-PAM50 genes.
Sensitivity: HER2 with low-confidence positives dropped, and with every low-confidence call dropped.
Run time: about 1 h on 4 cores.
"""
import sys, time, json, warnings
from pathlib import Path
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy import stats
from sklearn.base import BaseEstimator, ClassifierMixin, TransformerMixin, clone
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.svm import LinearSVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.cross_decomposition import PLSRegression
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold, GridSearchCV
from sklearn.metrics import roc_auc_score, balanced_accuracy_score, matthews_corrcoef, f1_score, accuracy_score, cohen_kappa_score

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.common3 import PROC, RES3, SEED, load_labels, write_params3, dump  # noqa: E402

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
OUT = RES3 / "classify"
OUT.mkdir(parents=True, exist_ok=True)
N_JOBS = 4
OUTER_SPLITS, OUTER_REPEATS, INNER_SPLITS = 5, 3, 5
N_VAR = 5000          # top-variance genes kept inside each training fold (before scaling and SelectKBest)
MIN_MEAN = 1.0        # expression filter: mean log2(RSEM+1) >= 1 in the training fold
K_GRID = [50, 200, 1000]
N_BOOT = 2000
import os  # noqa: E402
if os.environ.get("SMOKE"):          # quick pipeline check: 1 repeat, small bootstrap
    OUTER_REPEATS, N_BOOT = 1, 50


# ---------------------------------------------------------------- pipeline pieces
class ExprVarFilter(BaseEstimator, TransformerMixin):
    """Keep genes with training-fold mean >= min_mean, then the n_var most variable of them (fit on the fold only)."""
    def __init__(self, min_mean=MIN_MEAN, n_var=N_VAR):
        self.min_mean, self.n_var = min_mean, n_var

    def fit(self, X, y=None):
        X = np.asarray(X)
        mu, var = X.mean(0), X.var(0)
        ok = np.where(mu >= self.min_mean)[0]
        self.idx_ = ok[np.argsort(var[ok])[::-1][: self.n_var]]
        return self

    def transform(self, X):
        return np.asarray(X)[:, self.idx_]


class PLSDA(BaseEstimator, ClassifierMixin):
    """PLS-DA: PLS regression on the 0/1 label; the score is the predicted y.
    Threshold: midway between the two class means of the training-fold scores (a plain 0.5 cut would favour the
    majority class when classes are unbalanced; plan deviation, see README)."""
    def __init__(self, n_components=2):
        self.n_components = n_components

    def fit(self, X, y):
        self.classes_ = np.unique(y)
        self.pls_ = PLSRegression(n_components=self.n_components, scale=False).fit(X, y.astype(float))
        s = self.pls_.predict(X).ravel()
        self.thr_ = 0.5 * (s[y == 1].mean() + s[y == 0].mean())
        return self

    def decision_function(self, X):
        return self.pls_.predict(X).ravel() - self.thr_

    def predict(self, X):
        return (self.decision_function(X) > 0).astype(int)


def make_model(name):
    pre = [("filter", ExprVarFilter()), ("scale", StandardScaler()), ("select", SelectKBest(f_classif))]
    if name == "L1-logistic":
        m = LogisticRegression(l1_ratio=1.0, solver="liblinear", class_weight="balanced", max_iter=5000, random_state=SEED)
        grid = {"select__k": K_GRID, "model__C": [0.01, 0.1, 1.0]}
    elif name == "linear SVM":
        m = LinearSVC(class_weight="balanced", max_iter=20000, random_state=SEED, dual="auto")
        grid = {"select__k": K_GRID, "model__C": [0.0005, 0.005, 0.05]}
    elif name == "Random Forest":
        m = RandomForestClassifier(n_estimators=500, class_weight="balanced", min_samples_leaf=2, n_jobs=1,
                                   random_state=SEED)
        grid = {"select__k": K_GRID}
    elif name == "PLS-DA":
        m = PLSDA()
        grid = {"select__k": K_GRID, "model__n_components": [2, 5, 10]}
    else:
        raise ValueError(name)
    return Pipeline(pre + [("model", m)]), grid


def score_of(est, X):
    if hasattr(est, "predict_proba") and not isinstance(est[-1], LinearSVC):
        return est.predict_proba(X)[:, 1]
    return est.decision_function(X)


def fold_metrics(y, s, yhat):
    return {"AUC": roc_auc_score(y, s), "balanced_acc": balanced_accuracy_score(y, yhat),
            "MCC": matthews_corrcoef(y, yhat), "F1": f1_score(y, yhat)}


def youden_threshold(x, y):
    """Threshold on a single gene maximising sensitivity + specificity - 1 on the training data (higher = positive)."""
    best, thr = -1, np.median(x)
    P, N = (y == 1).sum(), (y == 0).sum()
    for t in np.quantile(x, np.linspace(0.01, 0.99, 197)):
        j = ((x > t) & (y == 1)).sum() / P + ((x <= t) & (y == 0)).sum() / N - 1
        if j > best:
            best, thr = j, t
    return thr


# ---------------------------------------------------------------- CV runners
def outer_folds(y, seed=SEED):
    return list(RepeatedStratifiedKFold(n_splits=OUTER_SPLITS, n_repeats=OUTER_REPEATS, random_state=seed).split(np.zeros(len(y)), y))


def run_fold(model_name, X, y, tr, te, fold_id):
    pipe, grid = make_model(model_name)
    inner = StratifiedKFold(n_splits=INNER_SPLITS, shuffle=True, random_state=SEED + fold_id)
    gs = GridSearchCV(pipe, grid, scoring="roc_auc", cv=inner, n_jobs=1, refit=True)
    gs.fit(X[tr], y[tr])
    est = gs.best_estimator_
    s, yhat = score_of(est, X[te]), est.predict(X[te])
    met = fold_metrics(y[te], s, yhat)
    # genes the fitted model uses (selected by SelectKBest; for L1-logistic only those with non-zero weight)
    filt, sel = est.named_steps["filter"], est.named_steps["select"]
    genes_idx = filt.idx_[sel.get_support()]
    if model_name == "L1-logistic":
        genes_idx = genes_idx[np.abs(est.named_steps["model"].coef_.ravel()) > 0]
    return {"fold": fold_id, "metrics": met, "score": s, "pred": yhat, "test_idx": te, "best_params": gs.best_params_,
            "genes_idx": genes_idx, "eligible_idx": filt.idx_}


def run_baseline_fold(x, y, tr, te, fold_id):
    thr = youden_threshold(x[tr], y[tr])
    s = x[te]
    return {"fold": fold_id, "metrics": fold_metrics(y[te], s, (s > thr).astype(int)), "score": s,
            "pred": (s > thr).astype(int), "test_idx": te, "best_params": {"threshold": float(thr)}}


def nadeau_bengio(a, b, n_train, n_test):
    """Corrected resampled t-test (Nadeau & Bengio 2003) on paired per-fold scores a - b; returns (mean diff, t, p)."""
    d = np.asarray(a) - np.asarray(b)
    J = len(d)
    var = d.var(ddof=1)
    if var == 0:
        return float(d.mean()), np.nan, np.nan
    t = d.mean() / np.sqrt((1 / J + n_test / n_train) * var)
    p = 2 * stats.t.sf(abs(t), df=J - 1)
    return float(d.mean()), float(t), float(p)


def pooled_scores(res, n):
    """Per-repeat pooled out-of-fold scores, shape (repeats, n)."""
    S = np.full((OUTER_REPEATS, n), np.nan)
    for r in res:
        S[r["fold"] // OUTER_SPLITS, r["test_idx"]] = r["score"]
    return S


def fast_auc(y, s):
    """ROC-AUC via the Mann-Whitney rank formula (ties get mid-ranks); equals roc_auc_score, much faster in the bootstrap."""
    r = stats.rankdata(s)
    n1 = y.sum(); n0 = len(y) - n1
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def pooled_auc(S, y, idx=None):
    idx = np.arange(len(y)) if idx is None else idx
    return float(np.mean([fast_auc(y[idx], S[k, idx]) for k in range(S.shape[0])]))


def boot_ci(S, y, rng, S2=None):
    """95% bootstrap CI of the pooled OOF AUC (mean over repeats); with S2, CI and p of AUC(S) - AUC(S2)."""
    n = len(y)
    vals = []
    for _ in range(N_BOOT):
        idx = rng.integers(0, n, n)
        if y[idx].min() == y[idx].max():
            continue
        a = pooled_auc(S, y, idx)
        vals.append(a - pooled_auc(S2, y, idx) if S2 is not None else a)
    vals = np.array(vals)
    lo, hi = np.percentile(vals, [2.5, 97.5])
    if S2 is None:
        return float(lo), float(hi)
    p = 2 * min((vals <= 0).mean(), (vals >= 0).mean())
    return float(lo), float(hi), float(max(p, 1 / len(vals)))


# ---------------------------------------------------------------- data
t0 = time.time()
rng = np.random.default_rng(SEED)
lab = load_labels()
Xdf = pd.read_parquet(PROC / "tcga_hiseq_tumour_log2.parquet").loc[lab.index]
genes = Xdf.columns.to_numpy()
amp = pd.read_csv(PROC / "erbb2_amplicon_genes.tsv", sep="\t")
pam50 = set(pd.read_csv(PROC / "pam50_centroids_hgnc.tsv", sep="\t", index_col=0).index)
receptor = {"ESR1", "PGR", "ERBB2"}
amp15 = set(amp.loc[amp.within_1500kb, "hgnc"].dropna())
amp20 = set(amp.loc[amp.within_2000kb, "hgnc"].dropna())
ABL = {"a_all_genes": set(), "b_no_receptor_amplicon1.5Mb_PAM50": receptor | amp15 | pam50,
       "c_no_receptor_amplicon2Mb_PAM50": receptor | amp20 | pam50}
abl_removed = {k: int(np.isin(genes, list(v)).sum()) for k, v in ABL.items()}
print("genes removed per ablation (present in HiSeqV2):", abl_removed, flush=True)

TARGETS = {"ER": ("ER", "ESR1", 1), "HER2": ("her2", "ERBB2", 1)}   # (label column, baseline gene, sign)
MODELS = ["L1-logistic", "linear SVM", "Random Forest", "PLS-DA"]
CONFIGS = [("ER", "a_all_genes"), ("ER", "b_no_receptor_amplicon1.5Mb_PAM50"),
           ("HER2", "a_all_genes"), ("HER2", "b_no_receptor_amplicon1.5Mb_PAM50"), ("HER2", "c_no_receptor_amplicon2Mb_PAM50")]
if os.environ.get("SMOKE"):
    CONFIGS, MODELS, OUT = [("ER", "a_all_genes"), ("HER2", "b_no_receptor_amplicon1.5Mb_PAM50")], ["L1-logistic", "PLS-DA"], OUT / "smoke"
    OUT.mkdir(exist_ok=True)


def target_data(target, labels=lab):
    col, gene, sign = TARGETS[target]
    m = labels[col].isin(["Positive", "Negative"])
    ids = labels.index[m]
    y = (labels.loc[ids, col] == "Positive").astype(int).to_numpy()
    return ids, y, gene, sign


fold_rows, summary_rows, oof_rows, feat_rows, params_rows = [], [], [], [], []
store = {}


def evaluate(target, abl, models, labels=lab, tag="main"):
    ids, y, gene, sign = target_data(target, labels)
    keep = ~np.isin(genes, list(ABL[abl]))
    X = Xdf.loc[ids].to_numpy()[:, keep]
    g_used = genes[keep]
    folds = outer_folds(y)
    n_test = np.mean([len(te) for _, te in folds]); n_train = len(y) - n_test
    xb = sign * Xdf.loc[ids, gene].to_numpy()
    base = [run_baseline_fold(xb, y, tr, te, i) for i, (tr, te) in enumerate(folds)]
    Sb = pooled_scores(base, len(y))
    out = {"baseline": (base, Sb)}
    for r in base:
        fold_rows.append({"tag": tag, "target": target, "ablation": abl, "model": f"single gene ({gene})", "fold": r["fold"], **r["metrics"]})
    for mname in models:
        t1 = time.time()
        res = Parallel(n_jobs=N_JOBS)(delayed(run_fold)(mname, X, y, tr, te, i) for i, (tr, te) in enumerate(folds))
        S = pooled_scores(res, len(y))
        out[mname] = (res, S)
        for r in res:
            fold_rows.append({"tag": tag, "target": target, "ablation": abl, "model": mname, "fold": r["fold"], **r["metrics"]})
            params_rows.append({"tag": tag, "target": target, "ablation": abl, "model": mname, "fold": r["fold"], **r["best_params"]})
        # how often each gene is used across the 15 outer fits, and how often it was eligible
        cnt = pd.Series(np.concatenate([g_used[r["genes_idx"]] for r in res])).value_counts()
        elig = pd.Series(np.concatenate([g_used[r["eligible_idx"]] for r in res])).value_counts()
        for g_, c in cnt.items():
            feat_rows.append({"tag": tag, "target": target, "ablation": abl, "model": mname, "gene": g_,
                              "folds_selected": int(c), "folds_eligible": int(elig.get(g_, 0))})
        print(f"  {tag} {target} {abl} {mname}: AUC {np.mean([r['metrics']['AUC'] for r in res]):.3f} "
              f"({time.time() - t1:.0f}s)", flush=True)
    # summaries
    for mname, (res, S) in out.items():
        fm = pd.DataFrame([r["metrics"] for r in res])
        lo, hi = boot_ci(S, y, rng)
        row = {"tag": tag, "target": target, "ablation": abl, "model": mname if mname != "baseline" else f"single gene ({gene})",
               "n": len(y), "n_pos": int(y.sum()), "AUC_fold_mean": fm.AUC.mean(), "AUC_fold_sd": fm.AUC.std(),
               "AUC_pooled": pooled_auc(S, y), "AUC_CI_lo": lo, "AUC_CI_hi": hi,
               "balanced_acc": fm.balanced_acc.mean(), "MCC": fm.MCC.mean(), "F1": fm.F1.mean()}
        if mname != "baseline":
            d, t, p = nadeau_bengio(fm.AUC, [r["metrics"]["AUC"] for r in base], n_train, n_test)
            dl, dh, pb = boot_ci(S, y, rng, Sb)
            row.update({"dAUC_vs_baseline": d, "NB_t": t, "NB_p": p, "boot_dAUC_lo": dl, "boot_dAUC_hi": dh, "boot_p": pb})
            dmcc = nadeau_bengio(fm.MCC, [r["metrics"]["MCC"] for r in base], n_train, n_test)
            row.update({"dMCC_vs_baseline": dmcc[0], "NB_p_MCC": dmcc[2]})
        summary_rows.append(row)
        sub = res
        for r in sub:
            for i, s, pr in zip(r["test_idx"], r["score"], r["pred"]):
                oof_rows.append((tag, target, abl, row["model"], r["fold"], ids[i], int(y[i]), float(s), int(pr)))
    store[(tag, target, abl)] = (ids, y, out)
    return out


# ---------------------------------------------------------------- main grid
for target, abl in CONFIGS:
    print(f"== {target} / {abl}", flush=True)
    evaluate(target, abl, MODELS)

# ablation b (and c) vs a, same model, same folds
abl_tests = []
for target in TARGETS:
    if ("main", target, "a_all_genes") not in store:
        continue
    ids, y, out_a = store[("main", target, "a_all_genes")]
    folds = outer_folds(y)
    n_test = np.mean([len(te) for _, te in folds]); n_train = len(y) - n_test
    for abl in [a for t, a in CONFIGS if t == target and a != "a_all_genes"]:
        if ("main", target, abl) not in store:
            continue
        out_b = store[("main", target, abl)][2]
        for m in MODELS:
            a_auc = [r["metrics"]["AUC"] for r in out_a[m][0]]
            b_auc = [r["metrics"]["AUC"] for r in out_b[m][0]]
            d, t, p = nadeau_bengio(b_auc, a_auc, n_train, n_test)
            dl, dh, pb = boot_ci(out_b[m][1], y, rng, out_a[m][1])
            abl_tests.append({"target": target, "ablation": abl, "model": m, "dAUC_b_minus_a": d, "NB_t": t, "NB_p": p,
                              "boot_lo": dl, "boot_hi": dh, "boot_p": pb})
pd.DataFrame(abl_tests).to_csv(OUT / "ablation_tests.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- controls
controls = {}
for target in TARGETS:
    ids, y, gene, _ = target_data(target)
    X = Xdf.loc[ids].to_numpy()
    yp = np.random.default_rng(SEED + 7).permutation(y)
    folds = list(StratifiedKFold(5, shuffle=True, random_state=SEED).split(X, yp))
    res = Parallel(n_jobs=N_JOBS)(delayed(run_fold)("L1-logistic", X, yp, tr, te, i) for i, (tr, te) in enumerate(folds))
    fm = pd.DataFrame([r["metrics"] for r in res])
    controls[target] = {"shuffled_labels_L1_logistic": {"AUC": fm.AUC.mean(), "MCC": fm.MCC.mean(), "MCC_sd": fm.MCC.std(),
                                                         "balanced_acc": fm.balanced_acc.mean()},
                        "majority_class": {"accuracy": float(max(y.mean(), 1 - y.mean())), "balanced_acc": 0.5, "MCC": 0.0,
                                           "AUC": 0.5}}
    print("control", target, controls[target], flush=True)

# ---------------------------------------------------------------- HER2 label-confidence sensitivity
sens = {}
for filt in (["drop_low_pos"] if os.environ.get("SMOKE") else ["drop_low_pos", "drop_low_all"]):
    labf = load_labels(her2_filter=filt)
    for abl in (["a_all_genes"] if os.environ.get("SMOKE") else ["a_all_genes", "b_no_receptor_amplicon1.5Mb_PAM50"]):
        print(f"== sensitivity {filt} HER2 / {abl}", flush=True)
        evaluate("HER2", abl, ["L1-logistic", "linear SVM"], labels=labf, tag=f"sens_{filt}")

# ER-/HER2- vs rest (for METABRIC transfer), L1-logistic, all genes
lab_t = lab.copy()
lab_t["tn"] = np.where(lab_t.erher2.eq("ERneg_HER2neg"), "Positive",
                       np.where(lab_t.ER.notna() & lab_t.her2.notna(), "Negative", None))
TARGETS["ERneg_HER2neg"] = ("tn", "ESR1", -1)   # baseline: low ESR1
print("== ER-/HER2- vs rest", flush=True)
evaluate("ERneg_HER2neg", "a_all_genes", ["L1-logistic"], labels=lab_t, tag="transfer_target")

# ---------------------------------------------------------------- secondary benchmark: PAM50 (circular)
pam = {}
for ncls in ([4] if os.environ.get("SMOKE") else [4, 5]):
    for abl in ["a_all_genes", "no_PAM50_genes"]:
        sub = lab[lab.PAM50_pancan.notna() & (lab.PAM50_pancan != "Normal" if ncls == 4 else lab.PAM50_pancan.notna())]
        yc = sub.PAM50_pancan.to_numpy()
        keep = ~np.isin(genes, list(pam50)) if abl == "no_PAM50_genes" else np.ones(len(genes), bool)
        X = Xdf.loc[sub.index].to_numpy()[:, keep]
        folds = list(StratifiedKFold(5, shuffle=True, random_state=SEED).split(X, yc))

        def pam_fold(tr, te):
            pipe = Pipeline([("filter", ExprVarFilter()), ("scale", StandardScaler()), ("select", SelectKBest(f_classif, k=1000)),
                             ("model", LogisticRegression(l1_ratio=1.0, solver="saga", C=0.1, max_iter=3000,
                                                          class_weight="balanced", random_state=SEED))])
            pipe.fit(X[tr], yc[tr])
            p = pipe.predict(X[te])
            return {"accuracy": accuracy_score(yc[te], p), "macro_F1": f1_score(yc[te], p, average="macro"),
                    "kappa": cohen_kappa_score(yc[te], p)}
        res = Parallel(n_jobs=N_JOBS)(delayed(pam_fold)(tr, te) for tr, te in folds)
        fm = pd.DataFrame(res)
        pam[f"{ncls}class_{abl}"] = {"n": len(yc), **fm.mean().to_dict(), "accuracy_sd": fm.accuracy.std()}
        print("PAM50", ncls, abl, pam[f"{ncls}class_{abl}"], flush=True)

# ---------------------------------------------------------------- write
pd.DataFrame(fold_rows).to_csv(OUT / "fold_metrics.tsv", sep="\t", index=False)
summ = pd.DataFrame(summary_rows)
summ.to_csv(OUT / "summary.tsv", sep="\t", index=False)
pd.DataFrame(params_rows).to_csv(OUT / "best_params_per_fold.tsv", sep="\t", index=False)
fe = pd.DataFrame(feat_rows)
fe.to_csv(OUT / "feature_selection_frequency.tsv", sep="\t", index=False)
oof = pd.DataFrame(oof_rows, columns=["tag", "target", "ablation", "model", "fold", "sample", "y", "score", "pred"])
oof.to_parquet(OUT / "oof_predictions.parquet")
report = {"genes_removed_per_ablation": abl_removed, "controls": controls, "pam50_secondary_benchmark": pam,
          "n_per_target": {t: {"n": int(len(target_data(t, lab_t)[1])), "pos": int(target_data(t, lab_t)[1].sum())}
                           for t in TARGETS},
          "runtime_min": round((time.time() - t0) / 60, 1)}
dump(report, OUT / "classify_report.json")
write_params3(OUT, {"outer": f"RepeatedStratifiedKFold({OUTER_SPLITS}x{OUTER_REPEATS}, seed {SEED})",
                    "inner": f"StratifiedKFold({INNER_SPLITS}), GridSearchCV scoring roc_auc",
                    "pipeline": "ExprVarFilter(mean log2>=%s, top %d variance) -> StandardScaler -> SelectKBest(f_classif)" % (MIN_MEAN, N_VAR),
                    "k_grid": K_GRID, "models": MODELS, "class_weight": "balanced",
                    "L1_C": [0.01, 0.1, 1.0], "SVM_C": [0.0005, 0.005, 0.05], "RF": "500 trees, min_samples_leaf 2, tune k only",
                    "PLSDA_components": [2, 5, 10], "PLSDA_threshold": "midpoint of training-fold class mean scores",
                    "baseline": "single gene, Youden threshold on the training fold",
                    "tests": "Nadeau-Bengio corrected resampled t on 15 fold AUCs; paired bootstrap (2000) on pooled OOF AUC",
                    "ablations": {k: sorted(v) if len(v) < 400 else len(v) for k, v in ABL.items()}})
print(f"done in {(time.time() - t0) / 60:.1f} min")
