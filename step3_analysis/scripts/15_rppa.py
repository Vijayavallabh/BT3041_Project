"""Step 3.5: protein-level validation with TCGA RPPA (same patients, another assay; plan section 3F.2-3).

1. Mann-Whitney TNBC vs HR+/HER2- for every RPPA antibody (131), BH over antibodies; effect = difference in
   medians (TNBC - HR+) and rank-biserial r.
2. Concordance: Spearman between the DE shrunk log2FC (main run) and the RPPA median difference over genes with a
   total-protein antibody (phospho/cleaved forms and isoform-ambiguous antibodies left out; one value per gene).
   Null: permute the TNBC/HR+ labels within the RPPA set 1,000 times, recompute the RPPA effects, correlate with the
   fixed DE log2FC. The same for HER2+ vs HR+/HER2- (DE run `her2`).
3. Error analysis: tumours whose out-of-fold ER call (L1-logistic, all genes, majority over the 3 repeats)
   disagrees with IHC ER. Their RPPA ER-alpha, ESR1 mRNA and PAM50 separate likely IHC errors from ER-low biology.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.common3 import PROC, RES3, SEED, load_labels, bh, write_params3, dump  # noqa: E402

OUT = RES3 / "rppa"
OUT.mkdir(parents=True, exist_ok=True)
N_PERM = 1000
rng = np.random.default_rng(SEED)
lab = load_labels()
R = pd.read_parquet(PROC / "tcga_rppa.parquet")
R = R.loc[R.index.intersection(lab.index)]
amap = pd.read_csv(PROC / "rppa_antibody_map.tsv", sep="\t").set_index("antibody")
X = pd.read_parquet(PROC / "tcga_hiseq_tumour_log2.parquet")
report = {"rppa_tumours_supervised": len(R)}


def mw_table(a_ids, b_ids, label):
    rows = []
    for ab in R.columns:
        a, b = R.loc[a_ids, ab].dropna(), R.loc[b_ids, ab].dropna()
        u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        rows.append({"contrast": label, "antibody": ab, "gene": amap.gene.get(ab), "total_protein": bool(amap.total_protein.get(ab, False)),
                     "isoform_ambiguous": bool(amap.isoform_ambiguous.get(ab, False)), "n_a": len(a), "n_b": len(b),
                     "median_diff": a.median() - b.median(), "rank_biserial": 2 * u / (len(a) * len(b)) - 1, "p": p})
    t = pd.DataFrame(rows)
    t["q_BH"] = bh(t.p)
    return t


def concordance(a_ids, b_ids, de_run, label):
    de = pd.read_csv(RES3 / "de" / f"{de_run}_results.tsv", sep="\t", index_col=0).dropna(subset=["symbol"])
    de = de.drop_duplicates("symbol").set_index("symbol")
    abs_ = amap[amap.total_protein & ~amap.isoform_ambiguous & amap.gene.isin(de.index)].index.intersection(R.columns)
    genes = amap.loc[abs_, "gene"]
    sub = R.loc[list(a_ids) + list(b_ids), abs_]
    is_a = np.r_[np.ones(len(a_ids), bool), np.zeros(len(b_ids), bool)]

    def effects(mask):
        e = sub[mask].median() - sub[~mask].median()
        return e.groupby(genes).mean()           # one value per gene
    eff = effects(is_a)
    lfc = de.loc[eff.index, "log2FC_shrunk"]
    rho = stats.spearmanr(lfc, eff).statistic
    null = np.array([stats.spearmanr(lfc, effects(rng.permutation(is_a)).loc[eff.index]).statistic for _ in range(N_PERM)])
    p_perm = (1 + (np.abs(null) >= abs(rho)).sum()) / (N_PERM + 1)
    tab = pd.DataFrame({"gene": eff.index, "rppa_median_diff": eff.values, "de_log2FC_shrunk": lfc.values,
                        "de_DE": de.loc[eff.index, "DE"].values, "de_padj_lfc1": de.loc[eff.index, "padj_lfc1"].values})
    tab.to_csv(OUT / f"concordance_{label}.tsv", sep="\t", index=False)
    de_genes = tab[tab.de_DE]
    return {"genes": len(eff), "antibodies": len(abs_), "spearman": rho, "perm_p": p_perm, "null_95pct_abs": float(np.percentile(np.abs(null), 95)),
            "sign_agreement_all": float((np.sign(tab.rppa_median_diff) == np.sign(tab.de_log2FC_shrunk)).mean()),
            "DE_genes_with_antibody": len(de_genes),
            "sign_agreement_DE_genes": float((np.sign(de_genes.rppa_median_diff) == np.sign(de_genes.de_log2FC_shrunk)).mean()) if len(de_genes) else None}


tn = lab.index[lab.group.eq("TNBC")].intersection(R.index)
hr = lab.index[lab.group.eq("HR+/HER2-")].intersection(R.index)
h2 = lab.index[lab.group.eq("HER2+")].intersection(R.index)
report["n"] = {"TNBC": len(tn), "HR+/HER2-": len(hr), "HER2+": len(h2)}
mw = pd.concat([mw_table(tn, hr, "TNBC vs HR+/HER2-"), mw_table(h2, hr, "HER2+ vs HR+/HER2-")])
mw.to_csv(OUT / "rppa_mannwhitney.tsv", sep="\t", index=False)
m1 = mw[mw.contrast == "TNBC vs HR+/HER2-"]
report["mw_TNBC_vs_HR"] = {"antibodies": len(m1), "q_lt_0.05": int((m1.q_BH < 0.05).sum()),
                           "top_higher_in_TNBC": m1[m1.q_BH < 0.05].sort_values("rank_biserial", ascending=False).antibody.head(8).tolist(),
                           "top_lower_in_TNBC": m1[m1.q_BH < 0.05].sort_values("rank_biserial").antibody.head(8).tolist(),
                           "selected": m1.set_index("antibody").loc[[a for a in ["ERALPHA", "PR", "GATA3", "AR", "BCL2", "CYCLINB1", "EGFR",
                                                                                 "CAVEOLIN1", "HER2", "INPP4B"] if a in m1.antibody.values],
                                                                    ["median_diff", "rank_biserial", "q_BH"]].round(4).to_dict("index")}
m2 = mw[mw.contrast == "HER2+ vs HR+/HER2-"].set_index("antibody")
report["mw_HER2_vs_HR"] = {"q_lt_0.05": int((m2.q_BH < 0.05).sum()), "HER2": m2.loc["HER2", ["median_diff", "rank_biserial", "q_BH"]].to_dict(),
                           "HER2PY1248": m2.loc["HER2PY1248", ["median_diff", "rank_biserial", "q_BH"]].to_dict()}
report["concordance_TNBC_vs_HR"] = concordance(tn, hr, "main", "TNBC_vs_HR")
report["concordance_HER2_vs_HR"] = concordance(h2, hr, "her2", "HER2_vs_HR")
print(report["concordance_TNBC_vs_HR"], report["concordance_HER2_vs_HR"], flush=True)

# ------------------------------------------------ error analysis of the ER classifier
oof = pd.read_parquet(RES3 / "classify" / "oof_predictions.parquet")
o = oof[(oof.tag == "main") & (oof.target == "ER") & (oof.ablation == "a_all_genes") & (oof.model == "L1-logistic")]
call = o.groupby("sample").agg(y=("y", "first"), pred_votes=("pred", "sum"), n=("pred", "size"), score=("score", "mean"))
call["pred"] = (call.pred_votes * 2 > call.n).astype(int)
call["error"] = np.where((call.y == 1) & (call.pred == 0), "IHC ER+ called ER-",
                         np.where((call.y == 0) & (call.pred == 1), "IHC ER- called ER+",
                                  np.where(call.y == 1, "correct ER+", "correct ER-")))
call["ESR1_mRNA"] = X.loc[call.index, "ESR1"]
call["RPPA_ERalpha"] = R.ERALPHA.reindex(call.index)
call["PR_IHC"] = lab.loc[call.index, "PR"]
call["PAM50"] = lab.loc[call.index, "PAM50_pancan"]
call["group"] = lab.loc[call.index, "group"]
call.to_csv(OUT / "er_error_analysis.tsv", sep="\t")
ea = {}
for k, g in call.groupby("error"):
    ea[k] = {"n": len(g), "with_RPPA": int(g.RPPA_ERalpha.notna().sum()), "median_ESR1": g.ESR1_mRNA.median(),
             "median_RPPA_ERalpha": g.RPPA_ERalpha.median(), "PAM50": g.PAM50.value_counts().to_dict(),
             "PR_positive": int((g.PR_IHC == "Positive").sum())}
# do misclassified IHC ER- tumours look like ER+ at the protein level? compare with correct ER+ and correct ER-
fp = call[call.error == "IHC ER- called ER+"].RPPA_ERalpha.dropna()
fn = call[call.error == "IHC ER+ called ER-"].RPPA_ERalpha.dropna()
cp = call[call.error == "correct ER+"].RPPA_ERalpha.dropna()
cn = call[call.error == "correct ER-"].RPPA_ERalpha.dropna()
thr_hi = cn.quantile(0.95)
ea["tests"] = {"ERneg_called_pos_vs_correct_neg_p": stats.mannwhitneyu(fp, cn).pvalue if len(fp) > 2 else None,
               "ERpos_called_neg_vs_correct_pos_p": stats.mannwhitneyu(fn, cp).pvalue if len(fn) > 2 else None,
               "ERalpha_95pct_of_correct_ERneg": thr_hi,
               "ERneg_called_pos_above_that": f"{int((fp > thr_hi).sum())}/{len(fp)}",
               "ERpos_called_neg_below_median_of_correct_ERneg": f"{int((fn < cn.median()).sum())}/{len(fn)}",
               "ERpos_called_neg_below_5pct_of_correct_ERpos": f"{int((fn < cp.quantile(0.05)).sum())}/{len(fn)}"}
report["er_error_analysis"] = ea
print(ea, flush=True)
dump(report, OUT / "rppa_report.json")
write_params3(OUT, {"test": "Mann-Whitney two-sided, BH over 131 antibodies", "concordance": "Spearman DE shrunk LFC vs RPPA median "
                    "difference; total-protein, non-ambiguous antibodies, mean per gene", "permutations": N_PERM,
                    "error_analysis": "L1-logistic all genes, majority vote over 3 repeats of outer CV"})
