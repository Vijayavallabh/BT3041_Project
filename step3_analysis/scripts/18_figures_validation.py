"""Step 3.8: figures for classification, protein validation, METABRIC and survival (report + 16:9 slides 8 and 9).

fig13_classification    ROC (pooled out-of-fold) for ER and HER2, and AUC with 95% CI for every model x ablation
fig14_rppa              DE vs RPPA concordance (with permutation null) and the ER error analysis
fig15_metabric          transfer AUC (TCGA CV vs METABRIC) and DE replication (effect-size scatter)
fig16_survival          Kaplan-Meier: TCGA PFI by clinical group; METABRIC breast-cancer death (10 y) by IHC group and
                        by the frozen classifier's predicted ER-/HER2- group
slide8_classification, slide9_validation, slideS1_survival (backup)   16:9 composites
Run after scripts 12-16.
"""
import sys, json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve
from lifelines import KaplanMeierFitter
from adjustText import adjust_text

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.common3 import (RES3, PROC, GROUP_COL, PAM50_COL, INK, INK2, MUTED, GRID, UP_COL, DOWN_COL, SLIDE, style,  # noqa: E402
                         slide_style, save3, load_labels)

CL, RP, TR, SV = RES3 / "classify", RES3 / "rppa", RES3 / "transfer", RES3 / "survival"
summ = pd.read_csv(CL / "summary.tsv", sep="\t")
oof = pd.read_parquet(CL / "oof_predictions.parquet")
crep = json.load(open(CL / "classify_report.json"))
abt = pd.read_csv(CL / "ablation_tests.tsv", sep="\t")
rrep = json.load(open(RP / "rppa_report.json"))
mrep = json.load(open(TR / "metabric_report.json"))
srep = json.load(open(SV / "survival_report.json"))
MODELS = ["L1-logistic", "linear SVM", "Random Forest", "PLS-DA"]
ABL = {"a_all_genes": "all genes", "b_no_receptor_amplicon1.5Mb_PAM50": "− receptor, amplicon ±1.5 Mb, PAM50",
       "c_no_receptor_amplicon2Mb_PAM50": "− receptor, amplicon ±2 Mb, PAM50"}
SHORT = {"L1-logistic": "L1-LR", "linear SVM": "SVM", "Random Forest": "RF", "PLS-DA": "PLS-DA"}
ABL_SHORT = {"a_all_genes": "all genes", "b_no_receptor_amplicon1.5Mb_PAM50": "ablated", "c_no_receptor_amplicon2Mb_PAM50": "ablated 2 Mb"}
COL_A, COL_B, COL_BASE = "#3d3c38", "#e07b39", "#2a78d6"
main = summ[summ.tag == "main"]


def best_model(target, abl="a_all_genes"):
    d = main[(main.target == target) & (main.ablation == abl) & main.model.isin(MODELS)]
    return d.sort_values("AUC_pooled", ascending=False).iloc[0]


def roc_panel(ax, target, fs=9, abl_b="b_no_receptor_amplicon1.5Mb_PAM50"):
    gene = {"ER": "ESR1", "HER2": "ERBB2"}[target]
    ba, bb = best_model(target), best_model(target, abl_b)
    base = main[(main.target == target) & (main.ablation == "a_all_genes") & main.model.str.startswith("single")].iloc[0]
    curves = [(base.model, "a_all_genes", COL_BASE, f"{gene} alone", base), (ba.model, "a_all_genes", COL_A, f"{SHORT[ba.model]}, all genes", ba),
              (bb.model, abl_b, COL_B, f"{SHORT[bb.model]}, ablated", bb)]
    for model, abl, col, lab, row in curves:
        o = oof[(oof.tag == "main") & (oof.target == target) & (oof.ablation == abl) & (oof.model == model) & (oof.fold < 5)]
        f, t, _ = roc_curve(o.y, o.score)
        ax.plot(f, t, color=col, lw=2, label=f"{lab}: {row.AUC_pooled:.3f}")
    ax.plot([0, 1], [0, 1], color=GRID, lw=1, zorder=0)
    ax.set_xlabel("1 − specificity"); ax.set_ylabel("sensitivity")
    ax.legend(loc="lower right", fontsize=fs, handlelength=1.2, title="AUC (pooled out-of-fold)", title_fontsize=fs)
    ax.set_xlim(-0.01, 1); ax.set_ylim(0, 1.01)


def auc_dots(ax, fs=9):
    rows = []
    for target in ["ER", "HER2"]:
        for abl in ["a_all_genes", "b_no_receptor_amplicon1.5Mb_PAM50"] + (["c_no_receptor_amplicon2Mb_PAM50"] if target == "HER2" else []):
            for m in MODELS:
                r = main[(main.target == target) & (main.ablation == abl) & (main.model == m)].iloc[0]
                rows.append((target, abl, m, r))
    y = 0
    ticks, labels = [], []
    for target in ["ER", "HER2"]:
        base = main[(main.target == target) & main.model.str.startswith("single")].iloc[0]
        sub = [r for r in rows if r[0] == target]
        y0 = y
        for t, abl, m, r in sub:
            col = COL_A if abl == "a_all_genes" else COL_B
            ax.errorbar(r.AUC_pooled, y, xerr=[[r.AUC_pooled - r.AUC_CI_lo], [r.AUC_CI_hi - r.AUC_pooled]], fmt="o", color=col,
                        ms=5, lw=1.4, capsize=0, mfc=col if abl != "c_no_receptor_amplicon2Mb_PAM50" else "white")
            ticks.append(y); labels.append(f"{SHORT[m]} · {ABL_SHORT[abl]}")
            y += 1
        ax.axvspan(base.AUC_CI_lo, base.AUC_CI_hi, ymin=0, ymax=1, color=COL_BASE, alpha=0.0)
        ax.plot([base.AUC_pooled] * 2, [y0 - 0.5, y - 0.5], color=COL_BASE, lw=2)
        ax.fill_betweenx([y0 - 0.5, y - 0.5], base.AUC_CI_lo, base.AUC_CI_hi, color=COL_BASE, alpha=0.12, lw=0)
        ax.text(base.AUC_CI_lo - 0.003, y0 - 0.75, f"{target}: {base.model.replace('single gene ', '').strip('()')} alone {base.AUC_pooled:.3f}",
                color=COL_BASE, fontsize=fs, ha="right", va="bottom", fontweight="bold")
        y += 1.4
    ax.set_yticks(ticks); ax.set_yticklabels(labels, fontsize=fs - 0.5)
    ax.invert_yaxis()
    ax.set_xlabel("ROC-AUC, pooled out-of-fold (95% CI)")


def concordance_panel(ax, fs=9):
    t = pd.read_csv(RP / "concordance_TNBC_vs_HR.tsv", sep="\t")
    c = rrep["concordance_TNBC_vs_HR"]
    col = np.where(t.de_DE, np.where(t.de_log2FC_shrunk > 0, UP_COL, DOWN_COL), "#b5b4ab")
    ax.scatter(t.de_log2FC_shrunk, t.rppa_median_diff, s=22, c=col, lw=0)
    ax.axhline(0, color=GRID, lw=0.8); ax.axvline(0, color=GRID, lw=0.8)
    texts = []
    for g in ["ESR1", "PGR", "GATA3", "AR", "BCL2", "INPP4B", "CCNB1", "EGFR", "CAV1", "ERBB2", "CDH1", "KIT", "ASNS"]:
        if g in set(t.gene):
            r = t[t.gene == g].iloc[0]
            texts.append(ax.text(r.de_log2FC_shrunk, r.rppa_median_diff, g, fontsize=fs - 1, fontstyle="italic"))
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.5))
    ax.set_xlabel("mRNA: shrunk log2FC (DESeq2)"); ax.set_ylabel("protein: RPPA median difference")
    ax.text(0.97, 0.03, f"Spearman ρ = {c['spearman']:.2f}\npermutation p = {c['perm_p']:.3f}\n{c['genes']} genes with a\ntotal-protein antibody",
            transform=ax.transAxes, va="bottom", ha="right", fontsize=fs, color=INK2)


def error_panel(ax, fs=9):
    e = pd.read_csv(RP / "er_error_analysis.tsv", sep="\t", index_col=0).dropna(subset=["RPPA_ERalpha"])
    order = ["correct ER-", "IHC ER+ called ER-", "IHC ER- called ER+", "correct ER+"]
    rng = np.random.default_rng(1)
    for i, k in enumerate(order):
        d = e[e.error == k]
        cols = [PAM50_COL.get(p, "#b5b4ab") for p in d.PAM50]
        ax.scatter(i + rng.uniform(-0.28, 0.28, len(d)), d.RPPA_ERalpha, s=10 if "correct" in k else 26, c=cols, lw=0,
                   alpha=0.45 if "correct" in k else 0.95)
        ax.plot([i - 0.32, i + 0.32], [d.RPPA_ERalpha.median()] * 2, color=INK, lw=2)
        ax.text(i, ax.get_ylim()[1] if False else e.RPPA_ERalpha.max() + 0.25, f"n={len(d)}", ha="center", fontsize=fs - 1, color=INK2)
    ax.set_xticks(range(4))
    ax.set_xticklabels(["ER−\ncorrect", "ER+,\ncalled\nER−", "ER−,\ncalled\nER+", "ER+\ncorrect"], fontsize=fs - 1)
    ax.set_xlabel("IHC ER status vs out-of-fold call", fontsize=fs)
    ax.set_ylabel("RPPA ERα (protein)")
    from matplotlib.lines import Line2D
    h = [Line2D([], [], marker="o", ls="", color=PAM50_COL[p], label=p) for p in ["LumA", "LumB", "Her2", "Basal", "Normal"]]
    ax.legend(handles=h, fontsize=fs - 2, ncol=1, loc="upper left", bbox_to_anchor=(0, 0.93), handletextpad=0.1, borderaxespad=0.1,
              title="PAM50", title_fontsize=fs - 2, frameon=False)


def transfer_panel(ax, fs=9):
    t = pd.DataFrame(mrep["transfer"])
    t = t[t.scaling == "z_within_cohort"].set_index("model")
    order = ["ER", "ER_ablated", "HER2", "HER2_ablated", "ERneg_HER2neg"]
    names = {"ER": "ER (IHC)", "ER_ablated": "ER, ablated", "HER2": "HER2 (SNP6 gain)", "HER2_ablated": "HER2, ablated",
             "ERneg_HER2neg": "ER−/HER2− vs rest"}
    base = mrep["single_gene_baselines"]
    for i, m in enumerate(order):
        r = t.loc[m]
        ax.plot([r.TCGA_CV_AUC, r.AUC], [i, i], color=GRID, lw=2, zorder=0)
        ax.scatter([r.TCGA_CV_AUC], [i], s=50, facecolors="white", edgecolors=COL_A, lw=1.5, zorder=3)
        ax.errorbar(r.AUC, i, xerr=[[r.AUC - r.AUC_CI[0]], [r.AUC_CI[1] - r.AUC]], fmt="o", color=COL_B, ms=7, capsize=0, zorder=4)
        bk = m.replace("_ablated", "")
        ax.scatter([base[bk]["AUC"]], [i], marker="|", s=160, color=COL_BASE, lw=2.5, zorder=5)
    ax.set_yticks(range(len(order))); ax.set_yticklabels([names[m] for m in order], fontsize=fs)
    ax.invert_yaxis()
    ax.set_xlabel("ROC-AUC")
    from matplotlib.lines import Line2D
    h = [Line2D([], [], marker="o", ls="", mfc="white", mec=COL_A, label="TCGA (CV)"), Line2D([], [], marker="o", ls="", color=COL_B, label="METABRIC (frozen model)"),
         Line2D([], [], marker="|", ls="", color=COL_BASE, ms=12, mew=2.5, label="METABRIC, receptor gene alone")]
    ax.legend(handles=h, fontsize=fs - 1.5, loc="upper left", frameon=False, borderaxespad=0.1)


def km_panel(ax, df, t, e, g, order, colors, fs=9, xmax=None, xlabel="years"):
    for lvl in order:
        d = df[df[g] == lvl]
        if len(d) == 0:
            continue
        k = KaplanMeierFitter().fit(d[t], d[e], label=f"{lvl} ({len(d)}, {int(d[e].sum())} ev.)")
        k.plot_survival_function(ax=ax, ci_show=False, color=colors[lvl], lw=2)
    ax.set_xlabel(xlabel); ax.set_ylabel("event-free probability")
    ax.legend(fontsize=fs - 1, loc="lower left", frameon=False)
    if xmax:
        ax.set_xlim(0, xmax)
    ax.set_ylim(0, 1.02)


lab = load_labels()
ml = pd.read_csv(PROC / "metabric_labels.tsv", sep="\t", index_col=0)
pr = pd.read_csv(TR / "metabric_predictions.tsv", sep="\t", index_col=0)


def mb_trunc(yrs=10):
    m = ml[ml.OS_MONTHS.notna() & ml.BC_death.notna()].join(pr[["pred_ERneg_HER2neg"]])
    t = m.OS_MONTHS / 12
    m = m.assign(years=np.minimum(t, yrs), BCD=((m.BC_death == 1) & (t <= yrs)).astype(int))
    m["group"] = m.contrast_group.map({"ER+/HER2-notgained": "ER+/HER2 not gained", "HER2-gained": "HER2 gained",
                                       "ER-/HER2-notgained": "ER−/HER2 not gained"})
    m["pred"] = m.pred_ERneg_HER2neg.map({1: "predicted ER−/HER2−", 0: "predicted other"})
    return m


MB_COL = {"ER+/HER2 not gained": GROUP_COL["HR+/HER2-"], "HER2 gained": GROUP_COL["HER2+"], "ER−/HER2 not gained": GROUP_COL["TNBC"],
          "predicted ER−/HER2−": GROUP_COL["TNBC"], "predicted other": "#7d7c74"}
TC_COL = {"HR+/HER2-": GROUP_COL["HR+/HER2-"], "HER2+": GROUP_COL["HER2+"], "TNBC": GROUP_COL["TNBC"]}


def fmt_p(p):
    return f"{p:.2g}" if p >= 1e-4 else f"{p:.0e}"


# ====================================================== report figures
style()
fig, axs = plt.subplots(1, 3, figsize=(10, 3.6), gridspec_kw={"width_ratios": [1, 1, 1.25], "wspace": 0.45})
roc_panel(axs[0], "ER", fs=6.5); axs[0].set_title("IHC ER", loc="left", fontsize=9)
roc_panel(axs[1], "HER2", fs=6.5); axs[1].set_title("IHC/FISH HER2", loc="left", fontsize=9)
auc_dots(axs[2], fs=7); axs[2].set_title("All models (nested CV, 5 × 3 folds)", loc="left", fontsize=9)
save3(fig, "fig13_classification")

fig, axs = plt.subplots(1, 2, figsize=(7.5, 3.6), gridspec_kw={"wspace": 0.35})
concordance_panel(axs[0], fs=7); axs[0].set_title("TNBC vs HR+/HER2−: mRNA vs protein", loc="left", fontsize=9)
error_panel(axs[1], fs=7.5); axs[1].set_title("ER classifier errors vs ERα protein", loc="left", fontsize=9)
save3(fig, "fig14_rppa")

fig, axs = plt.subplots(1, 2, figsize=(7.5, 3.4), gridspec_kw={"wspace": 0.55, "width_ratios": [1.1, 1]})
transfer_panel(axs[0], fs=7); axs[0].set_title("Classifier transfer to METABRIC", loc="left", fontsize=9)
de = pd.read_csv(RES3 / "de" / "matched_erher2_results.tsv", sep="\t").dropna(subset=["symbol"]).drop_duplicates("symbol").set_index("symbol")
mo = pd.read_csv(TR / "metabric_de_ols.tsv", sep="\t", index_col=0)
com = de.index.intersection(mo.index)
rp = mrep["de_replication"]
axs[1].scatter(de.loc[com, "log2FC_shrunk"], mo.loc[com, "beta_ERneg"], s=2, c=np.where(de.loc[com, "DE"], "#3d3c38", "#d9d8d2"), lw=0, rasterized=True)
axs[1].axhline(0, color=GRID, lw=0.8); axs[1].axvline(0, color=GRID, lw=0.8)
axs[1].set_xlabel("TCGA RNA-seq log2FC (shrunk)"); axs[1].set_ylabel("METABRIC array difference (log2)")
axs[1].set_title("ER−/HER2− vs ER+/HER2−: DE replication", loc="left", fontsize=9)
axs[1].text(0.03, 0.97, f"ρ = {rp['spearman_effect_all_shared']:.2f} (all genes)\n{rp['replicated_same_sign_q05']:,}/{rp['TCGA_DEGs_on_array']:,} TCGA DEGs "
            f"({rp['fraction']:.0%})\nreplicate, same sign, q < 0.05", transform=axs[1].transAxes, va="top", fontsize=7, color=INK2)
save3(fig, "fig15_metabric")

fig, axs = plt.subplots(1, 3, figsize=(10, 3.4), gridspec_kw={"wspace": 0.3})
tp = lab[lab.group.isin(TC_COL)].dropna(subset=["PFI", "PFI.time"]).assign(years=lambda x: x["PFI.time"] / 365.25)
km_panel(axs[0], tp, "years", "PFI", "group", list(TC_COL), TC_COL, fs=7, xmax=10)
axs[0].set_title(f"TCGA PFI by clinical group (log-rank p = {fmt_p(srep['tcga_logrank']['TCGA_PFI_group']['p'])})", loc="left", fontsize=8.5)
m10 = mb_trunc(10)
km_panel(axs[1], m10.dropna(subset=["group"]), "years", "BCD", "group", ["ER+/HER2 not gained", "HER2 gained", "ER−/HER2 not gained"], MB_COL, fs=7)
axs[1].set_title(f"METABRIC breast-cancer death, 10 y (p = {fmt_p(srep['metabric_logrank']['MB_BCD_10y_group']['p'])})", loc="left", fontsize=8.5)
km_panel(axs[2], m10.dropna(subset=["pred"]), "years", "BCD", "pred", ["predicted other", "predicted ER−/HER2−"], MB_COL, fs=7)
axs[2].set_title(f"METABRIC by frozen TCGA classifier (p = {fmt_p(srep['metabric_logrank']['MB_BCD_10y_pred_group']['p'])})", loc="left", fontsize=8.5)
save3(fig, "fig16_survival")

# ====================================================== slide 8
slide_style()
ba_er, ba_h = best_model("ER"), best_model("HER2")
bb_er, bb_h = best_model("ER", "b_no_receptor_amplicon1.5Mb_PAM50"), best_model("HER2", "b_no_receptor_amplicon1.5Mb_PAM50")
base_er = main[(main.target == "ER") & main.model.str.startswith("single")].iloc[0]
base_h = main[(main.target == "HER2") & main.model.str.startswith("single")].iloc[0]
ctrl = crep["controls"]
fig = plt.figure(figsize=SLIDE)
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.3], left=0.05, right=0.985, top=0.73, bottom=0.11, wspace=0.62)
a0, a1, a2 = (fig.add_subplot(gs[i]) for i in range(3))
roc_panel(a0, "ER", fs=10.5); a0.set_title("Clinical ER (IHC)", loc="left")
roc_panel(a1, "HER2", fs=10.5); a1.set_title("Clinical HER2 (IHC/FISH)", loc="left")
auc_dots(a2, fs=10.5); a2.set_title("4 models × gene sets", loc="left")


def vs_base(r):
    return f"ΔAUC {r.dAUC_vs_baseline:+.3f} (p {r.NB_p:.2f})"


abl_er = abt[(abt.target == "ER") & abt.ablation.str.startswith("b_")]
abl_h = abt[(abt.target == "HER2") & abt.ablation.str.startswith("b_")]
any_beats = ((main.NB_p < 0.05) & (main.dAUC_vs_baseline > 0) & main.target.isin(["ER", "HER2"])).any()
head = ("Models match but do not beat the receptor gene alone" if not any_beats else
        "Some models beat the receptor gene alone")
fig.text(0.05, 0.965, f"{head} (ESR1 {base_er.AUC_pooled:.2f}, ERBB2 {base_h.AUC_pooled:.2f}).",
         fontsize=17, fontweight="bold", va="top")
hb = main[(main.target == "HER2") & main.ablation.str.startswith("b_") & main.model.isin(MODELS)].AUC_pooled
fig.text(0.05, 0.912, f"Without receptor genes, 17q12 amplicon and PAM50: ER unchanged ({bb_er.AUC_pooled:.2f}), "
         f"HER2 down to {hb.min():.2f}–{hb.max():.2f}.",
         fontsize=15, fontweight="bold", va="top", color=INK2)
fig.text(0.05, 0.862, f"Best model vs gene alone, Nadeau-Bengio corrected t on 15 paired folds: ER {SHORT[ba_er.model]} {vs_base(ba_er)}; "
         f"HER2 {SHORT[ba_h.model]} {vs_base(ba_h)}.\nShuffled labels: MCC {ctrl['ER']['shuffled_labels_L1_logistic']['MCC']:+.2f} (ER), "
         f"{ctrl['HER2']['shuffled_labels_L1_logistic']['MCC']:+.2f} (HER2). Filtering, scaling and feature selection refitted inside every training fold.",
         fontsize=11, color=INK2, va="top")
fig.text(0.05, 0.02, f"ER {int(base_er.n_pos)}+/{int(base_er.n - base_er.n_pos)}−, HER2 {int(base_h.n_pos)}+/"
         f"{int(base_h.n - base_h.n_pos)}−. Nested CV: outer 5-fold × 3, inner 5-fold. Ablated: minus ESR1/PGR/ERBB2, "
         f"genes ±1.5 Mb (or 2 Mb) of ERBB2, and PAM50.", fontsize=10.5, color=INK2)
save3(fig, "slide8_classification", exact=True)

# ====================================================== slide 9
c = rrep["concordance_TNBC_vs_HR"]
t_er = [r for r in mrep["transfer"] if r["model"] == "ER" and r["scaling"] == "z_within_cohort"][0]
t_h = [r for r in mrep["transfer"] if r["model"] == "HER2" and r["scaling"] == "z_within_cohort"][0]
t_tn = [r for r in mrep["transfer"] if r["model"] == "ERneg_HER2neg" and r["scaling"] == "z_within_cohort"][0]
fig = plt.figure(figsize=SLIDE)
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.15], left=0.06, right=0.985, top=0.70, bottom=0.12, wspace=0.42)
a0, a1, a2 = (fig.add_subplot(gs[i]) for i in range(3))
concordance_panel(a0, fs=11); a0.set_title("Protein (RPPA), same tumours", loc="left")
error_panel(a1, fs=11); a1.set_title("Where the ER model errs", loc="left")
transfer_panel(a2, fs=11); a2.set_title("METABRIC (1,980 tumours)", loc="left")
ea = rrep["er_error_analysis"]
fig.text(0.06, 0.965, f"Holds in protein and in METABRIC: ER AUC {t_er['AUC']:.2f}, HER2 {t_h['AUC']:.2f}, "
         f"ER−/HER2− {t_tn['AUC']:.2f}.", fontsize=17, fontweight="bold", va="top")
fig.text(0.06, 0.912, f"mRNA effects match protein (ρ = {c['spearman']:.2f}); {mrep['de_replication']['fraction']:.0%} of TCGA DEGs replicate in METABRIC.",
         fontsize=15, fontweight="bold", va="top", color=INK2)
fig.text(0.06, 0.862, f"IHC ER− tumours called ER+ have ER+-level ERα protein ({ea['tests']['ERneg_called_pos_above_that']} above the ER− 95th percentile): likely IHC false negatives. "
         f"DEGs: same sign, q < 0.05.\nMETABRIC labels: ER by IHC, HER2 by SNP6 copy number; models frozen on TCGA before any METABRIC label was read.\n"
         "HER2 transfers better than it cross-validates: copy-number labels track mRNA more closely than IHC/FISH calls.",
         fontsize=11, color=INK2, va="top")
save3(fig, "slide9_validation", exact=True)

# ====================================================== backup slide: survival
fig = plt.figure(figsize=SLIDE)
gs = fig.add_gridspec(1, 3, left=0.05, right=0.985, top=0.70, bottom=0.11, wspace=0.28)
a0, a1, a2 = (fig.add_subplot(gs[i]) for i in range(3))
km_panel(a0, tp, "years", "PFI", "group", list(TC_COL), TC_COL, fs=11, xmax=10)
a0.set_title(f"TCGA PFI (log-rank p = {fmt_p(srep['tcga_logrank']['TCGA_PFI_group']['p'])})", loc="left", fontsize=14)
km_panel(a1, m10.dropna(subset=["group"]), "years", "BCD", "group", ["ER+/HER2 not gained", "HER2 gained", "ER−/HER2 not gained"], MB_COL, fs=11)
a1.set_title(f"METABRIC, BC death (p = {fmt_p(srep['metabric_logrank']['MB_BCD_10y_group']['p'])})", loc="left", fontsize=14)
km_panel(a2, m10.dropna(subset=["pred"]), "years", "BCD", "pred", ["predicted other", "predicted ER−/HER2−"], MB_COL, fs=11)
a2.set_title(f"Frozen TCGA model (p = {fmt_p(srep['metabric_logrank']['MB_BCD_10y_pred_group']['p'])})", loc="left", fontsize=14)
for a_ in (a0, a1, a2):
    a_.set_ylim(0.5, 1.01)
tc = [x for x in srep["tcga_cox"] if x["label"] == "TCGA Cox PFI"][0]
mc = srep["metabric_cox"][0]
sp = [x for x in srep["metabric_cox"] if "split" in x["label"]]
fig.text(0.05, 0.965, f"Triple-negative disease relapses early: TNBC PFI HR {tc['HR']['TNBC']['exp(coef)']:.1f} in TCGA, replicated in METABRIC.",
         fontsize=17, fontweight="bold", va="top")
fig.text(0.05, 0.912, (f"METABRIC ER−/HER2− hazard crosses: HR {sp[0]['HR_0_5y']['ERneg_HER2notgained']:.1f} in years 0–5, "
                       f"{sp[0]['HR_5_10y']['ERneg_HER2notgained']:.1f} in years 5–10.") if sp else "",
         fontsize=15, fontweight="bold", va="top", color=INK2)
fig.text(0.05, 0.862, f"TCGA Cox PFI ~ group + age + stage: n {tc['n']}, {tc['events']} events, 5 parameters; TNBC HR {tc['HR']['TNBC']['exp(coef)']:.2f} "
         f"({tc['HR']['TNBC']['exp(coef) lower 95%']:.2f}–{tc['HR']['TNBC']['exp(coef) upper 95%']:.2f}), HER2+ {tc['HR']['HER2pos']['exp(coef)']:.2f}.\n"
         f"METABRIC Cox ~ group + age + grade + size: n {mc['n']}, {mc['events']} breast-cancer deaths in 10 years; reference ER+/HER2 not gained. PH fails for ER−/HER2− in METABRIC, so the "
         "effect is split at 5 years.\nExploratory; TCGA median follow-up is about 2.4 years.", fontsize=11, color=INK2, va="top")
save3(fig, "slideS1_survival", exact=True)
