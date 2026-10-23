"""Step 5: the five report figures, sized for an A4 page with 2 cm margins (6.69 in wide), text >= 6.5 pt.

Fig 1  exploration and clustering   UMAP by clinical group, k = 2 cluster composition, ARI with each reference
Fig 2  differential expression      volcano (TNBC vs HR+/HER2-) and Hallmark NES with the correlation-aware checks
Fig 3  classification                ROC for ER and HER2 (gene alone, best model, ablated) and AUC of every model
Fig 4  validation                    mRNA vs protein, ER errors vs ERalpha protein, METABRIC transfer
Fig 5  survival                      TCGA PFI and METABRIC breast-cancer death (KM), Cox hazard ratios

Panels come from step 3's own plotting functions (scripts/17_figures_de.py, scripts/18_figures_validation.py), loaded
up to their "report figures" marker, so they read the same result files and stay in step 2/3 colours. Nothing is
typed by hand: if the labels change, rerun step 2/3, then this script, then build_report.py.
Run with the step 3 environment (step3_analysis/requirements.txt).
"""
import json, sys
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
PLAN = HERE.parent
S2, S3 = PLAN / "step2_preprocessing", PLAN / "step3_analysis"
FIG = HERE / "figures"
FIG.mkdir(exist_ok=True)
W = 6.69  # text width in inches (A4, 2 cm margins)
MARKER = "# ====================================================== report figures"


def load_step3(script):
    """Execute a step 3 figure script up to its report-figure marker and return its namespace (data + panel functions)."""
    path = S3 / "scripts" / script
    src = path.read_text()
    ns = {"__file__": str(path), "__name__": "step3_panels"}
    exec(compile(src.split(MARKER)[0], str(path), "exec"), ns)
    return ns


de = load_step3("17_figures_de.py")
va = load_step3("18_figures_validation.py")
plt = de["plt"]
INK, INK2, MUTED, GRID, GROUP_COL = (de[k] for k in ("INK", "INK2", "MUTED", "GRID", "GROUP_COL"))
de["style"]()
plt.rcParams.update({"font.size": 7.5, "axes.titlesize": 8.5, "axes.labelsize": 7.5, "xtick.labelsize": 7,
                     "ytick.labelsize": 7, "legend.fontsize": 6.8, "lines.linewidth": 1.5})


def tag(ax, letter, x=-0.16, y=1.06):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="left", color=INK)


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight", pad_inches=0.02)
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print("saved", name)


# ----------------------------------------------------------------------------------------- Fig 1
lab_all = pd.read_csv(S2 / "data/processed/tcga_labels.tsv", sep="\t", index_col=0)
lab_all = lab_all[~lab_all.exclude]
um = pd.read_csv(S2 / "results/explore/umap_nn15.tsv", sep="\t", index_col=0).join(lab_all[["clinical_group"]])
cl = pd.read_csv(S2 / "data/processed/tcga_clusters.tsv", sep="\t", index_col=0)
agr = pd.read_csv(S2 / "results/cluster/cluster_agreement.tsv", sep="\t")
crep = json.load(open(S2 / "results/cluster/cluster_report.json"))
ceil = json.load(open(S2 / "results/cluster/ari_ceiling_k2.json"))
GNAME = {"HR+/HER2-": "HR+/HER2−", "HER2+": "HER2+", "TNBC": "TNBC", "Incomplete": "incomplete"}

fig, axs = plt.subplots(1, 3, figsize=(W, 2.35), gridspec_kw={"width_ratios": [1.15, 0.95, 1.25], "wspace": 0.55})
ax = axs[0]
for g in ["Incomplete", "HR+/HER2-", "HER2+", "TNBC"]:
    d = um[um.clinical_group == g]
    ax.scatter(d.UMAP1, d.UMAP2, s=2.2, lw=0, c=GROUP_COL.get(g, "#cfcec7"), label=f"{GNAME[g]} ({len(d)})", rasterized=True)
ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
ax.legend(loc="lower right", bbox_to_anchor=(1.04, -0.02), markerscale=3, handletextpad=0.1, borderaxespad=0, fontsize=6.2)
ax.set_title("UMAP, 2,000 genes", loc="left")
tag(ax, "a", x=-0.1)

ax = axs[1]
k2 = cl["kmeans_k2"].reindex(lab_all.index)
basal = k2.value_counts().idxmin()                     # the small (basal-like) cluster
tab = pd.crosstab(lab_all.clinical_group, k2.eq(basal)).reindex(["HR+/HER2-", "HER2+", "TNBC", "Incomplete"])
y = np.arange(len(tab))[::-1]
for yi, (g, r) in zip(y, tab.iterrows()):
    n = r.sum()
    ax.barh(yi, r[True] / n, color=GROUP_COL.get(g, "#b5b4ab"), height=0.62)
    ax.barh(yi, r[False] / n, left=r[True] / n, color="#ebeae4", height=0.62)
    inside = r[True] / n > 0.6
    ax.text(r[True] / n + (-0.02 if inside else 0.02), yi, f"{r[True]}/{n}", va="center", ha="right" if inside else "left",
            fontsize=6.5, color="white" if inside else INK2)
ax.set_yticks(y); ax.set_yticklabels([GNAME[g] for g in tab.index])
ax.set_xlim(0, 1); ax.set_xlabel("share in basal-like cluster")
ax.set_title(f"k = 2 (PAC {crep['k_scan']['2']['PAC']:.2f})", loc="left")
tag(ax, "b", x=-0.5)

ax = axs[2]
refs = [("ER status (IHC)", "ER (IHC)"), ("Clinical group", "clinical group"), ("PAM50 (no Normal-like)", "PAM50")]
runs = [("kmeans_k2", "k-means, k = 2", "#5b4fd6"), ("kmeans_k5", "k-means, k = 5", "#2c8a3c")]
for i, (ref, nm) in enumerate(refs):
    for j, (run, rn, col) in enumerate(runs):
        r = agr[(agr.clustering == run) & (agr.reference == ref)].iloc[0]
        ax.errorbar(r.ARI, i + 0.18 * (2 * j - 1), xerr=[[r.ARI - r.ARI_lo], [r.ARI_hi - r.ARI]], fmt="o", ms=3.5, color=col,
                    lw=1.2, label=rn if i == 0 else None)
for i, (ref, _) in enumerate(refs):
    if ref in ceil:
        v = ceil[ref]["best_ARI_with_2_merged_groups"]
        ax.plot([v, v], [i - 0.35, i + 0.35], color=INK2, lw=1.6, label="best possible, k = 2" if i == 0 else None)
ax.set_yticks(range(3)); ax.set_yticklabels([nm for _, nm in refs]); ax.set_ylim(3.55, -0.5)
ax.set_xlim(0, 1.02); ax.set_xlabel("adjusted Rand index (95% CI)")
ax.legend(loc="lower left", fontsize=6.2, handletextpad=0.3, borderaxespad=0.1)
ax.set_title("Agreement of clusters", loc="left")
tag(ax, "c", x=-0.42)
save(fig, "report_fig1_exploration")

# ----------------------------------------------------------------------------------------- Fig 2
fig, axs = plt.subplots(1, 2, figsize=(W, 3.0), gridspec_kw={"width_ratios": [1.05, 1], "wspace": 0.75})
de["volcano"](axs[0], de["main"], ["FOXC1", "EN1", "KRT5", "SOX10", "FABP7"], ["FOXA1", "GATA3", "AGR3", "XBP1", "TFF1"],
              de["CONTROLS"], fs=6.5, point=2.2)
axs[0].set_ylabel("−log10 FDR (H0: |log2FC| ≤ 1)")
ms = de["er"]["de_summary"]["main"]
axs[0].set_title("TNBC vs HR+/HER2−", loc="left")
axs[0].text(0.99, 0.9, f"{ms['DE_up']:,} higher\nin TNBC", transform=axs[0].transAxes, ha="right", va="top", fontsize=6.5, color=de["UP_COL"])
axs[0].text(0.02, 0.2, f"{ms['DE_down']:,} lower\nin TNBC", transform=axs[0].transAxes, ha="left", va="top", fontsize=6.5, color=de["DOWN_COL"])
tag(axs[0], "a", x=-0.2)
g_all = de["G"]
gm = g_all[(g_all.run == "main") & (g_all.FDR_q < 0.05)]
keep = pd.concat([gm[gm.NES < 0], gm[gm.NES > 0].nlargest(13, "NES")])
de["G"] = pd.concat([g_all[g_all.run != "main"], keep])      # top 13 positive + every negative set
de["nes_panel"](axs[1], fs=6.3)
de["G"] = g_all
axs[1].set_xlabel("NES (GSEA prerank)")
axs[1].set_title("Hallmarks, prerank FDR < 0.05", loc="left")
tag(axs[1], "b", x=-0.75)
save(fig, "report_fig2_de_pathways")

# ----------------------------------------------------------------------------------------- Fig 3
fig, axs = plt.subplots(1, 3, figsize=(W, 2.45), gridspec_kw={"width_ratios": [1, 1, 0.95], "wspace": 0.45})
va["roc_panel"](axs[0], "ER", fs=5.8); axs[0].set_title("IHC ER", loc="left")
va["roc_panel"](axs[1], "HER2", fs=5.8); axs[1].set_title("IHC/FISH HER2", loc="left")
for a in axs[:2]:
    leg = a.get_legend(); leg.set_title("AUC"); leg.get_title().set_fontsize(6)
va["auc_dots"](axs[2], fs=6); axs[2].set_title("All models, nested CV", loc="left")
axs[2].set_yticklabels([t.get_text().split(" · ")[0] for t in axs[2].get_yticklabels()], fontsize=6.3)
axs[2].set_xlabel("ROC-AUC (95% CI)")
for a, l in zip(axs, "abc"):
    tag(a, l, x=-0.28 if l != "c" else -0.42)
save(fig, "report_fig3_classification")

# ----------------------------------------------------------------------------------------- Fig 4
fig, axs = plt.subplots(1, 3, figsize=(W, 2.5), gridspec_kw={"width_ratios": [1, 1, 1.05], "wspace": 0.55})
va["concordance_panel"](axs[0], fs=6.3); axs[0].set_title("mRNA vs protein", loc="left")
axs[0].set_xlabel("mRNA log2FC (TNBC vs HR+)"); axs[0].set_ylabel("RPPA median difference")
c_tn = va["rrep"]["concordance_TNBC_vs_HR"]
axs[0].texts[-1].set_text(f"ρ = {c_tn['spearman']:.2f}\nperm. p = {c_tn['perm_p']:.3f}\n{c_tn['genes']} genes")
axs[0].texts[-1].set_position((0.97, 0.04)); axs[0].texts[-1].set_ha("right"); axs[0].texts[-1].set_va("bottom")
va["error_panel"](axs[1], fs=6.8); axs[1].set_title("ER errors", loc="left")
old = axs[1].get_legend()
hl, ll = old.legend_handles, [t.get_text() for t in old.get_texts()]
old.remove()
axs[1].legend(hl, ll, loc="upper center", bbox_to_anchor=(0.5, -0.36), ncol=5, fontsize=5.8, handletextpad=0.05,
              columnspacing=0.5, frameon=False, title=None)
axs[1].set_ylim(top=axs[1].get_ylim()[1] + 0.8)
axs[1].set_xlabel("IHC ER vs out-of-fold call", fontsize=6.5)
va["transfer_panel"](axs[2], fs=6.6); axs[2].set_title("Transfer to METABRIC", loc="left")
axs[2].legend_.remove()
axs[2].set_yticklabels(["ER", "ER, ablated", "HER2", "HER2, ablated", "ER−/HER2−"])
from matplotlib.lines import Line2D  # noqa: E402
h = [Line2D([], [], marker="o", ls="", mfc="white", mec=va["COL_A"], label="TCGA CV"),
     Line2D([], [], marker="o", ls="", color=va["COL_B"], label="METABRIC"),
     Line2D([], [], marker="|", ls="", color=va["COL_BASE"], ms=8, mew=2, label="gene alone")]
axs[2].legend(handles=h, fontsize=6, loc="lower left", frameon=False, borderaxespad=0.1, handletextpad=0.2)
for a, l in zip(axs, "abc"):
    tag(a, l, x=-0.3 if l != "c" else -0.62)
save(fig, "report_fig4_validation")

# ----------------------------------------------------------------------------------------- Fig 5
srep = va["srep"]
fig, axs = plt.subplots(1, 3, figsize=(W, 2.3), gridspec_kw={"width_ratios": [1, 1, 0.85], "wspace": 0.4})
lab = va["lab"]
tp = lab[lab.group.isin(va["TC_COL"])].dropna(subset=["PFI", "PFI.time"]).assign(years=lambda x: x["PFI.time"] / 365.25)
va["km_panel"](axs[0], tp, "years", "PFI", "group", list(va["TC_COL"]), va["TC_COL"], fs=6.8, xmax=10)
axs[0].set_title("TCGA, progression-free", loc="left", fontsize=7.5)
axs[0].get_legend().set_title(f"log-rank p = {va['fmt_p'](srep['tcga_logrank']['TCGA_PFI_group']['p'])}", prop={"size": 6.3})
m10 = va["mb_trunc"](10)
va["km_panel"](axs[1], m10.dropna(subset=["group"]), "years", "BCD", "group",
               ["ER+/HER2 not gained", "HER2 gained", "ER−/HER2 not gained"], va["MB_COL"], fs=6.8)
axs[1].set_title("METABRIC, BC death", loc="left", fontsize=7.5)
axs[1].get_legend().set_title(f"log-rank p = {va['fmt_p'](srep['metabric_logrank']['MB_BCD_10y_group']['p'])}", prop={"size": 6.3})
for a in axs[:2]:
    a.set_ylim(0.3, 1.02); a.set_ylabel("event-free probability")
    a.get_legend()._legend_box.align = "left"
    for t in a.get_legend().get_texts():
        t.set_text(t.get_text().replace("HER2 not gained", "HER2 not gain.").replace(" ev.)", " ev)"))
# hazard ratios: TCGA Cox (PFI) and METABRIC split at 5 years
ax = axs[2]
tc = srep["tcga_cox"][0]["HR"]
mbc = [c for c in srep["metabric_cox"] if c["label"].startswith("METABRIC Cox breast-cancer death, 10 y")][0]["HR"]
split = [c for c in srep["metabric_cox"] if "split at 5 y" in c["label"]][0]
rows = [("TNBC, TCGA", tc["TNBC"]["exp(coef)"], tc["TNBC"]["exp(coef) lower 95%"], tc["TNBC"]["exp(coef) upper 95%"], GROUP_COL["TNBC"]),
        ("HER2+, TCGA", tc["HER2pos"]["exp(coef)"], tc["HER2pos"]["exp(coef) lower 95%"], tc["HER2pos"]["exp(coef) upper 95%"], GROUP_COL["HER2+"]),
        ("HER2 gain, MB", mbc["HER2_gained"]["exp(coef)"], mbc["HER2_gained"]["exp(coef) lower 95%"], mbc["HER2_gained"]["exp(coef) upper 95%"], GROUP_COL["HER2+"]),
        ("ER−/HER2− ≤5 y, MB", split["HR_0_5y"]["ERneg_HER2notgained"], None, None, GROUP_COL["TNBC"]),
        ("ER−/HER2− >5 y, MB", split["HR_5_10y"]["ERneg_HER2notgained"], None, None, GROUP_COL["TNBC"])]
for i, (nm, hr, lo, hi, col) in enumerate(rows):
    if lo is None:
        ax.plot(hr, i, "D", ms=4, color=col)
    else:
        ax.errorbar(hr, i, xerr=[[hr - lo], [hi - hr]], fmt="o", ms=4, color=col, lw=1.2)
ax.axvline(1, color=MUTED, lw=0.8, ls="--")
ax.set_xscale("log"); ax.set_xlim(0.35, 6); ax.set_xticks([0.5, 1, 2, 4]); ax.set_xticklabels(["0.5", "1", "2", "4"])
ax.xaxis.set_minor_formatter(plt.NullFormatter())
for i, r in enumerate(rows):
    ax.text(0.37, i - 0.3, r[0], fontsize=6.3, color=INK2, va="bottom", ha="left")
ax.set_yticks([]); ax.spines["left"].set_visible(False); ax.set_ylim(len(rows) - 0.5, -0.9)
ax.set_xlabel("adjusted hazard ratio")
ax.set_title("Cox models", loc="left", fontsize=7.5)
for a, l in zip(axs, "abc"):
    tag(a, l, x=-0.25 if l != "c" else -0.08)
save(fig, "report_fig5_survival")
