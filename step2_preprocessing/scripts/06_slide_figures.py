"""Step 6: two 16:9 slide composites with large type, drawn from saved tables (runs in seconds).

slide4_exploration: UMAP by clinical group + four marker violins (ESR1, ERBB2, FOXC1, MKI67)
slide5_clustering:  k-means k=5 vs PAM50 contingency + ARI of the label-free clusters against each reference
Both are saved at exactly 13.33 x 7.5 in (16:9), so they drop onto a widescreen slide at 100% scale.
Run from the project root:  python scripts/06_slide_figures.py   (after 02, 03 and 03b)
"""
import json, sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.common import FIG, PROC, RES, GROUP_ORDER, GROUP_COL, CLUSTER_COL, METHOD_COL, INK, INK2, MUTED, style, save, drop_excluded

style()
plt.rcParams.update({"font.size": 14, "axes.titlesize": 15, "axes.labelsize": 13, "xtick.labelsize": 12, "ytick.labelsize": 12})
SLIDE = (40 / 3, 7.5)
lab = drop_excluded(pd.read_csv(PROC / "tcga_labels.tsv", sep="\t", index_col=0))
X = pd.read_parquet(PROC / "tcga_hiseq_tumour_log2.parquet").loc[lab.index]
ex = json.load(open(RES / "explore/explore_report.json"))
mt = pd.read_csv(RES / "explore/marker_tests.tsv", sep="\t")
U = pd.read_csv(RES / "explore/umap_nn15.tsv", sep="\t", index_col=0).loc[lab.index]
NAME = {"HR+/HER2-": "HR+/HER2−", "HER2+": "HER2+", "TNBC": "TNBC"}

# ------------------------------------------------ slide 4
isl = ex["umap_island"]["15"]
tn_in, tn_all = (int(v) for v in isl["TNBC_in_island"].split("/"))
pct = tn_in / tn_all
fig = plt.figure(figsize=SLIDE)
gs = fig.add_gridspec(2, 3, width_ratios=[1.45, 1, 1], left=0.035, right=0.985, top=0.80, bottom=0.08, wspace=0.28, hspace=0.5)
ax = fig.add_subplot(gs[:, 0])
cg = lab.clinical_group
bg = ~cg.isin(GROUP_ORDER)
ax.scatter(U.UMAP1[bg], U.UMAP2[bg], s=10, c="#d9d8d2", lw=0, rasterized=True, label=f"labels incomplete ({bg.sum()})")
for g in GROUP_ORDER:
    m = cg == g
    ax.scatter(U.UMAP1[m], U.UMAP2[m], s=14, c=GROUP_COL[g], lw=0.2, edgecolors="white", rasterized=True, label=f"{NAME[g]} ({m.sum()})")
ax.legend(loc="upper left", fontsize=11, markerscale=1.8, handletextpad=0.2, borderaxespad=0.2)
# annotate the island in the empty lower-right corner, with a leader line up to the island
I = U[U.island]
ax.annotate(f"basal-like island, {isl['n_island']} tumours:\n{tn_in} TNBC ({pct:.0%} of all TNBC)\n"
            f"{isl['HR+/HER2- in island']} HR+/HER2−, {isl['HER2+ in island']} HER2+\n{isl['Incomplete in island']} incomplete",
            xy=(I.UMAP1.median(), I.UMAP2.min()), xytext=(0.98, 0.42), textcoords="axes fraction",
            ha="right", va="top", fontsize=11, color=INK, arrowprops=dict(arrowstyle="-", color=MUTED, lw=1, shrinkB=4))
ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
ax.set_title(f"UMAP of {len(lab):,} tumours (top 2,000 genes)", loc="left")
sv = lab[lab.supervised_set & lab.clinical_group.isin(GROUP_ORDER)]
nG = sv.clinical_group.value_counts()
for i, gname in enumerate(["ESR1", "ERBB2", "FOXC1", "MKI67"]):
    a = fig.add_subplot(gs[i // 2, 1 + i % 2])
    d = pd.DataFrame({"y": X.loc[sv.index, gname], "g": sv.clinical_group})
    sns.violinplot(data=d, x="g", y="y", order=GROUP_ORDER, hue="g", palette=GROUP_COL, legend=False, inner="quartile", saturation=1,
                   cut=0, linewidth=0.8, density_norm="width", ax=a)
    a.set_title(gname, fontstyle="italic", loc="left"); a.set_xlabel(""); a.set_ylabel("log2 expr." if i % 2 == 0 else "")
    a.set_xticks(range(3)); a.set_xticklabels(["HR+/\nHER2−", "HER2+", "TNBC"], fontsize=11)
    con = "HER2+ vs HR+/HER2-" if gname == "ERBB2" else "TNBC vs HR+/HER2-"
    r = mt[(mt.gene == gname) & (mt.contrast == con)].iloc[0]
    lo, hi = a.get_ylim(); a.set_ylim(lo, hi + (hi - lo) * 0.25)
    a.text(0.98, 0.97, f"{NAME[con.split(' vs')[0]]} vs HR+: r = {r.rank_biserial_r:+.2f}", transform=a.transAxes,
           ha="right", va="top", fontsize=11, color=INK2)
sel = mt.set_index(["gene", "contrast"]).loc[[("ESR1", "TNBC vs HR+/HER2-"), ("ERBB2", "HER2+ vs HR+/HER2-"),
                                               ("FOXC1", "TNBC vs HR+/HER2-"), ("MKI67", "TNBC vs HR+/HER2-")], "q_BH"]
assert (sel < 1e-10).all()
fig.text(0.035, 0.965, f"Expression sets TNBC apart: a basal-like island holds {pct:.0%} of TNBC.", fontsize=17,
         fontweight="bold", va="top")
fig.text(0.035, 0.915, "HER2+ tumours mostly overlap HR+/HER2− tumours; their clearest difference is ERBB2.", fontsize=15,
         fontweight="bold", va="top", color=INK2)
fig.text(0.035, 0.865, f"Violins: female supervised set (HR+/HER2− {nG['HR+/HER2-']}, HER2+ {nG['HER2+']}, TNBC {nG['TNBC']}); "
         "r = rank-biserial effect size; all four q < 1e-10 (Mann-Whitney, BH over 16 tests)", fontsize=11, color=INK2, va="top")
save(fig, "slide4_exploration", exact=True)

# ------------------------------------------------ slide 5
rep = json.load(open(RES / "cluster/cluster_report.json"))
k_star = rep["k_star_lowest_PAC"]
ag = pd.read_csv(RES / "cluster/cluster_agreement.tsv", sep="\t")
scan = pd.read_csv(RES / "cluster/k_selection.tsv", sep="\t", index_col=0)
ceil = json.load(open(RES / f"cluster/ari_ceiling_k{k_star}.json"))
pam2 = json.load(open(RES / "cluster/ari_pam50_two_groups_vs_refs.json"))
nmap = pd.read_csv(RES / "cluster/netanely_comparison.tsv", sep="\t").set_index("cluster").closest_Netanely2016_cluster
clus = pd.read_csv(PROC / "tcga_clusters.tsv", sep="\t", index_col=0).loc[lab.index]
fig = plt.figure(figsize=SLIDE)
gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1], left=0.17, right=0.98, top=0.71, bottom=0.25, wspace=0.42)
ax = fig.add_subplot(gs[0])
p4 = ["LumA", "LumB", "Her2", "Basal"]
t = pd.crosstab(clus.kmeans_k5.map(lambda c: f"C{c}"), lab.PAM50_pancan).reindex(columns=p4).fillna(0).astype(int)
sns.heatmap(t.div(t.sum(0), axis=1), annot=t.values, fmt="d", cmap="Blues", vmin=0, vmax=1, cbar=False, linewidths=2,
            linecolor="white", ax=ax, annot_kws={"size": 14})
ax.set_yticklabels([f"{c}  {nmap[c]}" for c in t.index], rotation=0, fontsize=12)
ax.set_xlabel("PAM50 subtype (comparator; Normal-like left out)"); ax.set_ylabel("")
a5 = ag[(ag.clustering == "kmeans_k5") & (ag.reference == "PAM50 (no Normal-like)")].iloc[0]
ax.set_title(f"k = 5 vs PAM50: ARI {a5.ARI:.2f} [{a5.ARI_lo:.2f}–{a5.ARI_hi:.2f}]", loc="left")
ax2 = fig.add_subplot(gs[1])
refs = ["ER status (IHC)", "Clinical group", "PAM50 (no Normal-like)"]
rlab = {"ER status (IHC)": "ER (IHC)", "Clinical group": "Clinical\ngroup", "PAM50 (no Normal-like)": "PAM50"}
off = {f"kmeans_k{k_star}": -0.13, "kmeans_k5": 0.13}
for j, (mth, nm) in enumerate([(f"kmeans_k{k_star}", f"k-means, k = {k_star}"), ("kmeans_k5", "k-means, k = 5")]):
    s = ag[ag.clustering == mth].set_index("reference").reindex(refs)
    y = np.arange(len(refs)) + off[mth]
    ax2.errorbar(s.ARI, y, xerr=[s.ARI - s.ARI_lo, s.ARI_hi - s.ARI], fmt="o", ms=9, lw=2, color=METHOD_COL[mth], label=nm)
for i, r in enumerate(refs):
    if r in pam2:
        ax2.plot(pam2[r], i, marker="D", ms=10, mfc=MUTED, mec=INK, mew=1.0, ls="none",
                 label="PAM50 as 2 groups\n(luminal vs Basal+Her2)" if i == 0 else None)
    b = ceil[r][f"best_ARI_with_{k_star}_merged_groups"]
    if b is not None:
        ax2.plot([b, b], [i + off[f"kmeans_k{k_star}"] - 0.11, i + off[f"kmeans_k{k_star}"] + 0.11], color=INK2, lw=2.5,
                 solid_capstyle="butt", label=f"best possible for k = {k_star}" if i == 0 else None)
ax2.set_yticks(range(len(refs))); ax2.set_yticklabels([rlab[r] for r in refs]); ax2.invert_yaxis()
ax2.set_xlim(0, 1); ax2.set_xlabel("Adjusted Rand index (95% bootstrap CI)")
ax2.legend(fontsize=11, loc="upper center", bbox_to_anchor=(0.45, -0.17), ncol=2)
ax2.set_title("Agreement with each reference", loc="left")
c1 = clus[f"kmeans_k{k_star}"] == 1
n1 = int(c1.sum())
tn1, tnN = int((c1 & (lab.clinical_group == "TNBC")).sum()), int((lab.clinical_group == "TNBC").sum())
er1, erN = int((c1 & (lab.ER == "Negative")).sum()), int((lab.ER == "Negative").sum())
b1, bN = int((c1 & (lab.PAM50_pancan == "Basal")).sum()), int((lab.PAM50_pancan == "Basal").sum())
assert k_star == 2, "the slide-5 headline is written for k* = 2; rewrite it if the most stable k changes"
fig.text(0.02, 0.965, f"Label-free clustering is most stable at k = 2: a basal-like group of {n1} tumours splits off",
         fontsize=16, fontweight="bold", va="top")
fig.text(0.02, 0.905, f"It holds {b1}/{bN} PAM50 Basal, {tn1}/{tnN} TNBC and {er1}/{erN} ER-negative tumours",
         fontsize=14, fontweight="bold", va="top", color=INK2)
fig.text(0.02, 0.85, "Finer clusters (k = 5) also recover a HER2-enriched group but split luminal tumours differently from PAM50",
         fontsize=12.5, va="top", color=INK2)
fig.text(0.02, 0.805, f"{len(clus):,} tumours, 2,000 most variable genes, PAM50 genes left out. k = 5 is less stable than k = 2 "
         f"(PAC {scan.PAC[5]:.2f} vs {scan.PAC[k_star]:.2f}; ARI {rep['sensitivity_ARI_vs_main']['top5000_excl_PAM50_k5']:.2f} "
         "against a 5,000-gene rerun).", fontsize=11, color=INK2, va="top")
fig.text(0.02, 0.015, "Heatmap: cell = count, shade = share of each PAM50 column. Row names: closest Netanely 2016 K = 5\n"
         "cluster by majority PAM50 class (a heuristic, not their assignment; they also clustered adjacent normals).", fontsize=11, color=INK2,
         linespacing=1.4)
save(fig, "slide5_clustering", exact=True)

for name in ["slide4_exploration", "slide5_clustering"]:
    w, h = plt.imread(FIG / f"{name}.png").shape[1::-1]
    assert abs(w / h - 16 / 9) < 0.01, f"{name} is not 16:9 ({w}x{h})"
print("done")
