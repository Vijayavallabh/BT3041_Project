"""Step 3.7: figures for differential expression and pathways (report figures + 16:9 slides 6 and 7).

fig11_volcano_de         volcano of TNBC vs HR+/HER2- (shrunk LFC vs threshold-test FDR), DEG counts across runs
fig12_hallmarks          Hallmark NES (prerank) with the two correlation-aware checks (phenotype permutation, CAMERA)
figS10_her2_volcano      HER2+ vs HR+/HER2-: the DE list is dominated by the 17q12 amplicon
figS11_nes_runs          Hallmark NES across all DE runs (robustness)
slide6_de, slide7_pathways   16:9 composites, text >= 11 pt
Every number drawn or written comes from results/ (no typed numbers). Run after scripts 10, 11 and 14.
"""
import sys, json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from adjustText import adjust_text

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.common3 import (RES3, PROC, GROUP_COL, INK, INK2, MUTED, GRID, UP_COL, DOWN_COL, SLIDE, style, slide_style,  # noqa: E402
                         save3, read_gmt)

DE, EN = RES3 / "de", RES3 / "enrichment"
er = json.load(open(EN / "enrichment_report.json"))
main = pd.read_csv(DE / "main_results.tsv", sep="\t", index_col=0)
G = pd.read_csv(EN / "gsea_prerank_all_runs.tsv", sep="\t")
cam = pd.read_csv(EN / "camera_hiseq.tsv", sep="\t").set_index("Term")
ph = pd.read_csv(EN / "gsea_phenotype_permutation_hiseq.tsv", sep="\t").set_index("Term")
mbr = json.load(open(RES3 / "transfer" / "metabric_report.json")) if (RES3 / "transfer" / "metabric_report.json").exists() else None
YCAP = 150
LABEL_UP = ["FOXC1", "EN1", "KRT5", "KRT17", "VGLL1", "ART3", "FABP7", "SOX10", "MMP7", "CLDN6", "COL9A1", "MPZ"]
LABEL_DN = ["FOXA1", "GATA3", "AGR3", "XBP1", "TFF1", "CA12", "SLC39A6", "MLPH"]
CONTROLS = ["ESR1", "PGR"]
RUN_NAME = {"main": "Main (all calls)", "age": "+ age in design", "idc_only": "IDC only", "strict_conf": "Strict labels",
            "matched_erher2": "ER−/HER2− vs ER+/HER2−"}


def volcano(ax, res, labels_up, labels_dn, controls, fs=9, cap=YCAP, lfc_lines=True, point=4):
    r = res.dropna(subset=["padj_lfc1"]).copy()
    r["y"] = -np.log10(r.padj_lfc1.clip(lower=1e-300))
    r["capped"] = r.y > cap
    r["y"] = r.y.clip(upper=cap)
    ns = ~r.DE
    ax.scatter(r.log2FC_shrunk[ns], r.y[ns], s=point * 0.6, c="#d9d8d2", lw=0, rasterized=True)
    for d, col in [("up", UP_COL), ("down", DOWN_COL)]:
        m = r.direction == d
        ax.scatter(r.log2FC_shrunk[m & ~r.capped], r.y[m & ~r.capped], s=point, c=col, lw=0, rasterized=True, alpha=0.8)
        ax.scatter(r.log2FC_shrunk[m & r.capped], r.y[m & r.capped], s=point * 3, c=col, marker="^", lw=0, rasterized=True)
    if lfc_lines:
        for x in (-1, 1):
            ax.axvline(x, color=MUTED, lw=0.8, ls="--", zorder=0)
    ax.axhline(-np.log10(0.05), color=MUTED, lw=0.8, ls=":", zorder=0)
    bysym = r.dropna(subset=["symbol"]).sort_values("y", ascending=False).drop_duplicates("symbol").set_index("symbol")
    texts = []
    for g in labels_up + labels_dn + controls:
        if g in bysym.index:
            x, y = bysym.loc[g, "log2FC_shrunk"], bysym.loc[g, "y"]
            ax.scatter([x], [y], s=point * 4, facecolors="none", edgecolors=INK, lw=0.7, zorder=4)
            texts.append(ax.text(x, y, g, fontsize=fs, fontstyle="italic", color=INK2 if g in controls else INK, zorder=5))
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.5), expand=(1.3, 1.6))
    ax.set_xlabel("shrunk log2 fold change")
    ax.set_ylabel(f"−log10 FDR (H0: |log2FC| ≤ 1; capped at {cap})")
    return r


def nes_panel(ax, run="main", fs=9, show_checks=True, max_sets=None):
    d = G[(G.run == run) & (G.FDR_q < 0.05)].sort_values("NES")
    if max_sets and len(d) > max_sets:
        d = pd.concat([d.head(max_sets // 2), d.tail(max_sets - max_sets // 2)]).drop_duplicates()
    y = np.arange(len(d))
    cols = [UP_COL if v > 0 else DOWN_COL for v in d.NES]
    ax.barh(y, d.NES, color=cols, height=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels([t.replace("_", " ").lower().capitalize().replace("Myc", "MYC").replace("E2f", "E2F").replace("G2m", "G2M")
                        .replace("Tnfa", "TNFα").replace("nfkb", "NF-κB").replace("Mtorc1", "mTORC1").replace("Il6", "IL6")
                        .replace("jak stat3", "JAK-STAT3").replace("Il2 stat5", "IL2-STAT5").replace("Kras", "KRAS")
                        .replace("Uv ", "UV ").replace("Wnt beta", "WNT β").replace("Dna", "DNA").replace("v1", "V1").replace("v2", "V2")
                        for t in d.Term], fontsize=fs)
    ax.axvline(0, color=INK2, lw=0.8)
    ax.set_xlabel("NES (GSEA prerank, Wald statistic)")
    if show_checks:
        xmax = d.NES.abs().max()
        for yi, t in zip(y, d.Term):
            c = cam.loc[t] if t in cam.index else None
            p = ph.loc[t] if t in ph.index else None
            mark_c = "●" if c is not None and c.q_BH < 0.05 else ("◐" if c is not None and c.p < 0.05 else "○")
            mark_p = "●" if p is not None and p.FDR_q < 0.05 else ("◐" if p is not None and p.nom_p < 0.05 else "○")
            ax.text(xmax * 1.12, yi, mark_c, ha="center", va="center", fontsize=fs, color=INK)
            ax.text(xmax * 1.32, yi, mark_p, ha="center", va="center", fontsize=fs, color=INK)
        ax.text(xmax * 1.12, len(d) - 0.3, "CAMERA", ha="center", va="bottom", fontsize=fs - 1, color=INK2, rotation=90)
        ax.text(xmax * 1.32, len(d) - 0.3, "label perm.", ha="center", va="bottom", fontsize=fs - 1, color=INK2, rotation=90)
        ax.set_xlim(-xmax * 1.08, xmax * 1.42)
    ax.set_ylim(-0.7, len(d) + (3.2 if show_checks else 0.3))
    return d


def runs_bar(ax, fs=9):
    rob = er["de_robustness_vs_main"]
    summ = er["de_summary"]
    order = ["main", "age", "idc_only", "strict_conf", "matched_erher2"]
    vals = [summ[r]["DE_lfc1_padj05"] for r in order]
    y = np.arange(len(order))[::-1]
    ax.barh(y, vals, color=["#3d3c38"] + ["#9a998f"] * (len(order) - 1), height=0.65)
    for yi, r, v in zip(y, order, vals):
        s = f"{v:,}" + ("" if r == "main" else f"  ({rob[r]['fraction']:.0%} of main kept)")
        ax.text(v + max(vals) * 0.02, yi, s, va="center", fontsize=fs, color=INK2)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{RUN_NAME[r]}\n{summ[r]['n_test']} vs {summ[r]['n_ref']}" for r in order], fontsize=fs)
    ax.set_xlim(0, max(vals) * 1.9)
    ax.set_xticks([])
    ax.spines["bottom"].set_visible(False)


# ====================================================== report figures
style()
fig, axs = plt.subplots(1, 2, figsize=(7.5, 3.9), gridspec_kw={"width_ratios": [1.35, 1], "wspace": 0.55})
volcano(axs[0], main, LABEL_UP[:8], LABEL_DN[:8], CONTROLS, fs=7, point=3)
s = er["de_summary"]["main"]
axs[0].set_title(f"TNBC ({s['n_test']}) vs HR+/HER2− ({s['n_ref']}), ~ histology + group", loc="left", fontsize=9)
x0 = (0 - axs[0].get_xlim()[0]) / (axs[0].get_xlim()[1] - axs[0].get_xlim()[0])
axs[0].text(x0 - 0.02, 0.9, f"{s['DE_down']:,} lower ←", color=DOWN_COL, transform=axs[0].transAxes, fontsize=7.5, va="center", ha="right")
axs[0].text(x0 + 0.02, 0.9, f"→ {s['DE_up']:,} higher", color=UP_COL, transform=axs[0].transAxes, fontsize=7.5, va="center", ha="left")
runs_bar(axs[1], fs=7)
axs[1].set_title("DEGs (|log2FC| > 1, FDR 5%) per run", loc="left", fontsize=9)
fig.text(0.01, -0.02, f"Plain Wald test without the fold-change threshold: {s['wald_padj05_no_threshold']:,} of {s['genes_tested']:,} genes "
         f"at FDR 5% ({s['wald_padj05_no_threshold'] / s['genes_tested']:.0%}).", fontsize=7, color=INK2)
save3(fig, "fig11_volcano_de")

fig, ax = plt.subplots(figsize=(6.2, 6.4))
nes_panel(ax, fs=7.5)
ax.set_title("Hallmarks at prerank FDR < 0.05 (TNBC vs HR+/HER2−)", loc="left", fontsize=9)
fig.text(0.01, 0.0, "Columns: correlation-aware checks on HiSeqV2 (● q < 0.05, ◐ nominal p < 0.05, ○ neither).\n"
         "CAMERA = competitive test with inter-gene correlation (Wu & Smyth 2012); label perm. = phenotype-permutation GSEA.",
         fontsize=7, color=INK2, va="top")
save3(fig, "fig12_hallmarks")

her2 = pd.read_csv(DE / "her2_results.tsv", sep="\t", index_col=0)
amp = pd.read_csv(PROC / "erbb2_amplicon_genes.tsv", sep="\t", index_col=0)
her2["amplicon_2Mb"] = her2.index.isin(amp.index[amp.within_2000kb])
hs = json.load(open(DE / "her2_summary.json"))
fig, ax = plt.subplots(figsize=(5.5, 4))
volcano(ax, her2, ["ERBB2", "GRB7", "STARD3", "PGAP3", "MIEN1", "PNMT", "TCAP", "PPP1R1B", "MED1", "CDK12"], ["CHGB", "PCSK1"], [], fs=7, point=3)
n_amp = int((her2.DE & her2.amplicon_2Mb).sum())
ax.set_title(f"HER2+ ({hs['n_test']}) vs HR+/HER2− ({hs['n_ref']}): {hs['DE_lfc1_padj05']} DEGs, {n_amp} within 2 Mb of ERBB2",
             loc="left", fontsize=8.5)
save3(fig, "figS10_her2_volcano")

nes = G.pivot(index="Term", columns="run", values="NES")
nes = nes.loc[nes["main"].sort_values().index]
fig, ax = plt.subplots(figsize=(6.5, 9))
im = ax.imshow(nes.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-3.2, vmax=3.2)
ax.set_yticks(range(len(nes))); ax.set_yticklabels(nes.index, fontsize=6)
ax.set_xticks(range(nes.shape[1])); ax.set_xticklabels(nes.columns, rotation=40, ha="right", fontsize=7)
plt.colorbar(im, ax=ax, shrink=0.4, label="NES")
ax.set_title("Hallmark NES across DE runs", loc="left", fontsize=9)
save3(fig, "figS11_nes_runs")

# ====================================================== slide 6
slide_style()
fig = plt.figure(figsize=SLIDE)
gs = fig.add_gridspec(1, 2, width_ratios=[1.35, 1], left=0.06, right=0.985, top=0.74, bottom=0.11, wspace=0.42)
ax = fig.add_subplot(gs[0])
volcano(ax, main, LABEL_UP, LABEL_DN, CONTROLS, fs=11, point=6)
ax.set_title("TNBC vs HR+/HER2− (PyDESeq2, ~ histology + group)", loc="left")
x0 = (0 - ax.get_xlim()[0]) / (ax.get_xlim()[1] - ax.get_xlim()[0])
ax.text(x0 - 0.02, 0.9, f"{s['DE_down']:,} lower ←", color=DOWN_COL, transform=ax.transAxes, fontsize=13, va="center", ha="right", fontweight="bold")
ax.text(x0 + 0.02, 0.9, f"→ {s['DE_up']:,} higher in TNBC", color=UP_COL, transform=ax.transAxes, fontsize=13, va="center", ha="left", fontweight="bold")
a2 = fig.add_subplot(gs[1])
runs_bar(a2, fs=11)
a2.set_title("DEGs in each sensitivity run", loc="left")
rob = er["de_robustness_vs_main"]
fig.text(0.06, 0.965, f"{s['DE_lfc1_padj05']:,} genes differ more than 2-fold between TNBC and HR+/HER2− tumours (FDR 5%).",
         fontsize=17, fontweight="bold", va="top")
fig.text(0.06, 0.912, f"Robust to dropping lobular tumours ({rob['idc_only']['fraction']:.0%} kept), adding age "
         f"({rob['age']['fraction']:.0%}) and strict labels ({rob['strict_conf']['fraction']:.0%}).",
         fontsize=15, fontweight="bold", va="top", color=INK2)
fig.text(0.06, 0.862, f"The 2-fold threshold is tested (lfc_null = 1), not applied afterwards: a plain Wald test calls {s['wald_padj05_no_threshold']:,} of "
         f"{s['genes_tested']:,} genes ({s['wald_padj05_no_threshold'] / s['genes_tested']:.0%}).\nESR1 and PGR are positive controls; "
         "beyond them: basal transcription factors (FOXC1, EN1, SOX10) up, the ER regulon (GATA3, FOXA1, XBP1, AGR3) down.",
         fontsize=11, color=INK2, va="top")
fig.text(0.06, 0.02, f"STAR counts, female supervised set; shrunk LFC (apeGLM-type prior); "
         f"HER2+ vs HR+/HER2−: {hs['DE_lfc1_padj05']} DEGs, {n_amp} of them within 2 Mb of ERBB2 (17q12 amplicon).",
         fontsize=11, color=INK2)
save3(fig, "slide6_de", exact=True)

# ====================================================== slide 7
fig = plt.figure(figsize=SLIDE)
gs = fig.add_gridspec(2, 2, width_ratios=[1.3, 0.75], height_ratios=[1, 1], left=0.20, right=0.86, top=0.74, bottom=0.08,
                      wspace=0.42, hspace=0.75)
ax = fig.add_subplot(gs[:, 0])
d = nes_panel(ax, fs=11, max_sets=16)
ax.set_title(f"Hallmarks at prerank FDR < 0.05 ({len(d)} of {int((G[G.run == 'main'].FDR_q < 0.05).sum())} shown)", loc="left", fontsize=14)
# cross-assay agreement 2x2
cx = er["cross_assay_agreement"][0]
a2 = fig.add_subplot(gs[0, 1])
tab = np.array([[cx["both"], cx["only_clinical"]], [cx["only_PAM50"], cx["neither"]]])
a2.imshow(np.log1p(tab), cmap="Greens", vmin=0, vmax=np.log1p(tab.max()) * 1.3)
for i in range(2):
    for j in range(2):
        a2.text(j, i, f"{tab[i, j]}", ha="center", va="center", fontsize=15, color=INK, fontweight="bold" if i == j == 0 else None)
a2.set_xticks([0, 1]); a2.set_xticklabels(["Basal", "other"])
a2.set_xlabel("PanCan PAM50 call", fontsize=11)
a2.set_yticks([0, 1]); a2.set_yticklabels(["TNBC", "HR+/HER2−\nor HER2+"])
a2.set_title(f"Clinical TNBC vs PAM50 Basal", loc="left", fontsize=14)
a2.text(1.08, 0.5, f"κ = {cx['kappa']:.2f}\n\n{cx['pct_clinical_called_PAM50']:.0f}% of TNBC\nare Basal\n(TCGA 2012: 86%)\n\n{cx['pct_PAM50_with_clinical']:.0f}% of Basal\nare TNBC\n(TCGA 2012: 80%)",
        transform=a2.transAxes, fontsize=11, color=INK2, va="center")
a3 = fig.add_subplot(gs[1, 1])
if mbr is not None:
    j = pd.read_csv(RES3 / "transfer" / "hallmark_nes_tcga_vs_metabric.tsv", sep="\t", index_col=0)
    a3.scatter(j.NES, j.NES_MB, s=28, c=[UP_COL if v > 0 else DOWN_COL for v in j.NES], lw=0)
    lim = max(j.NES.abs().max(), j.NES_MB.abs().max()) * 1.1
    a3.plot([-lim, lim], [-lim, lim], color=GRID, lw=1, zorder=0)
    a3.axhline(0, color=GRID, lw=0.8); a3.axvline(0, color=GRID, lw=0.8)
    for t in ["E2F_TARGETS", "ESTROGEN_RESPONSE_EARLY", "INTERFERON_GAMMA_RESPONSE"]:
        a3.annotate(t.split("_")[0].replace("INTERFERON", "IFNγ").replace("ESTROGEN", "Estrogen early"), (j.loc[t, "NES"], j.loc[t, "NES_MB"]),
                    xytext=(4, -12), textcoords="offset points", fontsize=10, color=INK2)
    rp = mbr["de_replication"]
    a3.set_title(f"METABRIC replicates: NES ρ = {rp['hallmark_NES_spearman']:.2f}", loc="left", fontsize=14)
    a3.set_xlabel("NES, TCGA ER−/HER2− vs ER+/HER2−", fontsize=11); a3.set_ylabel("NES, METABRIC", fontsize=11)
fig.text(0.02, 0.965, "TNBC ranks higher on proliferation, MYC and immune programmes, lower on estrogen response.",
         fontsize=17, fontweight="bold", va="top")
camE = cam.loc["ESTROGEN_RESPONSE_EARLY"]
fig.text(0.02, 0.912, "Correlation-aware checks agree on direction; only early estrogen response stays significant.", fontsize=15, fontweight="bold", va="top", color=INK2)
fig.text(0.02, 0.862, "NES: GSEA prerank on the DESeq2 Wald statistic (gene-permutation q-values are optimistic when genes are co-expressed).\n"
         "Marks: ● q < 0.05, ◐ nominal p < 0.05, ○ neither, for CAMERA (competitive test with inter-gene correlation) and label-permutation GSEA.",
         fontsize=11, color=INK2, va="top")
save3(fig, "slide7_pathways", exact=True)
