"""Step 5: sample-flow diagrams and the library-size / identity QC figure. Every number is read from the processed files.

- figures/slide3_sample_flow: 16:9 slide, all four data sets (HiSeqV2, STAR, RPPA, METABRIC), 11 pt text;
  its footnotes go to figures/slide3_sample_flow_notes.txt (speaker notes).
- figures/fig1_sample_flow: report figure (7.5 in wide, 9 pt text at that size), the two TCGA RNA data sets as boxes,
  RPPA and METABRIC counts and the rules as text below.
- figures/figS0_qc_library_identity.

Run from the project root:  python scripts/05_sample_flow.py   (after 01_preprocess.py)
"""
import json, sys, textwrap
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.common import PROC, RES, FIG, GROUP_COL, INK, INK2, MUTED, style, save, drop_excluded

style()
rep = json.load(open(RES / "preprocess/preprocess_report.json"))
lab = pd.read_csv(PROC / "tcga_labels.tsv", sep="\t", index_col=0)
ml = pd.read_csv(PROC / "metabric_labels.tsv", sep="\t", index_col=0)
ok = drop_excluded(lab)
sv = ok[ok.supervised_set]
ex = lab[lab.exclude]
n_id = int(ex.identity_unconfirmed.sum())
n_red = int(ex.cdr_redacted.sum())
n_both = int((ex.identity_unconfirmed & ex.cdr_redacted).sum())
star_ok = ok[ok.star_available]
sv_rppa = sv[sv.rppa_available]
hs = rep["hiseq_sample_types"]
st = rep["star_sample_types"]
n_norm_red = len(rep["normals_excluded_cdr_redacted"])
n_norm = rep["n_adjacent_normals_hiseq"] - n_norm_red
g = sv.clinical_group.value_counts()
gs = sv[sv.star_available].clinical_group.value_counts()
gr = sv_rppa.clinical_group.value_counts()
mg = ml.contrast_group.value_counts()
n_male = int((ok.sex == "Male").sum())
n_sex_na = int(ok.sex.isna().sum())
er_her2 = rep["ERneg_HER2neg_vs_ERpos_HER2neg_supervised_star"]
assert len(ok) - len(sv) == n_male, "supervised set should drop only male patients"

# box texts, shared by both layouts (one string per row; SIDE = what leaves the HiSeqV2 flow at that row)
HISEQ = [
    f"{sum(hs.values()):,} RNA-seq profiles:\n{hs['01']:,} primary, {hs['11']} normal,\n{hs.get('06', 0)} metastatic",
    f"{hs['01']:,} primary tumours\n(1 per patient)",
    f"{len(ok):,} tumours\n(all unsupervised work)",
    f"supervised set (female*)\n{len(sv):,} tumours\nER, PR, HER2 known: {len(sv) - g.get('Incomplete', 0):,}",
]
SIDE = {1: f"{n_norm} normals kept apart\n(+{n_norm_red} redacted, dropped)\n{hs.get('06', 0)} metastatic dropped",
        2: f"{len(ex)} excluded ({n_both} in both)\nidentity failed: {n_id}\nCDR redacted: {n_red}",
        3: f"{n_male} male patients:\nunsupervised only\n{g.get('Incomplete', 0)} incomplete labels"}
STAR = [
    f"{sum(st.values()):,} profiles:\n{st['01']:,} primary vials,\n{st['11']} normal, {st.get('06', 0)} metastatic",
    f"vial dedup: keep 01A\n({len(rep['star_vial_duplicates'])} patients)\n→ {rep['star_primary_tumours_after_dedup']:,} tumours",
    f"same exclusions → {len(star_ok):,}\n{rep['star_genes_passing_filter']:,} genes ≥10 counts\nin ≥{rep['star_min_samples']} tumours",
    f"main DE contrast\n(supervised set): TNBC {gs['TNBC']}\nvs HR+/HER2− {gs['HR+/HER2-']}",
    f"matched to METABRIC\n(PR ignored): ER−/HER2−\n{er_her2['ER-/HER2-']} vs ER+/HER2− {er_her2['ER+/HER2-']}",
]
RPPA = [
    f"{sum(rep['rppa_sample_types'].values())} samples × {rep['rppa_shape_antibodies_x_samples'][0]} antibodies\n"
    f"{rep['rppa_sample_types']['01']} primary tumours",
    f"{rep['rppa_antibodies_total_with_gene']} total-protein antibodies\nmapped to genes\n"
    f"({len(rep['rppa_antibodies_isoform_ambiguous'])} pan-isoform, flagged)",
    f"{len(sv_rppa)} supervised tumours\nwith RNA + protein: TNBC {gr['TNBC']}\nHER2+ {gr['HER2+']} · HR+/HER2− {gr['HR+/HER2-']}",
]
MB = [
    f"{len(ml):,} tumours, Illumina HT-12\n{rep['genes_shared_tcga_hiseq_and_metabric']:,} genes shared with TCGA",
    "labels from other assays only:\nER by IHC, HER2 by SNP6\ncopy number",
    f"both labels known: {rep['metabric_both_labels_known']:,}",
    f"ER+/HER2 not gained {mg['ER+/HER2-notgained']:,}\nER−/HER2 not gained {mg['ER-/HER2-notgained']}\nHER2 gained {mg['HER2-gained']}",
]
FOOT = [
    "HER2 rule: FISH first; otherwise IHC score (3+ positive, 0/1+ negative, 2+ without FISH = missing, IHC 3+ with negative FISH = missing); "
    "IHC status text only when no score; status/score contradictions = missing. Groups need ER, PR and HER2 known.",
    f"Identity check: a tumour's HiSeqV2 and STAR profiles must correlate at ρ ≥ 0.80, above the 99th percentile of between-patient ρ "
    f"({rep['identity_qc']['between_patient_rho_p99']:.2f}), and be each other's best match. "
    f"*Recorded sex, or XIST/RPS4Y1 expression when blank (n={n_sex_na}). METABRIC's ER−/HER2-not-gained group is not called TNBC (no IHC PR).",
]
VALID = (f"Validation cohorts (boxes on slide 3). TCGA RPPA: {len(sv_rppa)} supervised tumours with RNA and protein "
         f"(TNBC {gr['TNBC']}, HER2+ {gr['HER2+']}, HR+/HER2− {gr['HR+/HER2-']}); {rep['rppa_antibodies_total_with_gene']} total-protein "
         f"antibodies mapped to genes. METABRIC: {len(ml):,} tumours (Illumina HT-12, {rep['genes_shared_tcga_hiseq_and_metabric']:,} genes "
         f"shared with TCGA), labels from other assays only (ER by IHC, HER2 by SNP6 copy number): ER+/HER2 not gained "
         f"{mg['ER+/HER2-notgained']:,}, ER−/HER2 not gained {mg['ER-/HER2-notgained']}, HER2 gained {mg['HER2-gained']}.")


def flow(figsize, fs, name, report):
    """report=False: 16:9 slide with all four data sets. report=True: narrow figure, TCGA RNA columns plus text below."""
    fig, ax = plt.subplots(figsize=figsize)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.94 if report else 0.92, bottom=0.01)
    ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

    def box(x, y, w, h, text, fc="#f4f3ef", ec="#c3c2b7", bold=False, color=INK):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4,rounding_size=1.2", fc=fc, ec=ec, lw=1))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, color=color,
                fontweight="bold" if bold else "normal", linespacing=1.3)

    if report:
        xh, wh, xs, ws, xt, wt = 0.5, 30, 32, 30, 63.5, 36
        ys, hh, ytop = [83, 69, 55, 41, 27], 9, 96.5
    else:
        xh, wh, xs, ws, xt, wt = 0.5, 17, 18.6, 15, 34.8, 17
        ys, hh, ytop = [82, 64, 46, 28, 8], 10, 97

    def down(x, i):
        ax.annotate("", xy=(x, ys[i + 1] + hh + 0.5), xytext=(x, ys[i] - 0.4), arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.2))

    def column(x, w, texts, bold_row=None, arrow_last=False):
        for i, t in enumerate(texts):
            box(x, ys[i], w, hh, t, bold=(i == bold_row))
            if i + 1 < len(texts) or arrow_last:
                down(x + w / 2, i)

    heads = [(xh, xs + ws - xh, "TCGA HiSeqV2\n(exploration, clustering, classifiers)"), (xt, wt, "TCGA STAR counts\n(differential expression)")]
    if not report:
        xr, wr, xm, wm = 53, 20.5, 74.8, 24.7
        heads += [(xr, wr, "TCGA RPPA\n(protein validation)"), (xm, wm, "METABRIC\n(external validation)")]
    for x, w, t in heads:
        ax.text(x + w / 2, ytop, t, ha="center", va="center", fontsize=fs + 0.5, fontweight="bold", color=INK2)
    column(xh, wh, HISEQ, bold_row=2, arrow_last=True)  # the last HiSeqV2 arrow points to the group chips
    for i, t in SIDE.items():
        box(xs, ys[i], ws, hh, t, fc="white", ec="#e34948", color=INK2)
    cw = (xs + ws - xh - 2 * 1.6) / 3  # group chips span the HiSeqV2 and side columns
    for i, k in enumerate(["HR+/HER2-", "HER2+", "TNBC"]):
        box(xh + i * (cw + 1.6), ys[4], cw, hh, f"{k.replace('-', '−')}\n{g[k]}", fc=GROUP_COL[k], ec=GROUP_COL[k], color="white", bold=True)
    column(xt, wt, STAR)
    if report:
        wrap = int(figsize[0] * 72 / (fs * 0.52))  # characters that fit across the canvas at this font size
        txt = "\n".join(textwrap.fill(t, wrap) for t in [VALID] + FOOT)
        ax.text(0.5, ys[4] - 4, txt, fontsize=fs - 0.5, color=INK2, va="top", linespacing=1.4)
    else:
        column(xr, wr, RPPA)
        column(xm, wm, MB)
    fig.suptitle("Sample flow: from downloaded files to analysis groups", fontsize=fs + (2 if report else 4), fontweight="bold",
                 x=0.01, ha="left")
    save(fig, name, exact=True)


flow((7.5, 9), 9, "fig1_sample_flow", report=True)
flow((40 / 3, 7.5), 11, "slide3_sample_flow", report=False)
(FIG / "slide3_sample_flow_notes.txt").write_text("Speaker notes / footnotes for slide 3 (not drawn on the slide):\n\n" + "\n\n".join(FOOT) + "\n")
w_, h_ = Image.open(FIG / "slide3_sample_flow.png").size
assert abs(w_ / h_ - 16 / 9) < 0.01, f"slide3 is not 16:9 ({w_}x{h_})"

# library sizes and same-sample agreement
iq = pd.read_csv(RES / "preprocess/identity_qc.tsv", sep="\t", index_col=0)
p99 = rep["identity_qc"]["between_patient_rho_p99"]
failed = rep["identity_qc"]["failed"]
fig, axes = plt.subplots(1, 2, figsize=(11, 3.8))
ls = lab.star_library_size.dropna() / 1e6
axes[0].hist(ls, bins=40, color="#2a78d6", edgecolor="white", lw=0.5)
axes[0].axvline(0.25 * ls.median(), color=MUTED, ls="--", lw=1)
axes[0].text(0.25 * ls.median() + 1, axes[0].get_ylim()[1] * 0.9, "flag threshold\n(25% of median)", fontsize=7, color=INK2)
axes[0].set_xlabel("STAR total counts per tumour (millions)"); axes[0].set_ylabel("tumours")
axes[0].set_title(f"Library size: median {ls.median():.0f} M, {len(rep['star_low_library_flagged'])} flagged")
ax = axes[1]
okm = ~iq.index.isin(list(failed))
ax.scatter(iq.best_between_patient_rho[okm], iq.self_rho[okm], s=6, c="#2a78d6", lw=0, alpha=0.6, label="passed")
ax.scatter(iq.best_between_patient_rho[~okm], iq.self_rho[~okm], s=26, c="#e34948", lw=0, label=f"failed ({(~okm).sum()})")
lo = min(iq.self_rho.min(), iq.best_between_patient_rho.min()) - 0.02
ax.plot([lo, 1], [lo, 1], color=MUTED, lw=1, ls=":")
ax.axhline(0.8, color=MUTED, ls="--", lw=1); ax.axhline(p99, color=INK2, ls="--", lw=1)
ax.text(lo + 0.005, 0.8 - 0.02, "ρ = 0.80", fontsize=7, color=INK2, va="top")
ax.text(lo + 0.005, p99 + 0.005, f"99th pct between-patient ρ = {p99:.2f}", fontsize=7, color=INK2, va="bottom")
# labels for the three crowded failures are placed by hand (data coordinates) with leader lines, clear of the reference lines
LABEL_AT = {"TCGA-A7-A0DC-01": (0.74, 0.857), "TCGA-AC-A2QH-01": (0.74, 0.831), "TCGA-A7-A13G-01": (0.905, 0.855)}
for s_ in failed:
    x_, y_ = iq.loc[s_, "best_between_patient_rho"], iq.loc[s_, "self_rho"]
    tx, ty = LABEL_AT.get(s_, (x_ + 0.006, y_))
    ax.annotate(s_.replace("TCGA-", "").replace("-01", ""), xy=(x_, y_), xytext=(tx, ty), fontsize=6.5, color=INK2,
                ha="right" if tx < x_ else "left", va="center", zorder=5,
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.9),
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6, shrinkA=0, shrinkB=3) if s_ in LABEL_AT else None)
ax.set_xlabel("best Spearman ρ with another patient's STAR profile"); ax.set_ylabel("Spearman ρ with own STAR profile")
ax.legend(fontsize=7, loc="lower right")
ax.set_title(f"Sample identity: median own ρ {iq.self_rho.median():.2f}; {len(failed)} excluded")
fig.tight_layout()
save(fig, "figS0_qc_library_identity")
print("done")
