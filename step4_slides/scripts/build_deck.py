"""Step 4a: build the 10-slide talk (16:9 .pptx with speaker notes) from step 2 and step 3 outputs.

Every number on the text slides (1, 2, 10) and in the speaker notes is read from the result files, and slides 3-9 are
the slide-ready figures that steps 2 and 3 regenerate themselves. So if the team changes the labels:
    rerun step 2 / step 3 (step3_analysis/config.json, then bash run_all.sh), then
    python scripts/build_deck.py            # writes BT3041_talk.pptx (+ PDF and previews if LibreOffice is installed)
Needs: python-pptx (pip install python-pptx). Optional: LibreOffice (soffice) and pdftoppm for the PDF and previews.
"""
import json, shutil, subprocess, sys
from datetime import datetime, timezone
from pathlib import Path
import pandas as pd
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE

HERE = Path(__file__).resolve().parents[1]
PLAN = HERE.parent
S2, S3 = PLAN / "step2_preprocessing", PLAN / "step3_analysis"
OUT = HERE / "BT3041_talk.pptx"

# ----------------------------------------------------------------------------------------------- numbers
J = lambda p: json.load(open(p))
pre = J(S2 / "results/preprocess/preprocess_report.json")
exp = J(S2 / "results/explore/explore_report.json")
clu = J(S2 / "results/cluster/cluster_report.json")
de = J(S3 / "results/de/main_summary.json")
her2de = J(S3 / "results/de/her2_summary.json")
enr = J(S3 / "results/enrichment/enrichment_report.json")
cls = J(S3 / "results/classify/classify_report.json")
rp = J(S3 / "results/rppa/rppa_report.json")
mb = J(S3 / "results/transfer/metabric_report.json")
sv = J(S3 / "results/survival/survival_report.json")
summ = pd.read_csv(S3 / "results/classify/summary.tsv", sep="\t")
main = summ[summ.tag == "main"]
MODELS = ["L1-logistic", "linear SVM", "Random Forest", "PLS-DA"]


def auc(target, model, abl="a_all_genes"):
    return main[(main.target == target) & (main.ablation == abl) & (main.model == model)].iloc[0]


def best(target, abl="a_all_genes"):
    d = main[(main.target == target) & (main.ablation == abl) & main.model.isin(MODELS)]
    return d.sort_values("AUC_pooled", ascending=False).iloc[0]


def single(target):
    return main[(main.target == target) & main.model.str.startswith("single")].iloc[0]


def transfer(model):
    return [r for r in mb["transfer"] if r["model"] == model and r["scaling"] == "z_within_cohort"][0]


grp = pre["clinical_group_final_supervised"]
island = exp["umap_island"]["15"]                    # UMAP with 15 neighbours (slide 4)
tn_island = island["TNBC_in_island"]                     # "129/157"
a, b = (int(x) for x in tn_island.split("/"))
abl_h = main[(main.target == "HER2") & main.ablation.str.startswith("b_") & main.model.isin(MODELS)].AUC_pooled
abl_e = main[(main.target == "ER") & main.ablation.str.startswith("b_") & main.model.isin(MODELS)].AUC_pooled
N = dict(
    hr=grp["HR+/HER2-"], her2=grp["HER2+"], tnbc=grp["TNBC"], inc=grp["Incomplete"],
    n_tum=pre["n_tumours_final"], n_sup=pre["n_supervised_set_final"],
    island_pct=f"{a / b:.0%}", island_n=island["n_island"],
    k2_size=min(clu["cluster_sizes"]["kmeans_k2"]),
    deg=de["DE_lfc1_padj05"], deg_up=de["DE_up"], deg_dn=de["DE_down"], wald=de["wald_padj05_no_threshold"],
    genes_tested=de["genes_tested"], her2_deg=her2de["DE_lfc1_padj05"],
    rob_age=enr["de_robustness_vs_main"]["age"]["fraction"], rob_idc=enr["de_robustness_vs_main"]["idc_only"]["fraction"],
    rob_strict=enr["de_robustness_vs_main"]["strict_conf"]["fraction"],
    gsea_n=enr["gsea_main_sets_FDR_lt_0.05"]["positive_NES"] + enr["gsea_main_sets_FDR_lt_0.05"]["negative_NES"],
    camera_n=enr["camera_check"]["sets_q_lt_0.05"],
    esr1=single("ER").AUC_pooled, erbb2=single("HER2").AUC_pooled,
    best_er=best("ER"), best_h=best("HER2"),
    abl_er_lo=abl_e.min(), abl_er_hi=abl_e.max(), abl_h_lo=abl_h.min(), abl_h_hi=abl_h.max(),
    mcc_er=cls["controls"]["ER"]["shuffled_labels_L1_logistic"]["MCC"],
    mcc_h=cls["controls"]["HER2"]["shuffled_labels_L1_logistic"]["MCC"],
    rho=rp["concordance_TNBC_vs_HR"]["spearman"], rho_p=rp["concordance_TNBC_vs_HR"]["perm_p"],
    rppa_genes=rp["concordance_TNBC_vs_HR"]["genes"], rppa_q=rp["mw_TNBC_vs_HR"]["q_lt_0.05"],
    rppa_ab=rp["mw_TNBC_vs_HR"]["antibodies"],
    mb_er=transfer("ER")["AUC"], mb_h=transfer("HER2")["AUC"], mb_tn=transfer("ERneg_HER2neg")["AUC"],
    mb_h_abl=transfer("HER2_ablated")["AUC"],
    mb_esr1=mb["single_gene_baselines"]["ER"]["AUC"], mb_erbb2=mb["single_gene_baselines"]["HER2"]["AUC"],
    rep=mb["de_replication"]["fraction"], rep_nes=mb["de_replication"]["hallmark_NES_spearman"],
    pfi_p=sv["tcga_logrank"]["TCGA_PFI_group"]["p"],
    er_fn=rp["er_error_analysis"]["IHC ER- called ER+"]["n"],
    er_fn_prot=rp["er_error_analysis"]["tests"]["ERneg_called_pos_above_that"],     # e.g. "8/10"
    rppa_n=rp["rppa_tumours_supervised"],
)
SHORT = {"L1-logistic": "L1-logistic", "linear SVM": "linear SVM", "Random Forest": "Random Forest", "PLS-DA": "PLS-DA"}
f2, f3, pc = (lambda x: f"{x:.2f}"), (lambda x: f"{x:.3f}"), (lambda x: f"{x:.0%}")
c = lambda x: f"{x:,}"

# ----------------------------------------------------------------------------------------------- style
INK, INK2, MUTED, GRID, BG = (RGBColor.from_string(h) for h in ("0B0B0B", "52514E", "898781", "E1E0D9", "F4F3EF"))
COL = {"HR+/HER2−": RGBColor.from_string("2A78D6"), "HER2+": RGBColor.from_string("EB6834"),
       "TNBC": RGBColor.from_string("1BAF7A")}
ACCENT = RGBColor.from_string("3D3C38")
FONT = "Arial"
# Presenter for each plan workstream (M1–M5); edit here to reassign slides
TEAM = {"M1": "Vijayavallabh", "M2": "Vedant", "M3": "Harshit", "M4": "Navya", "M5": "Anitha"}
W, H = Inches(40 / 3), Inches(7.5)

prs = Presentation()
prs.core_properties.title = "BT3041 term project talk"
prs.core_properties.author = prs.core_properties.last_modified_by = "BT3041 team"
prs.core_properties.comments = ""
prs.core_properties.created = prs.core_properties.modified = datetime.now(timezone.utc).replace(tzinfo=None)
prs.slide_width, prs.slide_height = W, H
BLANK = prs.slide_layouts[6]


def text(slide, x, y, w, h, runs, size=16, color=INK, bold=False, align=PP_ALIGN.LEFT, anchor=MSO_ANCHOR.TOP,
         spacing=1.0, space_after=0):
    """runs: str, or list of paragraphs; a paragraph is a str or a list of (text, {bold, color, size, italic})."""
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.word_wrap = True
    tf.margin_left = tf.margin_right = Inches(0.04)
    tf.margin_top = tf.margin_bottom = Inches(0.02)
    tf.vertical_anchor = anchor
    paras = runs if isinstance(runs, list) else [runs]
    for i, para in enumerate(paras):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align
        p.line_spacing = spacing
        p.space_after = Pt(space_after)
        for seg in (para if isinstance(para, list) else [(para, {})]):
            t, o = seg if isinstance(seg, tuple) else (seg, {})
            r = p.add_run()
            r.text = t
            r.font.name = FONT
            r.font.size = Pt(o.get("size", size))
            r.font.bold = o.get("bold", bold)
            r.font.italic = o.get("italic", False)
            r.font.color.rgb = o.get("color", color)
    return tb


def box(slide, x, y, w, h, fill=BG, line=None, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.06):
    s = slide.shapes.add_shape(shape, Inches(x), Inches(y), Inches(w), Inches(h))
    s.fill.solid(); s.fill.fore_color.rgb = fill
    if line is None:
        s.line.fill.background()
    else:
        s.line.color.rgb = line; s.line.width = Pt(1)
    s.shadow.inherit = False
    if shape == MSO_SHAPE.ROUNDED_RECTANGLE:
        s.adjustments[0] = radius
    return s


def number(slide, i):
    text(slide, 12.78, 7.16, 0.5, 0.3, str(i), size=10, color=MUTED, align=PP_ALIGN.RIGHT)


def notes(slide, who, time, body):
    tf = slide.notes_slide.notes_text_frame
    tf.text = f"[{TEAM.get(who, who)} · {time}]\n" + body.strip()


def figure_slide(i, png, alt, who, time, body):
    s = prs.slides.add_slide(BLANK)
    pic = s.shapes.add_picture(str(png), 0, 0, width=W, height=H)
    pic._element.nvPicPr.cNvPr.set("descr", alt)
    number(s, i)
    notes(s, who, time, body)
    return s


def title(slide, t, sub=None):
    text(slide, 0.6, 0.38, 12.1, 0.7, t, size=28, bold=True)
    if sub:
        text(slide, 0.6, 1.08, 12.1, 0.5, sub, size=18, color=INK2, bold=True)


# ============================================================================================ slide 1
s = prs.slides.add_slide(BLANK)
box(s, 0, 0, 0.22, 7.5, fill=ACCENT, shape=MSO_SHAPE.RECTANGLE)
text(s, 0.75, 0.75, 11.8, 0.4, "BT3041 Analysis and Interpretation of Biological Data · Term project", size=16, color=INK2)
text(s, 0.75, 1.35, 11.8, 1.9, "What does tumour gene expression reveal about the clinical receptor groups of breast cancer?",
     size=35, bold=True, spacing=0.95)
text(s, 0.75, 2.85, 11.8, 0.9,
     f"TCGA breast cancer RNA-seq ({c(N['n_tum'])} tumours), checked against protein (RPPA) in the same tumours "
     f"and an independent cohort (METABRIC, 1,980 tumours)", size=20, color=INK2)
for j, (g, n) in enumerate([("HR+/HER2−", N["hr"]), ("HER2+", N["her2"]), ("TNBC", N["tnbc"])]):
    x = 0.75 + j * 2.35
    box(s, x, 4.2, 2.15, 1.0, fill=COL[g])
    text(s, x, 4.25, 2.15, 0.9, [[(g, {"bold": True, "size": 18})], [(f"n = {n}", {"size": 16})]],
         color=RGBColor(255, 255, 255), align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE)
text(s, 7.85, 4.2, 4.9, 1.0, "Groups from clinical IHC/FISH tests (ER, PR, HER2), not from expression",
     size=15, color=INK2, anchor=MSO_ANCHOR.MIDDLE)
# Roll numbers live in roll_numbers.json at the project root, which is kept out of the public repository.
ROLLS = PLAN / "roll_numbers.json"
roll = J(ROLLS) if ROLLS.exists() else {}
if roll:
    text(s, 0.75, 5.62, 11.8, 0.35, "Team 16", size=14, color=INK2)
    for j, n in enumerate(TEAM.values()):
        text(s, 0.75 + j * 2.35, 5.95, 2.15, 0.6, [[(n, {})], [(roll.get(n, ""), {"size": 13, "color": MUTED})]], size=16)
    foot_y = 6.78
else:
    text(s, 0.75, 6.3, 11.8, 0.4, "Team 16: " + " · ".join(TEAM.values()), size=16)
    foot_y = 6.72
text(s, 0.75, foot_y, 11.8, 0.4, "IIT Madras · November 2026 · Code: github.com/Vijayavallabh/BT3041_Project", size=14, color=MUTED)
notes(s, "M1", "0:20", f"""
Good morning. We ask what tumour gene expression reveals about the three clinical receptor groups of breast cancer.
We use {c(N['n_tum'])} TCGA tumours, and we check every finding in protein from the same tumours and in an independent
cohort, METABRIC.""")

# ============================================================================================ slide 2
s = prs.slides.add_slide(BLANK)
title(s, "Receptor tests guide treatment. What does expression add?")
text(s, 0.6, 1.3, 5.6, 0.4, "Background", size=18, bold=True, color=INK2)
text(s, 0.6, 1.75, 5.6, 2.9, [
    "• Pathologists test ER and PR (IHC) and HER2 (IHC/FISH); the result picks endocrine therapy, anti-HER2 drugs or chemotherapy.",
    "• TNBC (all three negative) has no targeted therapy and relapses early.",
    "• PAM50 subtypes are computed from expression itself, so \"does expression separate PAM50?\" is circular. We use the clinical labels and keep PAM50 only as a comparator.",
], size=15, space_after=8)
text(s, 6.75, 1.3, 6.0, 0.4, "Research questions", size=18, bold=True, color=INK2)
qs = [("Q1", "Does label-free structure (PCA, UMAP, clustering) recover the receptor groups?", "slides 4–5"),
      ("Q2", "Which genes and pathways differ between TNBC and HR+/HER2− (DE + FDR, GSEA)?", "slides 6–7"),
      ("Q3", "Can expression predict clinical ER and HER2 status better than the receptor gene alone?", "slide 8")]
for j, (q, t, ref) in enumerate(qs):
    y = 1.78 + j * 0.95
    box(s, 6.75, y, 0.62, 0.62, fill=ACCENT)
    text(s, 6.75, y, 0.62, 0.62, q, size=15, bold=True, color=RGBColor(255, 255, 255), align=PP_ALIGN.CENTER,
         anchor=MSO_ANCHOR.MIDDLE)
    text(s, 7.5, y - 0.05, 5.25, 0.85, [[(t + " ", {}), (ref, {"color": MUTED, "size": 13})]], size=15)
text(s, 0.6, 4.75, 12.1, 0.4, "Validation design (slide 9): every result is checked outside the data it came from",
     size=18, bold=True, color=INK2)
val = [("Labels", "Clinical IHC/FISH for ER, PR, HER2. PAM50 only as a comparator; PAM50 genes left out of clustering and ablations."),
       ("Protein", f"TCGA RPPA on the same tumours ({rp['rppa_tumours_supervised']} with RNA and protein): do mRNA changes show up in protein?"),
       ("Independent cohort", "METABRIC microarrays: ER by IHC, HER2 by DNA copy number. Models frozen on TCGA before any METABRIC label was read.")]
for j, (h, t) in enumerate(val):
    x = 0.6 + j * 4.1
    box(s, x, 5.25, 3.9, 1.75)
    text(s, x + 0.15, 5.33, 3.6, 1.6, [[(h, {"bold": True, "size": 16})], t], size=14, space_after=4)
number(s, 2)
notes(s, "M1", "0:45", """
Clinically, breast cancer is grouped by three tests: ER and PR by immunohistochemistry, and HER2 by IHC or FISH.
These tests decide treatment, and triple-negative tumours have no targeted therapy.
Many studies use the PAM50 subtypes instead, but PAM50 is computed from expression, so that would be circular.
We use the clinical results as labels and PAM50 only as a comparator.
Our three questions: does unsupervised structure recover the groups; which genes and pathways differ; and can
expression predict clinical status better than the receptor gene alone.
We validate in protein from the same tumours and in METABRIC.""")

# ============================================================================================ slides 3-9
F2, F3 = S2 / "figures", S3 / "figures"
figure_slide(3, F2 / "slide3_sample_flow.png",
             "Sample flow from downloaded TCGA HiSeqV2, STAR counts, RPPA and METABRIC files to the analysis groups.",
             "M1", "0:45", f"""
Data preparation. We keep one primary tumour per patient and exclude 9 that failed an RNA identity check or had
redacted clinical records, leaving {c(N['n_tum'])} tumours for unsupervised work.
For HER2 we use FISH first, then the IHC score, and set equivocal or contradictory results to missing.
Women with all three receptors known give {N['hr']} HR+/HER2−, {N['her2']} HER2+ and {N['tnbc']} triple-negative tumours.
Differential expression uses raw STAR counts, keeping genes with at least 10 counts in 157 tumours.
RPPA adds protein for {N['rppa_n']} tumours, and METABRIC adds 1,980 tumours with independent labels.""")
figure_slide(4, F2 / "slide4_exploration.png",
             "UMAP of tumours coloured by clinical group with a separate basal-like island, and violins of ESR1, ERBB2, FOXC1 and MKI67.",
             "M2", "0:40", f"""
First, exploration without labels. On the UMAP of the 2,000 most variable genes, one island separates completely.
It holds {N['island_pct']} of all triple-negative tumours. HER2-positive tumours, in orange, mostly sit among the
hormone-receptor-positive tumours, at one edge.
The violins show why: ESR1 is the strongest single difference between TNBC and HR+, with a rank-biserial effect of −0.95.
ERBB2 is the clearest marker of HER2+, while FOXC1, a basal transcription factor, and MKI67, proliferation, are higher in TNBC.
All four differences are significant after Benjamini–Hochberg correction.""")
figure_slide(5, F2 / "slide5_clustering.png",
             "k-means clustering compared with PAM50, and adjusted Rand index of k=2 and k=5 clusterings against ER, clinical group and PAM50.",
             "M2", "0:45", f"""
Next, clustering, with the PAM50 genes left out. Consensus clustering is most stable at k = 2, with a PAC of 0.02:
a basal-like group of {N['k2_size']} tumours splits off. Its agreement with IHC ER status is an adjusted Rand index of 0.67.
At k = 5 a HER2-enriched cluster appears, but it is the least stable, so we call it a tendency, and luminal tumours
split differently from PAM50 LumA and LumB, a known continuum.
So for question 1: expression recovers the ER-negative, basal axis strongly, HER2 only weakly.""")
figure_slide(6, F3 / "slide6_de.png",
             "Volcano plot of TNBC versus HR+/HER2− differential expression and DEG counts in each sensitivity run.",
             "M3", "0:50", f"""
Differential expression, TNBC versus HR+/HER2−, with PyDESeq2 on raw counts, adjusting for histology.
Instead of filtering by fold change after the test, we test the 2-fold threshold directly. A plain Wald test calls
{N['wald'] / N['genes_tested']:.0%} of genes significant, which is not informative; with the threshold test,
{c(N['deg'])} genes differ more than 2-fold at 5% FDR.
ESR1 and PGR are positive controls. Beyond them, the ER regulon (GATA3, FOXA1, AGR3) is lower in TNBC, and basal
transcription factors (FOXC1, EN1, SOX10) are higher.
On the right, the list is robust: it keeps {pc(N['rob_idc'])} of genes if we drop lobular tumours, {pc(N['rob_age'])} adding age,
and {pc(N['rob_strict'])} with strict labels. HER2+ vs HR+ gives only {N['her2_deg']} genes, led by genes on the 17q12 amplicon around ERBB2.""")
figure_slide(7, F3 / "slide7_pathways.png",
             "Hallmark NES bar chart with CAMERA and label-permutation checks, TNBC vs PAM50 Basal 2x2 table, and METABRIC NES agreement.",
             "M3", "0:50", f"""
Pathways. GSEA on the ranked genes finds {N['gsea_n']} Hallmark sets at FDR 5%. TNBC ranks higher on proliferation
(E2F targets, G2M), MYC and immune programmes, and lower on estrogen response.
These gene-permutation q-values are known to be optimistic when genes are co-expressed, so we ran two tests that
account for correlation: CAMERA and label-permutation GSEA. They agree on the direction of every set, but only
early estrogen response stays significant. So we report effect sizes, not q ≈ 0.
Top right: clinical TNBC agrees with PAM50 Basal with kappa 0.78, close to the 2012 TCGA paper. Bottom right: the same
pathway shifts replicate in METABRIC, with a correlation of {f2(N['rep_nes'])}.""")
figure_slide(8, F3 / "slide8_classification.png",
             "ROC curves for clinical ER and HER2 and AUC of four models by gene set compared with the receptor gene alone.",
             "M4", "1:05", f"""
Question 3: can expression predict clinical status? We trained four models (L1 logistic regression, linear SVM,
Random Forest and PLS-DA) in nested cross-validation, with every filtering, scaling and selection step inside the
training folds, so there is no leakage. Shuffled labels give an MCC of about zero.
The key result: models match, but do not beat, the receptor gene alone. ESR1 alone has an AUC of {f3(N['esr1'])},
the best model {f3(N['best_er'].AUC_pooled)}; ERBB2 alone {f3(N['erbb2'])}, the best model {f3(N['best_h'].AUC_pooled)}.
The corrected t-tests, which account for overlapping folds, give p = {f2(N['best_er'].NB_p)} and {f2(N['best_h'].NB_p)}.
Then we removed the receptor genes, every gene within 1.5 Mb of ERBB2, and the PAM50 genes. ER is unchanged, at
{f2(N['abl_er_lo'])}–{f2(N['abl_er_hi'])}: ER status is written across a whole regulon, not one gene.
HER2 drops to {f2(N['abl_h_lo'])}–{f2(N['abl_h_hi'])}: some HER2 signal lies outside the amplicon, largely the lower-ER,
proliferative character of HER2+ tumours.""")
figure_slide(9, F3 / "slide9_validation.png",
             "mRNA versus protein fold changes, RPPA ER-alpha by ER model error type, and METABRIC transfer AUCs.",
             "M5", "1:00", f"""
Validation. Left: for {N['rppa_genes']} genes with a protein antibody, the mRNA fold change predicts the protein difference
with Spearman rho {f2(N['rho'])}; a permutation test gives p = {N['rho_p']:.3f}.
Middle: where does the ER model err? Of the {N['er_fn']} tumours that IHC calls ER-negative but the model calls positive,
{N['er_fn_prot'].replace('/', ' of the ')} with protein data have ER-alpha protein at ER-positive levels. These are probably IHC false negatives, not
model errors.
Right: METABRIC. We froze the models on TCGA before reading any METABRIC label. The ER model reaches AUC {f2(N['mb_er'])},
HER2 {f2(N['mb_h'])}, and ER−/HER2− versus rest {f2(N['mb_tn'])}; HER2 transfers even better than it cross-validated, because
METABRIC's HER2 label is copy number. And {pc(N['rep'])} of the TCGA DEGs replicate with the same sign.""")

# ============================================================================================ slide 10
s = prs.slides.add_slide(BLANK)
title(s, "Conclusions")
cards = [
    ("TNBC", "TNBC is a distinct transcriptional state",
     f"A basal-like island holds {N['island_pct']} of TNBC; {c(N['deg'])} genes differ > 2-fold: ER regulon down, basal factors, "
     f"proliferation and MYC up. HER2+ differs far less ({N['her2_deg']} genes, led by the 17q12 amplicon)."),
    ("HR+/HER2−", "The receptor gene already says almost everything",
     f"Models match but do not beat ESR1 ({f2(N['esr1'])}) and ERBB2 ({f2(N['erbb2'])}). ER survives removing ESR1 "
     f"({f2(N['abl_er_hi'])}): a regulon. HER2 falls to {f2(N['abl_h_lo'])}–{f2(N['abl_h_hi'])} without the amplicon."),
    ("HER2+", "Findings hold beyond the training data",
     f"mRNA matches protein (ρ = {f2(N['rho'])}); METABRIC transfer AUC {f2(N['mb_er'])} (ER), {f2(N['mb_h'])} (HER2); "
     f"{pc(N['rep'])} of DEGs replicate. Several model 'errors' look like IHC false negatives."),
]
for j, (g, h, t) in enumerate(cards):
    x = 0.6 + j * 4.1
    box(s, x, 1.3, 3.9, 3.05)
    box(s, x, 1.3, 3.9, 0.14, fill=COL["TNBC"] if j == 0 else ACCENT, shape=MSO_SHAPE.RECTANGLE)
    text(s, x + 0.18, 1.55, 3.55, 2.75, [[(f"{j + 1}. {h}", {"bold": True, "size": 17})], t], size=14.5, space_after=8)
text(s, 0.6, 4.55, 6.2, 0.4, "Limitations", size=17, bold=True, color=INK2)
text(s, 0.6, 4.95, 6.2, 2.2, [
    "• Pathway q-values: only early estrogen response survives correlation-aware tests.",
    "• DE not adjusted for tumour purity; proliferation and immune shifts are partly composition.",
    f"• HER2 labels from status text alone are noisy (they cost about 0.05 AUC).",
    "• TCGA follow-up is short: survival is exploratory (100 progression events).",
], size=14, space_after=4)
text(s, 7.1, 4.55, 5.6, 0.4, "Done now / for the final report", size=17, bold=True, color=INK2)
text(s, 7.1, 4.95, 5.6, 2.2, [
    [("Done: ", {"bold": True}), ("data preparation, exploration, clustering, DE, pathways, classification, RPPA and METABRIC validation.", {})],
    [("Report: ", {"bold": True}), ("survival in TCGA and METABRIC, label-rule sensitivity, comparison with 21 published numbers.", {})],
    [("Future: ", {"bold": True}), ("purity-adjusted DE; proteomic (CPTAC) labels.", {})],
    [("Code: ", {"bold": True}), ("github.com/Vijayavallabh/BT3041_Project", {})],
], size=14, space_after=6)
number(s, 10)
notes(s, "M5", "0:45", f"""
Three take-home messages. First, triple-negative breast cancer is a distinct transcriptional state: ER regulon off,
basal and proliferation programmes on. HER2-positive tumours differ far less, led by the amplicon.
Second, the receptor gene already carries almost all the information about clinical status; ER signal is spread over
a regulon, HER2 depends more on the amplicon.
Third, these results hold in protein and in an independent cohort.
Limitations: pathway significance shrinks once gene correlation is counted, we did not adjust for purity, and TCGA
survival is exploratory. Thank you; we are happy to take questions.""")

prs.save(OUT)
print("wrote", OUT)

# speaker script for rehearsal: the notes of every slide, with the timing plan
rows, lines = [], []
for i, sl in enumerate(prs.slides, 1):
    head, _, body = sl.notes_slide.notes_text_frame.text.partition("\n")
    who, time = head.strip("[]").split(" · ")
    rows.append((i, who, time))
    lines += [f"## Slide {i} ({who}, {time})", "", " ".join(body.split()), ""]
secs = sum(int(t.split(":")[0]) * 60 + int(t.split(":")[1]) for _, _, t in rows)
res = sum(int(t.split(":")[0]) * 60 + int(t.split(":")[1]) for i, _, t in rows if 4 <= i <= 9)
md = ["# Speaker script", "", f"Total {secs // 60}:{secs % 60:02d} of the 8:00 limit; methods and results (slides 4–9) take "
      f"{res // 60}:{res % 60:02d}, {res / secs:.0%} of the talk (the guidelines ask for at least 40%). About 130 words a minute.",
      "", "| Slide | Presenter | Time |", "|---|---|---|"] + [f"| {i} | {w} | {t} |" for i, w, t in rows] + ["", *lines]
(HERE / "SPEAKER_SCRIPT.md").write_text("\n".join(md))
print("wrote", HERE / "SPEAKER_SCRIPT.md")

# ----------------------------------------------------------------------------------------------- PDF + previews
soffice = shutil.which("soffice") or shutil.which("libreoffice")
if soffice and "--no-pdf" not in sys.argv:
    subprocess.run([soffice, "--headless", "--convert-to", "pdf", "--outdir", str(HERE), str(OUT)], check=True,
                   capture_output=True)
    pdf = OUT.with_suffix(".pdf")
    print("wrote", pdf)
    if shutil.which("pdftoppm"):
        prev = HERE / "previews"
        prev.mkdir(exist_ok=True)
        for old in prev.glob("slide-*.png"):
            old.unlink()
        subprocess.run(["pdftoppm", "-png", "-r", "60", str(pdf), str(prev / "slide")], check=True)
        print("previews in", prev)
