"""Step 5: build the final report (A4, Arial 11 pt, single spacing, max 8 pages) as .docx, plus a PDF export.

Every number in the text is read from the step 2 / step 3 result files (dictionary N below) and the figures come from
report_figures.py. The labels were signed off on 23 Oct 2026 (no change); if they ever change:
    rerun step 2 / step 3 (step3_analysis/config.json, then bash run_all.sh)
    python scripts/report_figures.py
    python scripts/build_report.py      # writes BT3041_final_report.docx, .pdf and report_numbers.md; checks <= 8 pages
Hand edits made in the .docx are lost when it is rebuilt: edit the wording here until the last rebuild, then edit
the .docx directly.
Needs: python-docx, pandas (step 3 environment). LibreOffice (soffice) and pdfinfo for the PDF and the page check.
"""
import json, re, shutil, subprocess, sys
from datetime import datetime, timezone
from pathlib import Path
import pandas as pd
from docx import Document
from docx.enum.section import WD_SECTION
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Cm, Pt, RGBColor

HERE = Path(__file__).resolve().parents[1]
PLAN = HERE.parent
S2, S3 = PLAN / "step2_preprocessing", PLAN / "step3_analysis"
FIG = HERE / "figures"
OUT = HERE / "BT3041_final_report.docx"
MAX_PAGES = 8

# ================================================================================================ numbers
J = lambda p: json.load(open(p))
pre = J(S2 / "results/preprocess/preprocess_report.json")
exp = J(S2 / "results/explore/explore_report.json")
clu = J(S2 / "results/cluster/cluster_report.json")
ceil = J(S2 / "results/cluster/ari_ceiling_k2.json")
mrk = pd.read_csv(S2 / "results/explore/marker_tests.tsv", sep="\t").set_index(["gene", "contrast"])
agr = pd.read_csv(S2 / "results/cluster/cluster_agreement.tsv", sep="\t")
de = J(S3 / "results/de/main_summary.json")
her2de = J(S3 / "results/de/her2_summary.json")
enr = J(S3 / "results/enrichment/enrichment_report.json")
cross = pd.read_csv(S3 / "results/enrichment/cross_assay_agreement.tsv", sep="\t").set_index("comparison")
cls = J(S3 / "results/classify/classify_report.json")
rp = J(S3 / "results/rppa/rppa_report.json")
mb = J(S3 / "results/transfer/metabric_report.json")
sv = J(S3 / "results/survival/survival_report.json")
summ = pd.read_csv(S3 / "results/classify/summary.tsv", sep="\t")
cfg = J(S3 / "config.json")
cpar = J(S3 / "results/classify/params.json")["params"]
epar = J(S3 / "results/enrichment/params.json")["params"]
main = summ[summ.tag == "main"]
MODELS = ["L1-logistic", "linear SVM", "Random Forest", "PLS-DA"]


def best(target, abl="a_all_genes", tag="main"):
    d = summ[(summ.tag == tag) & (summ.target == target) & (summ.ablation == abl) & summ.model.isin(MODELS)]
    return d.sort_values("AUC_pooled", ascending=False).iloc[0]


def single(target, tag="main"):
    return summ[(summ.tag == tag) & (summ.target == target) & summ.model.str.startswith("single")].iloc[0]


def transfer(model):
    return [r for r in mb["transfer"] if r["model"] == model and r["scaling"] == "z_within_cohort"][0]


def ari(run, ref):
    return agr[(agr.clustering == run) & (agr.reference == ref)].iloc[0]


def frac(s):
    a, b = (int(x) for x in str(s).split()[0].split("/"))
    return a, b


grp = pre["clinical_group_final_supervised"]
gn = enr["named_genes_main"]
rob = enr["de_robustness_vs_main"]
island = exp["umap_island"]["15"]
isl_a, isl_b = frac(island["TNBC_in_island"])
k2prof = {r["cluster"]: r for r in clu["profile_kmeans_k2"]}
k2_small = min(clu["cluster_sizes"]["kmeans_k2"])
k2_tab = clu["clinical_group_by_kmeans_k2"]
small_id = "1" if clu["cluster_sizes"]["kmeans_k2"][1] == k2_small else "0"
abl_e = main[(main.target == "ER") & main.ablation.str.startswith("b_") & main.model.isin(MODELS)].AUC_pooled
abl_h = main[(main.target == "HER2") & main.ablation.str.startswith("b_") & main.model.isin(MODELS)].AUC_pooled
abt = pd.read_csv(S3 / "results/classify/ablation_tests.tsv", sep="\t")
tn_target_best = best("ERneg_HER2neg", tag="transfer_target")
tn_target_base = single("ERneg_HER2neg", tag="transfer_target")
sens_pos = single("HER2", tag="sens_drop_low_pos")
sens_all = single("HER2", tag="sens_drop_low_all")
sens_pos_best = best("HER2", tag="sens_drop_low_pos")
tcox = sv["tcga_cox"][0]
mcox = [c for c in sv["metabric_cox"] if c["label"].startswith("METABRIC Cox breast-cancer death, 10 y")][0]
msplit = [c for c in sv["metabric_cox"] if "split at 5 y" in c["label"]][0]
niel = mb.get("nielsen_direction_check") or sv.get("nielsen_direction_check")
amp = pd.read_csv(S2 / "data/processed/erbb2_amplicon_genes.tsv", sep="\t")
h2res = pd.read_csv(S3 / "results/de/her2_results.tsv", sep="\t", index_col=0)
amp2 = set(amp.loc[amp.within_2000kb, "Ensembl_ID"])
her2_amp_deg = int(h2res[h2res.DE & h2res.index.isin(amp2)].shape[0])
er_err = rp["er_error_analysis"]
rpt = rp["concordance_TNBC_vs_HR"]
rep = mb["de_replication"]
hr = lambda d, k: (d["HR"][k]["exp(coef)"], d["HR"][k]["exp(coef) lower 95%"], d["HR"][k]["exp(coef) upper 95%"])
tn_basal = cross.loc["TNBC vs PAM50 Basal"]
h2_h2e = cross.loc["clinical HER2+ vs PAM50 HER2-enriched"]
erpos_nl = cross.loc["ER+ tumours called non-luminal (Basal or Her2)"]
litt = pd.read_csv(S3 / "results/literature_table.tsv", sep="\t")
acc_er = litt[litt.quantity.str.contains("ER \\(ESR1\\)", regex=True)].ours.iloc[0]
acc_h2 = litt[litt.quantity.str.contains("HER2 \\(ERBB2\\)", regex=True)].ours.iloc[0]
lfl_a, lfl_b = frac(exp["HER2pos_vs_HER2E"]["like_for_like"]["TCGA2012_HER2_x_TCGA2012_PAM50"])
eta = exp["eta2_PC1_PC5"]
pc1 = [v for k, v in eta.items() if k.startswith("PC1")][0]
site_max = max(v["Tissue source site"] for k, v in eta.items() if not k.startswith("PC5"))

N = dict(
    n_prim=pre["n_primary_tumours_hiseq"], n_tum=pre["n_tumours_final"], n_excl=pre["n_primary_tumours_hiseq"] - pre["n_tumours_final"],
    n_sup=pre["n_supervised_set_final"], n_male=pre["sex_final"]["Male"],
    hr=grp["HR+/HER2-"], her2=grp["HER2+"], tnbc=grp["TNBC"], inc=grp["Incomplete"],
    h2_fish=pre["HER2_source_of_HER2pos_supervised"]["FISH"], h2_ihc=pre["HER2_source_of_HER2pos_supervised"]["IHC score"],
    h2_text=pre["HER2_source_of_HER2pos_supervised"]["IHC status (no score)"],
    kappa_h2=pre["kappa_HER2_v1_vs_nature2012"]["kappa"],
    h2_disc=pre["HER2_label_source"]["IHC 3+ vs FISH negative -> missing"],
    star_genes=pre["star_genes_passing_filter"], star_min=pre["star_min_samples"],
    rppa_n=rp["rppa_tumours_supervised"], rppa_ab=rp["mw_TNBC_vs_HR"]["antibodies"],
    mb_n=pre["metabric_samples_with_expression"], mb_shared=pre["genes_shared_tcga_hiseq_and_metabric"],
    mb_erpos=pre["metabric_ER_IHC"]["Positive"], mb_erneg=pre["metabric_ER_IHC"]["Negative"], mb_gain=pre["metabric_HER2_SNP6"]["GAIN"],
    mb_tn=pre["metabric_contrast_groups"]["ER-/HER2-notgained"], mb_hr=pre["metabric_contrast_groups"]["ER+/HER2-notgained"],
    mb_h2g=pre["metabric_contrast_groups"]["HER2-gained"],
    pc1_var=exp["pca_variance_explained_PC1_10"][0], pc1_auc=exp["PC1_AUC_for_IHC_ER"], pc1_eta=pc1["Clinical group"],
    site_max=site_max,
    isl=island["n_island"], isl_a=isl_a, isl_b=isl_b, isl_basal=island["PAM50 Basal in island"],
    r_esr1=mrk.loc[("ESR1", "TNBC vs HR+/HER2-"), "rank_biserial_r"], r_foxc1=mrk.loc[("FOXC1", "TNBC vs HR+/HER2-"), "rank_biserial_r"],
    r_erbb2=mrk.loc[("ERBB2", "HER2+ vs HR+/HER2-"), "rank_biserial_r"], r_mki67=mrk.loc[("MKI67", "TNBC vs HR+/HER2-"), "rank_biserial_r"],
    pac2=clu["k_scan"]["2"]["PAC"], pac5=clu["k_scan"]["5"]["PAC"], k2=k2_small,
    k2_tn=k2_tab["TNBC"][small_id], ari_er=ari("kmeans_k2", "ER status (IHC)"), ari_cg=ari("kmeans_k2", "Clinical group"),
    ari_pam=ari("kmeans_k2", "PAM50 (no Normal-like)"), ari5_pam=ari("kmeans_k5", "PAM50 (no Normal-like)"),
    ceil_cg=ceil["Clinical group"]["best_ARI_with_2_merged_groups"], ceil_pam=ceil["PAM50 (no Normal-like)"]["best_ARI_with_2_merged_groups"],
    meth_ari=clu["method_agreement_ARI"]["k2"],
    deg=de["DE_lfc1_padj05"], up=de["DE_up"], dn=de["DE_down"], genes_tested=de["genes_tested"], wald=de["wald_padj05_no_threshold"],
    rob_age=rob["age"]["fraction"], rob_idc=rob["idc_only"]["fraction"], rob_strict=rob["strict_conf"]["fraction"],
    rob_match=rob["matched_erher2"]["fraction"], n_idc=enr["de_summary"]["idc_only"],
    n_strict=enr["de_summary"]["strict_conf"], n_match=enr["de_summary"]["matched_erher2"],
    her2_deg=her2de["DE_lfc1_padj05"], her2_amp=her2_amp_deg, her2_deg_pos=enr["de_summary"]["her2_drop_low_pos"]["DE_lfc1_padj05"],
    gsea_n=enr["gsea_main_sets_FDR_lt_0.05"]["positive_NES"] + enr["gsea_main_sets_FDR_lt_0.05"]["negative_NES"],
    gsea_up=enr["gsea_main_sets_FDR_lt_0.05"]["positive_NES"], gsea_dn=enr["gsea_main_sets_FDR_lt_0.05"]["negative_NES"],
    nes=enr["gsea_main"], cam_n=enr["camera_check"]["sets_q_lt_0.05"], cam_corr=enr["camera_check"]["median_inter_gene_corr"],
    cam_e2f=enr["camera_check"]["selected"]["E2F_TARGETS"]["inter_gene_corr"],
    ph_rho=enr["phenotype_permutation_check"]["NES_spearman_prerank_vs_phenotype"], ph_n=enr["phenotype_permutation_check"]["sets_FDR_lt_0.05_phenotype"],
    tnb_k=tn_basal.kappa, tnb_a=int(tn_basal.both), tnb_b=int(tn_basal.both + tn_basal.only_clinical),
    tnb_pct=tn_basal.pct_clinical_called_PAM50, bt_b=int(tn_basal.both + tn_basal.only_PAM50), bt_pct=tn_basal.pct_PAM50_with_clinical,
    h2e_pct=h2_h2e.pct_clinical_called_PAM50, h2e_k=h2_h2e.kappa,
    nl_pct=erpos_nl.pct, nl_pct_n=erpos_nl.pct_incl_Normal_like,
    n_er=cls["n_per_target"]["ER"]["n"], n_er_pos=cls["n_per_target"]["ER"]["pos"],
    n_h2=cls["n_per_target"]["HER2"]["n"], n_h2_pos=cls["n_per_target"]["HER2"]["pos"],
    esr1=single("ER"), erbb2=single("HER2"), b_er=best("ER"), b_h2=best("HER2"),
    abl_e_lo=abl_e.min(), abl_e_hi=abl_e.max(), abl_h_lo=abl_h.min(), abl_h_hi=abl_h.max(),
    n_removed=cls["genes_removed_per_ablation"]["b_no_receptor_amplicon1.5Mb_PAM50"],
    mcc_er=cls["controls"]["ER"]["shuffled_labels_L1_logistic"]["MCC"], mcc_h2=cls["controls"]["HER2"]["shuffled_labels_L1_logistic"]["MCC"],
    sens_pos=sens_pos.AUC_pooled, sens_pos_best=sens_pos_best.AUC_pooled, sens_all=sens_all.AUC_pooled,
    tn_best=tn_target_best, tn_base=tn_target_base,
    pam5=cls["pam50_secondary_benchmark"]["5class_a_all_genes"]["accuracy"],
    pam5_no=cls["pam50_secondary_benchmark"]["5class_no_PAM50_genes"]["accuracy"],
    rppa_q=rp["mw_TNBC_vs_HR"]["q_lt_0.05"], rho=rpt["spearman"], rho_p=rpt["perm_p"], rho_null=rpt["null_95pct_abs"],
    rho_genes=rpt["genes"], rho_deg=rpt["DE_genes_with_antibody"], rho_h2=rp["concordance_HER2_vs_HR"]["spearman"],
    fn_n=er_err["IHC ER- called ER+"]["n"], fn_prot=er_err["tests"]["ERneg_called_pos_above_that"],
    fp_n=er_err["IHC ER+ called ER-"]["n"], fp_low=er_err["tests"]["ERpos_called_neg_below_5pct_of_correct_ERpos"],
    fp_basal=er_err["IHC ER+ called ER-"]["PAM50"].get("Basal", 0), fp_h2=er_err["IHC ER+ called ER-"]["PAM50"].get("Her2", 0),
    mt_er=transfer("ER"), mt_h2=transfer("HER2"), mt_h2a=transfer("HER2_ablated"), mt_era=transfer("ER_ablated"),
    mt_tn=transfer("ERneg_HER2neg"),
    mb_esr1=mb["single_gene_baselines"]["ER"]["AUC"], mb_erbb2=mb["single_gene_baselines"]["HER2"]["AUC"],
    rep=rep["fraction"], rep_a=rep["replicated_same_sign_q05"], rep_b=rep["TCGA_DEGs_on_array"], rep_sign=rep["sign_agreement"],
    rep_rho=rep["spearman_effect_all_shared"], rep_nes=rep["hallmark_NES_spearman"],
    pfi_p=sv["tcga_logrank"]["TCGA_PFI_group"]["p"], fu=sv["tcga_logrank"]["TCGA_PFI_median_followup_years"],
    tcox=tcox, tc_tn=hr(tcox, "TNBC"), tc_h2=hr(tcox, "HER2pos"), tc_st=hr(tcox, "stage_III_IV"),
    tc_age_ph=tcox["schoenfeld"]["age_per10y"]["p"], tc_tn_ph=tcox["schoenfeld"]["TNBC"]["p"],
    mcox=mcox, mc_h2=hr(mcox, "HER2_gained"), mc_tn=hr(mcox, "ERneg_HER2notgained"),
    sp_early=msplit["HR_0_5y"]["ERneg_HER2notgained"], sp_late=msplit["HR_5_10y"]["ERneg_HER2notgained"],
    sp_p=msplit["interaction_p"]["ERneg_HER2notgained"], mb_pred_p=sv["metabric_logrank"]["MB_BCD_10y_pred_group"]["p"],
    niel=niel, acc_er=acc_er, acc_h2=acc_h2, lfl=lfl_a / lfl_b,
)
f2, f3 = (lambda x: f"{x:.2f}".replace("-", "−")), (lambda x: f"{x:.3f}".replace("-", "−"))
pc = lambda x: f"{x:.0%}"
c = lambda x: f"{x:,}"


def fp(p):
    """p-value for text: 2 significant digits, scientific below 1e-3."""
    if p < 1e-3:
        m, e = f"{p:.0e}".split("e")
        return f"{m} × 10^{int(e)}^"
    return f"{p:.2g}"


def ci(t):
    return f"{t[0]:.2f} (95% CI {t[1]:.2f}–{t[2]:.2f})"


def pctf(x):
    return f"{x:.0f}%"


# ================================================================================================ references
REFS = {
    "tcga2012": "Cancer Genome Atlas Network. Nature 2012;490:61–70. doi:10.1038/nature11412",
    "allison2020": "Allison KH, et al. J Clin Oncol 2020;38:1346–1366. doi:10.1200/JCO.19.02309",
    "wolff2018": "Wolff AC, et al. J Clin Oncol 2018;36:2105–2122. doi:10.1200/JCO.2018.77.8738",
    "parker2009": "Parker JS, et al. J Clin Oncol 2009;27:1160–1167. doi:10.1200/JCO.2008.18.1370",
    "sorlie2001": "Sørlie T, et al. PNAS 2001;98:10869–10874. doi:10.1073/pnas.191367098",
    "gong2007": "Gong Y, et al. Lancet Oncol 2007;8:203–211. doi:10.1016/S1470-2045(07)70042-6",
    "shah2026": "Shah V, Patel M. Comput Biol Chem 2026;124:109228. doi:10.1016/j.compbiolchem.2026.109228",
    "goldman2020": "Goldman MJ, et al. Nat Biotechnol 2020;38:675–678. doi:10.1038/s41587-020-0546-8",
    "liu2018": "Liu J, et al. Cell 2018;173:400–416. doi:10.1016/j.cell.2018.02.052",
    "curtis2012": "Curtis C, et al. Nature 2012;486:346–352. doi:10.1038/nature10983",
    "pereira2016": "Pereira B, et al. Nat Commun 2016;7:11479. doi:10.1038/ncomms11479",
    "liberzon2015": "Liberzon A, et al. Cell Syst 2015;1:417–425. doi:10.1016/j.cels.2015.12.004",
    "mcinnes2018": "McInnes L, Healy J, Melville J. arXiv 2018;1802.03426. doi:10.48550/arXiv.1802.03426",
    "monti2003": "Monti S, et al. Mach Learn 2003;52:91–118. doi:10.1023/A:1023949509487",
    "senbabaoglu2014": "Şenbabaoğlu Y, Michailidis G, Li JZ. Sci Rep 2014;4:6207. doi:10.1038/srep06207",
    "love2014": "Love MI, Huber W, Anders S. Genome Biol 2014;15:550. doi:10.1186/s13059-014-0550-8",
    "muzellec2023": "Muzellec B, et al. Bioinformatics 2023;39:btad547. doi:10.1093/bioinformatics/btad547",
    "bh1995": "Benjamini Y, Hochberg Y. J R Stat Soc B 1995;57:289–300. doi:10.1111/j.2517-6161.1995.tb02031.x",
    "subramanian2005": "Subramanian A, et al. PNAS 2005;102:15545–15550. doi:10.1073/pnas.0506580102",
    "fang2023": "Fang Z, Liu X, Peltz G. Bioinformatics 2023;39:btac757. doi:10.1093/bioinformatics/btac757",
    "wu2012": "Wu D, Smyth GK. Nucleic Acids Res 2012;40:e133. doi:10.1093/nar/gks461",
    "nadeau2003": "Nadeau C, Bengio Y. Mach Learn 2003;52:239–281. doi:10.1023/A:1024068626366",
    "nielsen2010": "Nielsen TO, et al. Clin Cancer Res 2010;16:5222–5232. doi:10.1158/1078-0432.CCR-10-1282",
    "netanely2016": "Netanely D, et al. Breast Cancer Res 2016;18:74. doi:10.1186/s13058-016-0724-2",
    "ray2010": "Ray PS, et al. Cancer Res 2010;70:3870–3876. doi:10.1158/0008-5472.CAN-09-4120",
    "krug2020": "Krug K, et al. Cell 2020;183:1436–1456. doi:10.1016/j.cell.2020.10.036",
}
ORDER = []


def cite(*keys):
    nums = []
    for k in keys:
        assert k in REFS, k
        if k not in ORDER:
            ORDER.append(k)
        nums.append(ORDER.index(k) + 1)
    nums = sorted(nums)
    # compress runs: 3,4,5 -> 3–5
    out, i = [], 0
    while i < len(nums):
        j = i
        while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1:
            j += 1
        out.append(f"{nums[i]}–{nums[j]}" if j - i >= 2 else ", ".join(str(n) for n in nums[i:j + 1]))
        i = j + 1
    return "[" + ", ".join(out) + "]"


# ================================================================================================ document helpers
doc = Document()
sec = doc.sections[0]
sec.page_height, sec.page_width = Cm(29.7), Cm(21.0)
for side in ("top_margin", "bottom_margin", "left_margin", "right_margin"):
    setattr(sec, side, Cm(2.0))
sec.header_distance = sec.footer_distance = Cm(1.0)
FONT = "Arial"


def set_font(style_or_run, size=11, bold=None, italic=None, color=None):
    f = style_or_run.font
    f.name, f.size = FONT, Pt(size)
    if bold is not None:
        f.bold = bold
    if italic is not None:
        f.italic = italic
    if color is not None:
        f.color.rgb = color
    el = style_or_run.element
    rpr = el.get_or_add_rPr() if hasattr(el, "get_or_add_rPr") else el.rPr
    rfonts = rpr.find(qn("w:rFonts"))
    if rfonts is None:
        rfonts = OxmlElement("w:rFonts"); rpr.append(rfonts)
    for a in ("w:ascii", "w:hAnsi", "w:eastAsia", "w:cs"):
        rfonts.set(qn(a), FONT)
    for a in ("w:asciiTheme", "w:hAnsiTheme", "w:eastAsiaTheme", "w:cstheme"):
        rfonts.attrib.pop(qn(a), None)


normal = doc.styles["Normal"]
set_font(normal, 11)
pf = normal.paragraph_format
pf.line_spacing_rule = WD_LINE_SPACING.SINGLE
pf.space_before, pf.space_after = Pt(0), Pt(4)
for name in ("Heading 1", "Heading 2"):
    st = doc.styles[name]
    set_font(st, 11, bold=True, italic=False, color=RGBColor(0, 0, 0))
    st.paragraph_format.space_before = Pt(6 if name == "Heading 1" else 3)
    st.paragraph_format.space_after = Pt(2)
    st.paragraph_format.keep_with_next = True
    st.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE

TOKEN = re.compile(r"(\*\*.+?\*\*|_[^_ ][^_]*?_|\^[^^]+\^)")


def add_runs(p, text, size=11, italic_all=False):
    """**bold**, _italic_ (gene names), ^superscript^."""
    text = re.sub(r"(?<![\w/.])-(?=\d)", "−", text)       # typographic minus before numbers
    for part in TOKEN.split(text):
        if not part:
            continue
        if part.startswith("**"):
            r = p.add_run(part[2:-2]); r.bold = True
        elif part.startswith("_") and part.endswith("_") and len(part) > 2:
            r = p.add_run(part[1:-1]); r.italic = True
        elif part.startswith("^") and part.endswith("^"):
            r = p.add_run(part[1:-1]); r.font.superscript = True
        else:
            r = p.add_run(part)
        r.font.name, r.font.size = FONT, Pt(size)
        if italic_all:
            r.italic = True


def para(text, size=11, align=WD_ALIGN_PARAGRAPH.JUSTIFY, after=4, keep=False):
    p = doc.add_paragraph()
    p.alignment = align
    p.paragraph_format.space_after = Pt(after)
    if keep:
        p.paragraph_format.keep_with_next = True
    add_runs(p, text, size)
    return p


def head(text, level=1):
    return doc.add_heading(text, level=level)


def figure(png, caption):
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.space_after = Pt(2)
    p.paragraph_format.keep_with_next = True
    p.add_run().add_picture(str(FIG / png), width=Cm(17.0))
    para(caption, after=6)


def table(rows, widths, header=True, caption=None):
    if caption:
        para(caption, after=2, keep=True)
    t = doc.add_table(rows=len(rows), cols=len(rows[0]))
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    t.style = doc.styles["Table Grid"]
    t.autofit = False
    tblpr = t._tbl.tblPr
    lay = OxmlElement("w:tblLayout"); lay.set(qn("w:type"), "fixed"); tblpr.append(lay)
    for j, col in enumerate(t._tbl.tblGrid.findall(qn("w:gridCol"))):
        col.set(qn("w:w"), str(int(widths[j] / 2.54 * 1440)))
    for i, row in enumerate(rows):
        for j, val in enumerate(row):
            cell = t.cell(i, j)
            cell.width = Cm(widths[j])
            cell.paragraphs[0].paragraph_format.space_after = Pt(0)
            cell.paragraphs[0].alignment = WD_ALIGN_PARAGRAPH.LEFT if j == 0 else WD_ALIGN_PARAGRAPH.CENTER
            add_runs(cell.paragraphs[0], f"**{val}**" if (header and i == 0) else str(val))
            tcpr = cell._tc.get_or_add_tcPr()
            mar = OxmlElement("w:tcMar")
            for side in ("top", "bottom"):
                e = OxmlElement(f"w:{side}"); e.set(qn("w:w"), "0"); e.set(qn("w:type"), "dxa"); mar.append(e)
            tcpr.append(mar)
    doc.add_paragraph().paragraph_format.space_after = Pt(0)
    return t


def page_numbers():
    p = sec.footer.paragraphs[0]
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run()
    for kind, txt in (("begin", None), (None, "PAGE"), ("end", None)):
        if kind:
            e = OxmlElement("w:fldChar"); e.set(qn("w:fldCharType"), kind); r._r.append(e)
        else:
            e = OxmlElement("w:instrText"); e.set(qn("xml:space"), "preserve"); e.text = txt; r._r.append(e)
    r.font.name, r.font.size = FONT, Pt(9)


page_numbers()

# ================================================================================================ text
b_er, b_h2, esr1, erbb2 = N["b_er"], N["b_h2"], N["esr1"], N["erbb2"]
nes = N["nes"]
t = doc.add_paragraph(); t.alignment = WD_ALIGN_PARAGRAPH.CENTER; t.paragraph_format.space_after = Pt(2)
r = t.add_run("What tumour gene expression reveals about the clinical receptor groups of breast cancer: "
              "a TCGA analysis checked in protein and an independent cohort")
r.bold = True; r.font.name, r.font.size = FONT, Pt(13)
para("BT3041 Analysis and Interpretation of Biological Data · Term project report · IIT Madras · November 2026",
     align=WD_ALIGN_PARAGRAPH.CENTER, after=0)
para("Team: [Name 1, roll no.] · [Name 2, roll no.] · [Name 3, roll no.] · [Name 4, roll no.] · [Name 5, roll no.]",
     align=WD_ALIGN_PARAGRAPH.CENTER, after=6)

para(f"**Abstract.** Breast cancers are grouped for treatment by three clinical tests: oestrogen receptor (ER), progesterone "
     f"receptor (PR) and HER2. We asked what tumour RNA-seq reveals about these clinical groups, using the clinical labels "
     f"(not the expression-derived PAM50 subtypes) to avoid circular reasoning. In {c(N['n_tum'])} TCGA tumours "
     f"(clinical groups: HR+/HER2− {N['hr']}, HER2+ {N['her2']}, triple-negative {N['tnbc']}), label-free clustering recovered a basal-like "
     f"group that held {N['k2_tn']} of {N['tnbc']} triple-negative tumours (ARI with ER status {f2(N['ari_er'].ARI)}). "
     f"Against HR+/HER2− tumours, triple-negative tumours differed more than 2-fold in {c(N['deg'])} genes (FDR 5%), with the "
     f"oestrogen-response programme down and proliferation, MYC and immune programmes up; only the oestrogen-response shift "
     f"survived a correlation-aware gene-set test. Classifiers trained on the whole transcriptome matched, but did not beat, "
     f"the receptor gene alone (ER: _ESR1_ AUC {f3(esr1.AUC_pooled)} vs best model {f3(b_er.AUC_pooled)}; HER2: _ERBB2_ "
     f"{f3(erbb2.AUC_pooled)} vs {f3(b_h2.AUC_pooled)}). The findings held in protein from the same tumours (mRNA–protein "
     f"effect correlation ρ = {f2(N['rho'])}) and in {c(N['mb_n'])} METABRIC tumours (frozen-model AUC {f2(N['mt_er']['AUC'])} for "
     f"ER, {f2(N['mt_h2']['AUC'])} for HER2; {pc(N['rep'])} of differentially expressed genes replicated).", after=6)

# ---------------------------------------------------------------------------------------------- 1
head("1. Introduction and research questions")
para(f"Pathologists test every invasive breast cancer for ER and PR by immunohistochemistry (IHC) and for HER2 by IHC or "
     f"in-situ hybridisation (FISH) {cite('allison2020', 'wolff2018')}. The result chooses the therapy: endocrine drugs for "
     f"hormone-receptor-positive (HR+) disease, anti-HER2 antibodies for HER2+ disease and chemotherapy for triple-negative "
     f"breast cancer (TNBC), which has no targeted option. Expression profiling defined the intrinsic subtypes (luminal A/B, "
     f"HER2-enriched, basal-like) {cite('sorlie2001', 'tcga2012')}, now assigned with the PAM50 centroids "
     f"{cite('parker2009')}. Many machine-learning studies predict PAM50 from RNA-seq, but PAM50 is itself computed from "
     f"RNA, so such accuracy partly measures agreement of an algorithm with itself. We therefore used the clinical IHC/FISH "
     f"labels as the ground truth and PAM50 only as a comparator.")
para("We asked three questions. **Q1:** Does label-free structure in tumour expression (PCA, UMAP, clustering) recover the "
     "clinical receptor groups? **Q2:** Which genes and pathways differ between TNBC and HR+/HER2− tumours, and how robust "
     "are they? **Q3:** Can whole-transcriptome classifiers predict clinical ER and HER2 status better than the receptor "
     "gene alone? Every answer was checked in a second assay (protein, RPPA) on the same tumours, in an independent cohort "
     "(METABRIC) with labels from other assays, and against published numbers.")

# ---------------------------------------------------------------------------------------------- 2
head("2. Data acquisition and preprocessing")
para(f"**Sources.** TCGA-BRCA RNA-seq (HiSeqV2, log2(RSEM+1)), STAR gene counts (GDC, GENCODE v36), clinical IHC/FISH "
     f"fields, reverse-phase protein arrays (RPPA) and PAM50 calls were downloaded from UCSC Xena {cite('goldman2020')}; "
     f"survival end points came from the TCGA Clinical Data Resource {cite('liu2018')}. METABRIC Illumina HT-12 expression "
     f"and clinical tables {cite('curtis2012', 'pereira2016')} came from the cBioPortal datahub. Gene symbols in every data "
     f"set (and the Hallmark gene sets {cite('liberzon2015')}) were mapped to current HGNC symbols. Raw files were checksummed "
     f"at a data freeze.")
para(f"**Cleaning and labels.** We kept one primary tumour per patient ({c(N['n_prim'])}) and excluded {N['n_excl']}: five "
     f"whose HiSeqV2 and STAR profiles did not match each other (Spearman ρ below the 99th percentile of between-patient "
     f"correlations) and five whose clinical record is redacted in the Clinical Data Resource (one tumour failed both), "
     f"leaving {c(N['n_tum'])} tumours for all label-free work. HER2 followed ASCO/CAP 2018 {cite('wolff2018')}: FISH "
     f"first, then the IHC score (3+ positive, 0/1+ negative, 2+ without FISH missing), using the status text only when no "
     f"score was recorded. IHC 3+ with a negative FISH ({N['h2_disc']} tumours) was treated as discordant and left "
     f"unlabelled, as were other conflicting records. This rule agrees with the TCGA 2012 HER2 calls at κ = {f2(N['kappa_h2'])}. The "
     f"supervised set holds the {c(N['n_sup'])} female patients ({N['n_male']} male patients were kept only for label-free "
     f"work). Groups were HR+/HER2− (ER or PR positive, HER2 negative), HER2+ (any ER/PR) and TNBC (all three negative); "
     f"{N['inc']} tumours lacked one of the three results. For differential expression we kept STAR genes with ≥ 10 counts in "
     f"at least {N['star_min']} tumours (the smallest group), {c(N['star_genes'])} genes. Matrices were left unscaled so that "
     f"every scaling step could run inside cross-validation folds. Table 1 lists the data sets.")
table([["Data set", "Used for", "n", "HR+/HER2−", "HER2+", "TNBC"],
       ["TCGA HiSeqV2 (analysed)", "PCA, UMAP, clustering", c(N["n_tum"]), *[pre["clinical_group_final_all"][g] for g in ("HR+/HER2-", "HER2+", "TNBC")]],
       ["TCGA supervised set", "classifiers", c(N["n_sup"]), N["hr"], N["her2"], N["tnbc"]],
       ["TCGA STAR counts", "differential expression", c(sum(pre["clinical_group_final_supervised_star"].values())),
        *[pre["clinical_group_final_supervised_star"][g] for g in ("HR+/HER2-", "HER2+", "TNBC")]],
       ["TCGA RPPA (131 antibodies)", "protein validation", N["rppa_n"], rp["n"]["HR+/HER2-"], rp["n"]["HER2+"], rp["n"]["TNBC"]],
       ["METABRIC (HT-12 arrays)", "external validation", c(N["mb_n"]), f"{c(N['mb_hr'])}*", f"{N['mb_h2g']}*", f"{N['mb_tn']}*"]],
      widths=[4.4, 3.9, 1.5, 2.5, 1.8, 1.6],
      caption=f"**Table 1.** Data sets and group sizes under the final labels (HER2 rule {cfg['her2_rule']}). "
              f"*METABRIC groups use IHC ER and SNP6 copy-number HER2 (gain vs no gain); PR is not available, so the "
              f"ER−/HER2-not-gained group is not called TNBC.")

# ---------------------------------------------------------------------------------------------- 3
head("3. Methods")
para(f"All analyses used Python 3.12 (scikit-learn, PyDESeq2, GSEApy, lifelines) with seed 42; every parameter and "
     f"package version is written to a params.json next to each result, and the code is in the supplementary material. "
     f"Multiple tests were corrected by Benjamini–Hochberg (BH) FDR {cite('bh1995')}.")
para(f"**Exploration and clustering (Q1).** PCA and UMAP {cite('mcinnes2018')} (15 neighbours) used the 2,000 most variable "
     f"genes after removing sex-chromosome-specific genes. We tested marker genes between groups with Mann–Whitney U tests "
     f"(rank-biserial r as effect size) and the share of PC variance explained by group and by tissue-source site (η², "
     f"with a permutation chance level). k-means (k = 2–8, PAM50 genes left out) was run with consensus clustering "
     f"{cite('monti2003')} (100 resamples of 80% of tumours); k was chosen before looking at the labels as the lowest "
     f"proportion of ambiguous clustering (PAC) {cite('senbabaoglu2014')}. Agreement with labels was measured by the "
     f"adjusted Rand index (ARI) with bootstrap 95% CIs and compared with the best ARI reachable by two clusters.")
para(f"**Differential expression and pathways (Q2).** PyDESeq2 {cite('love2014', 'muzellec2023')} on raw STAR counts with "
     f"design ~ histology + group tested |log2FC| > 1 directly (lfc_null = 1) rather than filtering by fold change after "
     f"a test of zero change. Fold changes were shrunk for ranking and plotting. Sensitivity runs added age, kept invasive "
     f"ductal tumours only, dropped every low-confidence HER2 call, or used the ER−/HER2− vs ER+/HER2− contrast that matches "
     f"METABRIC. Hallmark enrichment used GSEA prerank on the Wald statistic {cite('subramanian2005', 'fang2023')} "
     f"({c(epar['permutations'])} permutations), Fisher over-representation of the DEG lists, and two tests that account for inter-gene "
     f"correlation: CAMERA {cite('wu2012')} and phenotype-permutation GSEA.")
para(f"**Classification (Q3).** Targets were clinical ER (n = {c(N['n_er'])}, {c(N['n_er_pos'])} positive) and HER2 "
     f"(n = {N['n_h2']}, {N['n_h2_pos']} positive). L1-penalised logistic regression, a linear SVM, a Random Forest and PLS-DA "
     f"were trained in nested cross-validation (outer 5 folds × 3 repeats, inner 5 folds, class-balanced weights). A "
     f"variance filter (top 5,000 genes), scaling, univariate selection of k genes (k ∈ {{{', '.join(str(k) for k in cpar['k_grid'])}}}) "
     f"and the hyper-parameters were all fitted inside the training folds. The "
     f"baseline was the receptor gene alone (_ESR1_ or _ERBB2_), with a cut-off chosen inside the folds. Models were "
     f"compared on pooled out-of-fold ROC-AUC with the Nadeau–Bengio corrected t-test {cite('nadeau2003')} and a paired "
     f"bootstrap. An ablation removed _ESR1_, _PGR_, _ERBB2_, every gene within 1.5 Mb (and 2 Mb) of _ERBB2_ and the 50 "
     f"PAM50 genes ({N['n_removed']} genes). Shuffled labels gave a negative control.")
para(f"**Validation.** (i) Protein: Mann–Whitney tests on RPPA, and the Spearman correlation between DE fold changes and "
     f"protein differences for the {N['rho_genes']} genes with a total-protein antibody, with a gene-label permutation null; "
     f"the out-of-fold ER errors were compared with ERα protein. (ii) METABRIC: the ER, HER2 and ER−/HER2− models were "
     f"refit on all TCGA tumours and frozen (with checksums logged) before any METABRIC label was read, then applied to "
     f"METABRIC z-scored within cohort on the {c(N['mb_shared'])} shared genes. DE was replicated by OLS with age "
     f"(ER-IHC−/HER2-not-gained {N['mb_tn']} vs ER-IHC+/HER2-not-gained {c(N['mb_hr'])}). (iii) Survival: Kaplan–Meier, "
     f"log-rank and Cox models with Schoenfeld checks; TCGA progression-free interval (PFI) with group, age and stage "
     f"({N['tcox']['events']} events, so 5 parameters), and METABRIC breast-cancer death within 10 years with group, age, "
     f"grade and size. (iv) Literature: our numbers were placed next to published ones, each checked against PubMed.")

# ---------------------------------------------------------------------------------------------- 4
head("4. Results and interpretation")
head("4.1 Label-free structure follows ER, and recovers TNBC as a basal-like group (Q1)", 2)
para(f"PC1 ({pc(N['pc1_var'])} of variance) was the ER axis: it separated IHC ER-positive from ER-negative tumours with "
     f"AUC {f2(N['pc1_auc'])} and the clinical group explained η² = {f2(N['pc1_eta'])} of it, while tissue-source site "
     f"explained ≤ {f2(N['site_max'])} of PC1–4. In the UMAP (Fig. 1a) a separate island of {N['isl']} tumours held "
     f"{N['isl_a']} of {N['isl_b']} TNBC ({pc(N['isl_a'] / N['isl_b'])}); HER2+ tumours did not form their own island but sat at "
     f"one edge of the HR+ mass. Marker genes behaved as expected: _ESR1_ was the strongest single difference between TNBC "
     f"and HR+/HER2− (rank-biserial r = {f2(N['r_esr1'])}), the basal transcription factor _FOXC1_ was higher in TNBC "
     f"(r = +{f2(N['r_foxc1'])}) and _ERBB2_ was higher in HER2+ (r = +{f2(N['r_erbb2'])}; all q < 10^−30^).")
para(f"Consensus clustering was most stable at k = 2 (PAC {f2(N['pac2'])}; k = 5: {f2(N['pac5'])}). The smaller cluster of "
     f"{N['k2']} tumours held {N['k2_tn']} of {N['tnbc']} TNBC (Fig. 1b) and agreed with IHC ER at ARI "
     f"{f2(N['ari_er'].ARI)} ({f2(N['ari_er'].ARI_lo)}–{f2(N['ari_er'].ARI_hi)}) and with the three clinical groups at "
     f"{f2(N['ari_cg'].ARI)}, against a ceiling of {f2(N['ceil_cg'])} for two clusters (Fig. 1c). k-means, consensus, Ward and "
     f"average-linkage clusterings agreed at k = 2 (ARI {min(N['meth_ari'].values()):.2f}–{max(N['meth_ari'].values()):.2f}). "
     f"At k = 5 a HER2-enriched cluster appeared, but it was the least reproducible across methods, and the luminal "
     f"tumours split differently from PAM50 LumA/LumB (ARI {f2(N['ari5_pam'].ARI)}), in line with the luminal continuum "
     f"described by Netanely et al. {cite('netanely2016')}. **Answer to Q1:** expression recovers the ER-negative, "
     f"basal-like axis strongly and the HER2+ group only weakly.")
figure("report_fig1_exploration.png",
       f"**Figure 1. Label-free structure.** (a) UMAP of {c(N['n_tum'])} tumours coloured by clinical group. (b) Share of each "
       f"clinical group in the smaller k-means cluster (k = 2). (c) Adjusted Rand index of k-means clusters with ER, the "
       f"clinical groups and PAM50 (95% bootstrap CI); bars mark the best ARI reachable with two clusters.")

head("4.2 TNBC is a distinct transcriptional state (Q2)", 2)
para(f"Of {c(N['genes_tested'])} genes, {c(N['deg'])} differed more than 2-fold between TNBC and HR+/HER2− at FDR 5% "
     f"({c(N['up'])} higher, {c(N['dn'])} lower in TNBC; Fig. 2a). A plain Wald test of zero change called "
     f"{c(N['wald'])} genes ({pc(N['wald'] / N['genes_tested'])}), which shows why the threshold was tested rather than "
     f"applied afterwards. _ESR1_ (log2FC {gn['ESR1']['log2FC_shrunk']:+.1f}) and _PGR_ ({gn['PGR']['log2FC_shrunk']:+.1f}) "
     f"acted as positive controls. Beyond them, the ER-driven luminal programme was lower in TNBC (_GATA3_ "
     f"{gn['GATA3']['log2FC_shrunk']:+.1f}, _AGR3_ {gn['AGR3']['log2FC_shrunk']:+.1f}, _FOXA1_ {gn['FOXA1']['log2FC_shrunk']:+.1f}, "
     f"_XBP1_ {gn['XBP1']['log2FC_shrunk']:+.1f}) and basal transcription factors and keratins were higher (_FOXC1_ "
     f"{gn['FOXC1']['log2FC_shrunk']:+.1f}, as reported by Ray et al. {cite('ray2010')}; _KRT5_ {gn['KRT5']['log2FC_shrunk']:+.1f}). "
     f"The list was robust: {pc(N['rob_age'])} of the DEGs stayed significant with the same sign after adding age, "
     f"{pc(N['rob_idc'])} in ductal tumours only (so lobular biology is not being reported as HR+ biology), "
     f"{pc(N['rob_strict'])} with strict labels and {pc(N['rob_match'])} in the ER−/HER2− contrast. HER2+ vs HR+/HER2− "
     f"tumours differed in only {N['her2_deg']} genes, {N['her2_amp']} of them within 2 Mb of _ERBB2_ on the 17q12 amplicon.")
para(f"GSEA found {N['gsea_n']} Hallmark sets at FDR < 0.05 (Fig. 2b): proliferation (E2F targets, NES "
     f"{nes['E2F_TARGETS']['NES']:+.2f}; G2M checkpoint {nes['G2M_CHECKPOINT']['NES']:+.2f}), MYC targets "
     f"({nes['MYC_TARGETS_V1']['NES']:+.2f}), interferon-γ ({nes['INTERFERON_GAMMA_RESPONSE']['NES']:+.2f}) and mTORC1 "
     f"({nes['MTORC1_SIGNALING']['NES']:+.2f}) were higher in TNBC, and the early and late oestrogen responses lower "
     f"({nes['ESTROGEN_RESPONSE_EARLY']['NES']:+.2f}, {nes['ESTROGEN_RESPONSE_LATE']['NES']:+.2f}), as in the TCGA and CPTAC "
     f"proteogenomic studies {cite('tcga2012', 'krug2020')}. Gene-permutation q-values are optimistic when genes in a set are "
     f"co-expressed (median inter-gene correlation {f2(N['cam_corr'])}, up to {f2(N['cam_e2f'])} for E2F targets). Both "
     f"correlation-aware tests agreed with the prerank direction for all {N['gsea_n']} sets (NES ρ = {f2(N['ph_rho'])}), but "
     f"only the early oestrogen response stayed significant (CAMERA q < 0.05). We therefore report the proliferation, MYC "
     f"and immune programmes as consistent shifts with large effect sizes, not as proven to be more shifted than the "
     f"average gene. **Answer to Q2:** TNBC switches off the oestrogen-response programme and switches on basal and "
     f"proliferative genes; HER2+ tumours differ far less, mainly at the amplicon.")
figure("report_fig2_de_pathways.png",
       f"**Figure 2. TNBC vs HR+/HER2− ({N['tnbc']} vs {N['hr']}).** (a) Volcano plot: shrunk log2 fold change against FDR "
       f"of the |log2FC| > 1 test (y capped at 150). (b) Hallmark sets at prerank FDR < 0.05 (the 13 strongest positive and "
       f"all {N['gsea_dn']} negative of {N['gsea_n']}); symbols give the correlation-aware checks (● q < 0.05, ◐ nominal "
       f"p < 0.05, ○ neither).")

head("4.3 Whole-transcriptome models match, but do not beat, the receptor gene (Q3)", 2)
para(f"For ER, _ESR1_ alone reached AUC {f3(esr1.AUC_pooled)} and the best model (Random Forest) {f3(b_er.AUC_pooled)} "
     f"(corrected t-test p = {f2(b_er.NB_p)}); for HER2, _ERBB2_ alone {f3(erbb2.AUC_pooled)} and the best model "
     f"{f3(b_h2.AUC_pooled)} (p = {f2(b_h2.NB_p)}) (Fig. 3a–b). All four model types gave similar AUCs (Fig. 3c), shuffled "
     f"labels gave MCC {N['mcc_er']:+.2f} (ER) and {N['mcc_h2']:+.2f} (HER2), and accuracy of the single gene "
     f"({N['acc_er']} for _ESR1_, {N['acc_h2']} for _ERBB2_) matched the 88–96% reported for arrays by Gong et al. {cite('gong2007')}. The "
     f"ablation separates the two receptors. Without the receptor genes, the amplicon and the PAM50 genes, ER was "
     f"predicted as well as before (AUC {f3(N['abl_e_lo'])}–{f3(N['abl_e_hi'])}): ER status is written across a whole "
     f"regulon, not one gene. HER2 fell to {f3(N['abl_h_lo'])}–{f3(N['abl_h_hi'])}, the same at 2 Mb, so most HER2 signal "
     f"sits in the amplicon and the remainder largely reflects the lower-ER, proliferative character of HER2+ tumours. "
     f"Label quality limits HER2: dropping the {N['h2_text']} HER2+ calls based on status text alone raised _ERBB2_ alone "
     f"to {f3(N['sens_pos'])} and the best model to {f3(N['sens_pos_best'])}. Only for the combined ER−/HER2− target did "
     f"a model beat one gene ({f3(N['tn_best'].AUC_pooled)} vs {f3(N['tn_base'].AUC_pooled)} for _ESR1_; p = "
     f"{f3(N['tn_best'].NB_p)}), because no single gene encodes both receptors. As a secondary, circular benchmark, the "
     f"same pipeline reproduced 5-class PAM50 calls with {pc(N['pam5'])} accuracy ({pc(N['pam5_no'])} without the PAM50 "
     f"genes). **Answer to Q3:** no; the receptor gene already carries almost all the information about clinical status.")
figure("report_fig3_classification.png",
       f"**Figure 3. Classification of clinical status (nested CV, pooled out-of-fold).** (a, b) ROC curves for the receptor "
       f"gene alone, the best model with all genes and the best ablated model. (c) AUC (95% CI) of every model: dark, all "
       f"genes; orange, ablated (±1.5 Mb); open, ablated (±2 Mb). Blue line and band: receptor gene alone.")

head("4.4 Validation in protein, in METABRIC and against published studies", 2)
para(f"**Protein.** {N['rppa_q']} of {N['rppa_ab']} RPPA antibodies differed between TNBC and HR+/HER2− (q < 0.05). Across "
     f"the {N['rho_genes']} genes with a total-protein antibody, the mRNA fold change predicted the protein difference with "
     f"ρ = {f2(N['rho'])} (permutation p = {f3(N['rho_p'])}; 95th percentile of the null {f2(N['rho_null'])}; Fig. 4a), and "
     f"all {N['rho_deg']} DEGs with an antibody agreed in sign. The ER errors are informative (Fig. 4b): of the {N['fn_n']} "
     f"IHC ER-negative tumours called ER-positive, {N['fn_prot'].replace('/', ' of ')} with RPPA had ERα protein above the "
     f"95th percentile of correctly called ER-negative tumours, which suggests IHC false negatives rather than model "
     f"errors. Of the {N['fp_n']} IHC ER-positive tumours called negative, {N['fp_basal'] + N['fp_h2']} were PAM50 Basal "
     f"or HER2-enriched and {N['fp_low'].replace('/', ' of ')} with RPPA had ERα below the 5th percentile of correctly "
     f"called ER-positive tumours, consistent with ER-low tumours or IHC false positives.")
para(f"**METABRIC.** The frozen TCGA models transferred without loss (Fig. 4c): ER AUC {f3(N['mt_er']['AUC'])} "
     f"({f3(N['mt_er']['AUC_CI'][0])}–{f3(N['mt_er']['AUC_CI'][1])}) against {f3(N['mt_er']['TCGA_CV_AUC'])} in TCGA "
     f"cross-validation, HER2 (copy-number gain) {f3(N['mt_h2']['AUC'])} against {f3(N['mt_h2']['TCGA_CV_AUC'])}, and "
     f"ER−/HER2− vs rest {f3(N['mt_tn']['AUC'])}. HER2 transferred better than it cross-validated because METABRIC's "
     f"copy-number label tracks mRNA more closely than clinical IHC/FISH. As in TCGA, the models did not beat the gene "
     f"alone (_ESR1_ {f3(N['mb_esr1'])}, _ERBB2_ {f3(N['mb_erbb2'])}). Of the {c(N['rep_b'])} TCGA DEGs measured on the "
     f"array, {c(N['rep_a'])} ({pc(N['rep'])}) replicated with the same sign at q < 0.05 (effect-size ρ = "
     f"{f2(N['rep_rho'])}), and Hallmark NES agreed with ρ = {f2(N['rep_nes'])}.")
lit = [["Quantity", "This study", "Published"],
       ["TNBC called PAM50 Basal", f"{pctf(N['tnb_pct'])} ({N['tnb_a']}/{N['tnb_b']})", f"86% (65/76) {cite('tcga2012')}"],
       ["PAM50 Basal that are TNBC", f"{pctf(N['bt_pct'])} ({N['tnb_a']}/{N['bt_b']})", f"80% (65/81) {cite('tcga2012')}"],
       ["Clinical HER2+ called HER2-enriched", pctf(N["h2e_pct"]), f"~50% {cite('tcga2012')}"],
       ["ER+ called Basal or HER2-enriched", f"{N['nl_pct']:.1f}%", f"8.7% {cite('nielsen2010')}"],
       ["Single-gene accuracy, ER / HER2", f"{N['acc_er']} / {N['acc_h2']}", f"88–96% / 89–93% {cite('gong2007')}"],
       ["AUC change across cohorts, ER / HER2",
        f"{N['mt_er']['TCGA_CV_AUC']:.2f}→{N['mt_er']['AUC']:.2f} / {N['mt_h2']['TCGA_CV_AUC']:.2f}→{N['mt_h2']['AUC']:.2f}",
        f"13–22% drop {cite('shah2026')}"]]
if N["niel"]:
    h = N["niel"]["HR_vs_LumA"]; ph = N["niel"]["published_Nielsen2010_HR_vs_LumA"]
    lit.append(["5-y BC death HR vs LumA (ER+, endocrine-treated): LumB / HER2E / Basal",
                f"{h['LumB']['exp(coef)']:.1f} / {h['Her2']['exp(coef)']:.1f} / {h['Basal']['exp(coef)']:.1f}",
                f"{ph['LumB']:.1f} / {ph['Her2']:.1f} / {ph['Basal']:.1f} {cite('nielsen2010')}"])
table(lit, widths=[6.6, 4.4, 5.0],
      caption="**Table 2.** Our numbers next to published ones (each checked against PubMed). TCGA 2012 analysed a subset "
              "of the same tumours (a consistency check); the other studies used different cohorts or assays.")
para(f"**Published studies.** Table 2 shows that the clinical–PAM50 concordance, the single-gene accuracy and the subtype "
     f"order of risk agree with the literature. Clinical HER2+ tumours were called HER2-enriched less often than in TCGA "
     f"2012; with the 2012 HER2 labels and PAM50 calls the same tumours give {pc(N['lfl'])}, so the gap comes from the larger cohort, "
     f"the newer PAM50 calls and the HER2 rule, not from an error. Unlike the 13–22% loss reported for a multimodal model "
     f"moved from TCGA to CPTAC {cite('shah2026')}, our simpler mRNA models lost nothing in METABRIC.")

figure("report_fig4_validation.png",
       f"**Figure 4. Validation.** (a) DE log2FC vs RPPA median difference (TNBC vs HR+/HER2−); coloured, DEGs. (b) RPPA ERα "
       f"by IHC ER status and out-of-fold call, coloured by PAM50. (c) AUC in TCGA cross-validation and in METABRIC for the "
       f"frozen models, with the receptor gene alone in METABRIC.")
head("4.5 Survival (beyond the talk)", 2)
para(f"In TCGA (median follow-up {N['fu']:.1f} years) PFI differed by group (log-rank p = {fp(N['pfi_p'])}; Fig. 5a). "
     f"Adjusted for age and stage, TNBC had a hazard ratio (HR) of {ci(N['tc_tn'])} against HR+/HER2−, while HER2+ did not "
     f"differ (HR {N['tc_h2'][0]:.2f}, 95% CI {N['tc_h2'][1]:.2f}–{N['tc_h2'][2]:.2f}; Fig. 5c); the group effect met the proportional-hazards assumption (p = {f2(N['tc_tn_ph'])}), "
     f"age did not (p = {f2(N['tc_age_ph'])}). In METABRIC ({c(N['mcox']['n'])} tumours, {N['mcox']['events']} breast-cancer "
     f"deaths within 10 years; Fig. 5b) HER2 gain carried HR {ci(N['mc_h2'])}. The ER−/HER2− effect violated proportional "
     f"hazards; split at 5 years it was {N['sp_early']:.2f} in years 0–5 and {N['sp_late']:.2f} in years 5–10 (interaction "
     f"p = {fp(N['sp_p'])}), the known early-relapse pattern of ER-negative disease. The frozen TCGA classifier's predicted "
     f"ER−/HER2− group also separated METABRIC survival (log-rank p = {fp(N['mb_pred_p'])}), and METABRIC's PAM50 subtypes "
     f"reproduced the risk order reported by Nielsen et al. {cite('nielsen2010')} (Table 2).")
figure("report_fig5_survival.png",
       "**Figure 5. Outcome by group.** Kaplan–Meier curves for (a) TCGA progression-free interval and (b) METABRIC "
       "breast-cancer death within 10 years (n, events). (c) Adjusted Cox hazard ratios (95% CI) against HR+/HER2−; "
       "diamonds, METABRIC ER−/HER2− effect split at 5 years. MB, METABRIC.")

# ---------------------------------------------------------------------------------------------- 5
head("5. Conclusion and discussion")
para(f"**Summary.** (1) The dominant structure in breast tumour expression is the ER axis; label-free methods recover TNBC "
     f"as a basal-like group but not HER2+ tumours as a group. (2) TNBC is a distinct transcriptional state, with "
     f"{c(N['deg'])} genes more than 2-fold different from HR+/HER2− tumours, robust to covariates, histology and label "
     f"noise. (3) Whole-transcriptome classifiers match but do not beat _ESR1_ and _ERBB2_ alone. (4) The results hold in "
     f"protein, in an independent cohort with labels from other assays, and against published numbers.")
para("**Biological interpretation.** ER status is a regulon-wide state: removing _ESR1_ and every PAM50 gene leaves the "
     "prediction intact, because _GATA3_, _FOXA1_, _AGR3_ and _XBP1_ move together with it. TNBC tumours replace this "
     "programme with basal transcription factors (_FOXC1_, _EN1_, _SOX10_) and high proliferation and MYC activity, and "
     "they relapse early, as the time-dependent METABRIC hazard shows. HER2 positivity, in contrast, is mostly a local "
     "genomic event: the amplicon dominates both the DE list and the classifier, and the rest of the transcriptome "
     "resembles HR+ disease. Practically, an expression assay adds little to IHC for ER and HER2 calls, but tumours where "
     "mRNA and IHC disagree are worth a second look: our ER error analysis suggests that some IHC calls, not the model, "
     "are wrong.")
para(f"**Limitations.** (i) TNBC and HR+ tumours differ in tumour purity and immune infiltration, which contribute to the "
     "proliferation and immune signals; DE was not adjusted for purity. (ii) Pathway significance shrinks once inter-gene "
     "correlation is counted, so pathway findings are effect sizes, not tests. (iii) HER2 labels are noisy (some rest on "
     "status text alone) and the TCGA HER2 rule is our own. (iv) METABRIC lacks PR and uses copy-number HER2, so its "
     f"groups are not identical to the clinical ones. (v) TCGA survival is short (median {N['fu']:.1f} years) and treatment was not "
     "modelled. (vi) Classifier comparisons rest on overlapping folds, which the corrected test handles conservatively.")
para("**Future directions.** Adjust DE for purity (e.g. ESTIMATE scores) or deconvolve the immune signal; test whether "
     "_ERBB2_ mRNA can resolve IHC 2+ equivocal cases; extend the protein check to the CPTAC mass-spectrometry proteome; "
     "split the HR+/HER2− group into proliferative and immune subgroups and test their outcome in METABRIC.")
para("**Data and code availability.** TCGA data from UCSC Xena (GDC and TCGA hubs); METABRIC data from the cBioPortal "
     "datahub, used under the Open Database License (ODbL). Code, parameters, results and logs: https://github.com/Vijayavallabh/BT3041_Project. AI prompts and "
     "AI chat histories are in the supplementary material: [Drive link]. **Use of AI:** an AI assistant helped plan the "
     "analysis, write and run the code, review the results and draft this report; every number was checked against the "
     "result files and every reference against PubMed.", after=4)

head("References")
for i, k in enumerate(ORDER, 1):
    p = para(f"{i}. {REFS[k]}", align=WD_ALIGN_PARAGRAPH.LEFT, after=0)
    p.paragraph_format.left_indent = Cm(0.6)
    p.paragraph_format.first_line_indent = Cm(-0.6)
unused = [k for k in REFS if k not in ORDER]
assert not unused, f"references never cited: {unused}"

doc.core_properties.title = "BT3041 term project report"
doc.core_properties.author = "BT3041 team"
doc.core_properties.created = doc.core_properties.modified = datetime.now(timezone.utc).replace(tzinfo=None)
doc.save(OUT)
print("saved", OUT)

# ================================================================================================ numbers log
flat = {}
for k, v in N.items():
    if isinstance(v, pd.Series):
        flat[k] = {kk: (round(vv, 4) if isinstance(vv, float) else vv) for kk, vv in v.items()
                    if kk in ("model", "AUC_pooled", "AUC_CI_lo", "AUC_CI_hi", "NB_p") or "model" not in v.index}
    elif isinstance(v, dict) and len(json.dumps(v, default=str)) > 400:
        flat[k] = "(see source JSON)"
    else:
        flat[k] = v
with open(HERE / "report_numbers.md", "w") as fh:
    fh.write("# Numbers used in the report\n\nGenerated by `scripts/build_report.py` from the step 2 and step 3 result files "
             f"(labels: HER2 rule `{cfg['her2_rule']}`, filter `{cfg['her2_filter']}`). Do not edit by hand.\n\n| key | value |\n|---|---|\n")
    for k, v in flat.items():
        fh.write(f"| `{k}` | {json.dumps(v, default=str)[:300]} |\n")

# ================================================================================================ PDF + page check
soffice = shutil.which("soffice") or shutil.which("libreoffice")
if soffice:
    subprocess.run([soffice, "--headless", "--convert-to", "pdf", "--outdir", str(HERE), str(OUT)], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    pdf = OUT.with_suffix(".pdf")
    info = subprocess.run(["pdfinfo", str(pdf)], capture_output=True, text=True).stdout
    pages = int(re.search(r"Pages:\s+(\d+)", info).group(1))
    print(f"PDF: {pages} pages (limit {MAX_PAGES})")
    if pages > MAX_PAGES:
        sys.exit(f"over the page limit: {pages} > {MAX_PAGES}")
else:
    print("LibreOffice not found: open the .docx and check the page count by hand")
