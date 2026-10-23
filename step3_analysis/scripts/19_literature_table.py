"""Step 3.9: literature comparison table (plan section 3F.1): every headline number next to its published counterpart.

Our numbers are read from the step 2 and step 3 reports (none typed by hand). Published numbers were checked against
the paper or its PubMed abstract on 21 Oct 2026 (DOIs below); `comparability` says how like-for-like each pair is.
Writes results/literature_table.tsv and results/literature_table.md.
"""
import sys, json
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.common3 import RES3, RES2  # noqa: E402

en = json.load(open(RES3 / "enrichment" / "enrichment_report.json"))
cl = pd.read_csv(RES3 / "classify" / "summary.tsv", sep="\t")
cr = json.load(open(RES3 / "classify" / "classify_report.json"))
oof = pd.read_parquet(RES3 / "classify" / "oof_predictions.parquet")
mb = json.load(open(RES3 / "transfer" / "metabric_report.json"))
rp = json.load(open(RES3 / "rppa" / "rppa_report.json"))
sv = json.load(open(RES3 / "survival" / "survival_report.json"))
ex2 = json.load(open(RES2 / "explore" / "explore_report.json"))
main = cl[cl.tag == "main"]


def acc_single(target):
    o = oof[(oof.tag == "main") & (oof.target == target) & oof.model.str.startswith("single")]
    return (o.y == o.pred).mean()


def best(target, abl="a_all_genes"):
    d = main[(main.target == target) & (main.ablation == abl) & ~main.model.str.startswith("single")]
    return d.sort_values("AUC_pooled", ascending=False).iloc[0]


ca = {d["comparison"]: d for d in en["cross_assay_agreement"]}
tn = ca["TNBC vs PAM50 Basal"]; h2 = ca["clinical HER2+ vs PAM50 HER2-enriched"]; erl = ca["ER+ tumours called non-luminal (Basal or Her2)"]
g = en["gsea_main"]
tr = {(r["model"], r["scaling"]): r for r in mb["transfer"]}
t_er, t_h = tr[("ER", "z_within_cohort")], tr[("HER2", "z_within_cohort")]
drop_er = (t_er["TCGA_CV_AUC"] - t_er["AUC"]) / t_er["TCGA_CV_AUC"]
drop_h = (t_h["TCGA_CV_AUC"] - t_h["AUC"]) / t_h["TCGA_CV_AUC"]
ni = sv["nielsen_direction_check"]["HR_vs_LumA"]
tc = [x for x in sv["tcga_cox"] if x["label"] == "TCGA Cox PFI"][0]
pam5 = cr["pam50_secondary_benchmark"].get("5class_a_all_genes", {})
ber, bh2 = best("ER"), best("HER2")

rows = [
    ("Clinical groups vs PAM50", "TNBC called PAM50 Basal", f"{tn['pct_clinical_called_PAM50']:.0f}% ({tn['both']}/{tn['both'] + tn['only_clinical']})",
     "86% (65/76)", "TCGA 2012", "10.1038/nature11412", "same cohort (subset), RNA-seq vs microarray PAM50"),
    ("Clinical groups vs PAM50", "PAM50 Basal that are TNBC", f"{tn['pct_PAM50_with_clinical']:.0f}% ({tn['both']}/{tn['both'] + tn['only_PAM50']})",
     "80% (65/81)", "TCGA 2012", "10.1038/nature11412", "same cohort (subset)"),
    ("Clinical groups vs PAM50", "Clinical HER2+ called HER2-enriched", f"{h2['pct_clinical_called_PAM50']:.0f}% ({h2['both']}/{h2['both'] + h2['only_clinical']})",
     "~50%", "TCGA 2012", "10.1038/nature11412", "gap explained in step 2 README (cohort size, PAM50 call, HER2 rule)"),
    ("Clinical groups vs PAM50", "ER+ tumours called non-luminal (Basal/HER2E)", f"{erl['pct']:.1f}% ({erl['pct_incl_Normal_like']:.1f}% with Normal-like)",
     "8.7% (10.7%)", "Nielsen 2010", "10.1158/1078-0432.CCR-10-1282", "different cohort (786 ER+ tamoxifen-treated, qRT-PCR PAM50)"),
    ("DE", "FOXC1 higher in TNBC/basal", f"log2FC {en['named_genes_main']['FOXC1']['log2FC_shrunk']:+.2f}, DE",
     "higher in basal-like", "Ray 2010", "10.1158/0008-5472.CAN-09-4120", "direction only"),
    ("Pathways", "E2F targets / G2M up in TNBC", f"NES {g['E2F_TARGETS']['NES']:+.2f} / {g['G2M_CHECKPOINT']['NES']:+.2f}",
     "up in basal (proteogenomic)", "Krug 2020", "10.1016/j.cell.2020.10.036", "direction only; different cohort (CPTAC)"),
    ("Pathways", "MYC targets up in TNBC", f"NES {g['MYC_TARGETS_V1']['NES']:+.2f}", "core MYC signature prominent in basal-like",
     "Chandriani 2009", "10.1371/journal.pone.0006693", "direction only"),
    ("Pathways", "Estrogen response down in TNBC", f"NES {g['ESTROGEN_RESPONSE_EARLY']['NES']:+.2f} (early), {g['ESTROGEN_RESPONSE_LATE']['NES']:+.2f} (late)",
     "ER/FOXA1/GATA3 luminal hub", "TCGA 2012; Liberzon 2015", "10.1038/nature11412; 10.1016/j.cels.2015.12.004", "direction only"),
    ("Pathways", "Correlation-adjusted enrichment", f"CAMERA: {en['camera_check']['sets_q_lt_0.05']} hallmark(s) at q < 0.05 (median inter-gene ρ {en['camera_check']['median_inter_gene_corr']:.2f})",
     "gene-permutation tests are very sensitive to inter-gene correlation", "Wu & Smyth 2012", "10.1093/nar/gks461", "method"),
    ("Classification", "Single receptor gene accuracy, ER (ESR1)", f"{acc_single('ER'):.0%}", "88–96%", "Gong 2007", "10.1016/S1470-2045(07)70042-6",
     "Affymetrix arrays, fixed cut-off; ours: cut-off by Youden inside CV folds"),
    ("Classification", "Single receptor gene accuracy, HER2 (ERBB2)", f"{acc_single('HER2'):.0%}", "89–93%", "Gong 2007", "10.1016/S1470-2045(07)70042-6", "as above"),
    ("Classification", "ER prediction MCC", f"{ber.MCC:.2f} ({ber.model}, all genes)", "0.80–0.83", "Chierici 2020", "10.3389/fonc.2020.01065",
     "TCGA BRCA-ER (n 381), multi-omics; ours mRNA only, n 1,024"),
    ("Classification", "In-distribution AUC, ER / HER2", f"{ber.AUC_pooled:.3f} / {bh2.AUC_pooled:.3f}", "0.99 / 0.98", "Shah & Patel 2026",
     "10.1016/j.compbiolchem.2026.109228", "multimodal (WSI + omics) upper bound"),
    ("Transfer", "Cross-cohort AUC change, ER / HER2", f"{t_er['TCGA_CV_AUC']:.3f}→{t_er['AUC']:.3f} / {t_h['TCGA_CV_AUC']:.3f}→{t_h['AUC']:.3f} (TCGA CV→METABRIC; no drop)",
     "13–22% drop (TCGA→CPTAC)", "Shah & Patel 2026",
     "10.1016/j.compbiolchem.2026.109228", "different target cohort and label assays"),
    ("Secondary (circular)", "5-class PAM50 accuracy from RNA-seq", f"{pam5.get('accuracy', np.nan):.1%} (L1 logistic, CV)", "88.42%", "Mostavi 2020",
     "10.1186/s12920-020-0677-2", "TCGA, CNN; agreement with a centroid rule, not ground truth"),
    ("Secondary (circular)", "Concordance with PAM50 calls, external tests", f"{pam5.get('accuracy', np.nan):.1%}", "~90%", "Cascianelli 2020",
     "10.1038/s41598-020-70832-2", "regularised multiclass logistic regression"),
    ("Protein", "ESR1 / ERBB2 mRNA vs protein (Spearman, step 2)",
     f"{ex2['mrna_protein_spearman']['selected']['ERALPHA']:.2f} / {ex2['mrna_protein_spearman']['selected']['HER2']:.2f} (RPPA)",
     "0.62–0.77 (mRNA vs clinical IHC/FISH status)", "Gong 2007", "10.1016/S1470-2045(07)70042-6", "protein assay vs clinical status; indicative"),
    ("Protein", "DE vs RPPA effect concordance", f"ρ = {rp['concordance_TNBC_vs_HR']['spearman']:.2f} (perm. p {rp['concordance_TNBC_vs_HR']['perm_p']:.3f})",
     "basal vs luminal separation at proteome level", "Mertins 2016", "10.1038/nature18003", "direction only"),
    ("Survival", "TCGA TNBC vs HR+ PFI hazard ratio", f"{tc['HR']['TNBC']['exp(coef)']:.2f} ({tc['HR']['TNBC']['exp(coef) lower 95%']:.2f}–{tc['HR']['TNBC']['exp(coef) upper 95%']:.2f})",
     "use PFI; OS needs caution in BRCA (short follow-up)", "Liu 2018 (TCGA-CDR)", "10.1016/j.cell.2018.02.052", "guidance"),
    ("Survival", "HR vs LumA, ER+ endocrine-treated, 5-y BC death: LumB / HER2E / Basal",
     f"{ni['LumB']['exp(coef)']:.2f} / {ni['Her2']['exp(coef)']:.2f} / {ni['Basal']['exp(coef)']:.2f} (METABRIC)", "1.99 / 3.65 / 17.71", "Nielsen 2010",
     "10.1158/1078-0432.CCR-10-1282", "direction and order; METABRIC published PAM50 calls"),
    ("Survival", "Subtype order of outcome", "see survival_report.json", "LumA best, basal worst", "Sørlie 2001", "10.1073/pnas.191367098", "direction"),
]
t = pd.DataFrame(rows, columns=["topic", "quantity", "ours", "published", "source", "doi", "comparability"])
t.to_csv(RES3 / "literature_table.tsv", sep="\t", index=False)
with open(RES3 / "literature_table.md", "w") as f:
    f.write("| Topic | Quantity | Ours | Published | Source | Comparability |\n|---|---|---|---|---|---|\n")
    for r in t.itertuples():
        f.write(f"| {r.topic} | {r.quantity} | {r.ours} | {r.published} | {r.source} ([doi](https://doi.org/{r.doi.split(';')[0].strip()})) | {r.comparability} |\n")
print(t[["quantity", "ours", "published", "source"]].to_string())
