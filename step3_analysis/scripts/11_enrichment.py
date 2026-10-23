"""Step 3.2: Hallmark enrichment, DE robustness across runs, and cross-assay agreement (plan section 3D).

- GSEA prerank (gseapy) on the 50 MSigDB Hallmarks, ranking genes by the plain Wald statistic of each DE run;
  1,000 gene-set permutations, seed 42. Reported as NES with leading-edge genes. Gene-permutation q-values are
  optimistic when genes in a set are correlated, so two checks that respect inter-gene correlation are run on HiSeqV2
  (all TNBC vs all HR+/HER2-): phenotype-permutation GSEA (1,000 label permutations), and CAMERA (Wu & Smyth 2012,
  competitive test with a variance inflation factor from the estimated inter-gene correlation of each set; design
  ~ group + histology).
- ORA: Fisher's exact test of the up and the down DE lists (main run) against each Hallmark, BH over the 50 sets.
  Universe = every gene tested in the DE run that maps to an HGNC symbol.
- DE robustness: how many main-run DEGs stay DE (same sign) in each sensitivity run; Spearman of shrunk LFCs.
- Cross-assay agreement of clinical groups with PanCan PAM50 (comparator): odds ratio (95% CI), Cohen's kappa, overlap.
"""
import sys, json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import gseapy as gp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.common3 import PROC, RES3, SEED, load_labels, read_gmt, bh, write_params3, dump  # noqa: E402

OUT = RES3 / "enrichment"
OUT.mkdir(parents=True, exist_ok=True)
DE = RES3 / "de"
HALL = read_gmt()
N_PERM = 1000
report = {}


def ranks(run):
    r = pd.read_csv(DE / f"{run}_results.tsv", sep="\t", index_col=0).dropna(subset=["wald_stat", "symbol"])
    # one value per symbol: the Ensembl ID with the largest |Wald|
    r = r.reindex(r.wald_stat.abs().sort_values(ascending=False).index)
    r = r[~r.symbol.duplicated()]
    return r.set_index("symbol").wald_stat.sort_values(ascending=False), r


def prerank(rnk, label):
    pre = gp.prerank(rnk=rnk, gene_sets=HALL, permutation_num=N_PERM, seed=SEED, threads=4, outdir=None,
                     min_size=15, max_size=500, verbose=False)
    d = pre.res2d.copy()
    d = d.rename(columns={"FDR q-val": "FDR_q", "NOM p-val": "nom_p", "FWER p-val": "FWER_p", "Lead_genes": "leading_edge"})
    d["Term"] = d.Term.str.replace("HALLMARK_", "", regex=False)
    for c in ["ES", "NES", "nom_p", "FDR_q", "FWER_p"]:
        d[c] = d[c].astype(float)
    d["run"] = label
    d["n_leading_edge"] = d.leading_edge.str.count(";") + 1
    return d.sort_values("NES")[["run", "Term", "ES", "NES", "nom_p", "FDR_q", "FWER_p", "Tag %", "n_leading_edge", "leading_edge"]]


# ------------------------------------------------ GSEA prerank on every DE run
runs = [p.name.replace("_results.tsv", "") for p in sorted(DE.glob("*_results.tsv"))]
gsea_all = []
for run in runs:
    rnk, _ = ranks(run)
    d = prerank(rnk, run)
    gsea_all.append(d)
    print(run, "top:", d.tail(3).Term.tolist(), "bottom:", d.head(3).Term.tolist(), flush=True)
G = pd.concat(gsea_all)
G.to_csv(OUT / "gsea_prerank_all_runs.tsv", sep="\t", index=False)
main = G[G.run == "main"].set_index("Term")
report["gsea_main"] = {t: {"NES": round(main.loc[t, "NES"], 2), "FDR_q": main.loc[t, "FDR_q"]}
                       for t in ["E2F_TARGETS", "G2M_CHECKPOINT", "MYC_TARGETS_V1", "MYC_TARGETS_V2",
                                 "ESTROGEN_RESPONSE_EARLY", "ESTROGEN_RESPONSE_LATE", "ANDROGEN_RESPONSE",
                                 "INTERFERON_GAMMA_RESPONSE", "MTORC1_SIGNALING"] if t in main.index}
report["gsea_main_sets_FDR_lt_0.05"] = {"positive_NES": int(((main.FDR_q < 0.05) & (main.NES > 0)).sum()),
                                        "negative_NES": int(((main.FDR_q < 0.05) & (main.NES < 0)).sum())}
# agreement of NES across runs
nes = G.pivot(index="Term", columns="run", values="NES")
report["NES_spearman_vs_main"] = {r: round(stats.spearmanr(nes["main"], nes[r]).statistic, 3) for r in nes if r != "main"}
nes.to_csv(OUT / "nes_matrix_runs.tsv", sep="\t")

# ------------------------------------------------ ORA (Fisher) on the main DE lists
_, rmain = ranks("main")
universe = set(rmain.symbol)
up, down = set(rmain.symbol[rmain.direction == "up"]), set(rmain.symbol[rmain.direction == "down"])
rows = []
for direction, lst in [("up_in_TNBC", up), ("down_in_TNBC", down)]:
    for term, gs in HALL.items():
        gs = set(gs) & universe
        a = len(lst & gs); b = len(lst - gs); c = len(gs - lst); d = len(universe) - a - b - c
        orr, p = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
        rows.append({"direction": direction, "Term": term.replace("HALLMARK_", ""), "overlap": a, "set_size_in_universe": len(gs),
                     "list_size": len(lst), "odds_ratio": orr, "p": p, "genes": ";".join(sorted(lst & gs))})
ora = pd.DataFrame(rows)
ora["q_BH"] = ora.groupby("direction").p.transform(lambda p: bh(p))
ora.sort_values(["direction", "p"]).to_csv(OUT / "ora_fisher_main.tsv", sep="\t", index=False)
report["ora_main"] = {"universe": len(universe), "up": len(up), "down": len(down),
                      "sets_q_lt_0.05": ora[ora.q_BH < 0.05].groupby("direction").Term.apply(list).to_dict()}

# ------------------------------------------------ phenotype-permutation GSEA check (HiSeqV2, all TNBC vs all HR+/HER2-)
lab = load_labels()
X = pd.read_parquet(PROC / "tcga_hiseq_tumour_log2.parquet")
rng = np.random.default_rng(SEED)
tn = lab.index[lab.group == "TNBC"]; hr = lab.index[lab.group == "HR+/HER2-"]
pick = list(tn) + list(hr)
n = (len(tn), len(hr))
E = X.loc[pick]
E = E.loc[:, E.mean() >= 1].T
cls = ["TNBC"] * len(tn) + ["HRpos"] * len(hr)
ph = gp.gsea(data=E, gene_sets=HALL, cls=cls, permutation_type="phenotype", permutation_num=N_PERM, seed=SEED,
             threads=4, outdir=None, min_size=15, max_size=500, method="signal_to_noise", verbose=False)
P = ph.res2d.copy()
P["Term"] = P.Term.str.replace("HALLMARK_", "", regex=False)
for c in ["NES", "FDR q-val", "NOM p-val"]:
    P[c] = P[c].astype(float)
P = P.rename(columns={"FDR q-val": "FDR_q", "NOM p-val": "nom_p"}).sort_values("NES")
P[["Term", "ES", "NES", "nom_p", "FDR_q", "Lead_genes"]].to_csv(OUT / "gsea_phenotype_permutation_hiseq.tsv", sep="\t", index=False)
Pm = P.set_index("Term")
both = main.join(Pm[["NES", "FDR_q"]], rsuffix="_pheno")
report["phenotype_permutation_check"] = {
    "n_per_group": n, "genes": int(E.shape[0]),
    "NES_spearman_prerank_vs_phenotype": round(stats.spearmanr(both.NES, both.NES_pheno).statistic, 3),
    "sets_FDR_lt_0.05_phenotype": int((Pm.FDR_q < 0.05).sum()), "sets_FDR_lt_0.05_prerank": int((main.FDR_q < 0.05).sum()),
    "sign_agreement_among_prerank_FDR_lt_0.05": float((np.sign(both.NES) == np.sign(both.NES_pheno))[both.FDR_q < 0.05].mean()),
    "selected": {t: {"NES": round(Pm.loc[t, "NES"], 2), "FDR_q": Pm.loc[t, "FDR_q"]}
                 for t in ["E2F_TARGETS", "G2M_CHECKPOINT", "MYC_TARGETS_V1", "ESTROGEN_RESPONSE_EARLY", "ESTROGEN_RESPONSE_LATE"]}}
print("phenotype check", report["phenotype_permutation_check"], flush=True)

# ------------------------------------------------ CAMERA (correlation-adjusted competitive test), HiSeqV2
def camera(Y, design, coef, sets, genes):
    """limma::camera with estimated inter-gene correlation (ordinary t statistics; df are large here)."""
    Q, Rm = np.linalg.qr(design)
    beta = np.linalg.solve(Rm, Q.T @ Y)
    resid = Y - design @ beta
    dfr = Y.shape[0] - design.shape[1]
    s2 = (resid ** 2).sum(0) / dfr
    XtXi = np.linalg.inv(design.T @ design)
    t = beta[coef] / np.sqrt(s2 * XtXi[coef, coef])
    U = resid / np.sqrt((resid ** 2).mean(0))          # unit mean square per gene
    G = len(t)
    gidx = {g: i for i, g in enumerate(genes)}
    out = []
    for term, gs in sets.items():
        idx = np.array([gidx[g] for g in gs if g in gidx])
        m = len(idx)
        if m < 15:
            continue
        vif = m * np.mean(U[:, idx].mean(1) ** 2)
        corr = (vif - 1) / (m - 1)
        mean_in = t[idx].mean()
        delta = G / (G - m) * (mean_in - t.mean())
        var_pooled = ((G - 1) * t.var(ddof=1) - delta ** 2 * m * (G - m) / G) / (G - 2)
        tt = delta / np.sqrt(var_pooled * (vif / m + 1 / (G - m)))
        out.append({"Term": term.replace("HALLMARK_", ""), "n_genes": m, "inter_gene_corr": corr, "VIF": vif,
                    "mean_t_in_set": mean_in, "direction": "up_in_TNBC" if delta > 0 else "down_in_TNBC",
                    "camera_t": tt, "p": 2 * stats.t.sf(abs(tt), G - 2)})
    d = pd.DataFrame(out)
    d["q_BH"] = bh(d.p)
    return d.sort_values("p")


cm = lab.loc[pick]
Dm = np.column_stack([np.ones(len(cm)), (cm.group == "TNBC").astype(float), (cm.histology == "ILC").astype(float),
                      (cm.histology == "Other").astype(float)])
cam = camera(E.T.to_numpy(dtype=float), Dm, 1, HALL, list(E.index))
cam.to_csv(OUT / "camera_hiseq.tsv", sep="\t", index=False)
cj = cam.set_index("Term").join(main[["NES", "FDR_q"]])
report["camera_check"] = {"sets_q_lt_0.05": int((cam.q_BH < 0.05).sum()),
                          "prerank_FDR05_sets_also_camera_q05_same_sign": f"{int(((cj.FDR_q < 0.05) & (cj.q_BH < 0.05) & (np.sign(cj.NES) == np.sign(cj.mean_t_in_set))).sum())}/{int((cj.FDR_q < 0.05).sum())}",
                          "median_inter_gene_corr": float(cam.inter_gene_corr.median()),
                          "selected": cam.set_index("Term").loc[[t for t in ["E2F_TARGETS", "G2M_CHECKPOINT", "MYC_TARGETS_V1",
                                                                             "ESTROGEN_RESPONSE_EARLY", "ESTROGEN_RESPONSE_LATE",
                                                                             "INTERFERON_GAMMA_RESPONSE"]],
                                                                ["inter_gene_corr", "camera_t", "q_BH"]].round(4).to_dict("index")}
print("camera", report["camera_check"], flush=True)

# ------------------------------------------------ DE robustness across runs
summ = {r: json.load(open(DE / f"{r}_summary.json")) for r in runs}
pd.DataFrame(summ).T.to_csv(OUT / "de_run_summary.tsv", sep="\t")
m = pd.read_csv(DE / "main_results.tsv", sep="\t", index_col=0)
rob = {}
for r in runs:
    if r == "main":
        continue
    o = pd.read_csv(DE / f"{r}_results.tsv", sep="\t", index_col=0)
    common = m.index.intersection(o.index)
    de_m = m.loc[common][m.loc[common].DE]
    same = (o.loc[de_m.index].DE & (np.sign(o.loc[de_m.index].log2FC) == np.sign(de_m.log2FC)))
    rob[r] = {"DE_in_run": int(o.DE.sum()), "main_DEGs_tested_in_run": int(len(de_m)),
              "main_DEGs_still_DE_same_sign": int(same.sum()), "fraction": round(same.mean(), 3),
              "spearman_shrunk_LFC_all_genes": round(stats.spearmanr(m.loc[common, "log2FC_shrunk"], o.loc[common, "log2FC_shrunk"]).statistic, 3),
              "sign_agreement_main_DEGs": round((np.sign(o.loc[de_m.index].log2FC) == np.sign(de_m.log2FC)).mean(), 3)}
report["de_robustness_vs_main"] = rob
report["de_summary"] = {r: {k: summ[r][k] for k in ["contrast", "formula", "n_test", "n_ref", "genes_tested", "DE_lfc1_padj05",
                                                     "DE_up", "DE_down", "wald_padj05_no_threshold"]} for r in runs}

# positive controls and named genes from the plan (main run)
ms = m.dropna(subset=["symbol"]).drop_duplicates("symbol").set_index("symbol")
named = ["ESR1", "PGR", "FOXA1", "GATA3", "XBP1", "AGR3", "MLPH", "SLC39A6", "BCL2", "FOXC1", "KRT5", "KRT6B", "KRT17", "MKI67",
         "BIRC5", "CCNB1", "UBE2C", "EN1", "SOX10", "VGLL1", "ART3", "FABP7", "TFF1", "TFF3", "CA12", "AR", "ERBB2"]
report["named_genes_main"] = {g: {"log2FC_shrunk": round(ms.loc[g, "log2FC_shrunk"], 2), "padj_lfc1": ms.loc[g, "padj_lfc1"],
                                  "DE": bool(ms.loc[g, "DE"])} for g in named if g in ms.index}
# top DEGs other than the receptor genes, by shrunk LFC and by significance
de = ms[ms.DE & ~ms.index.isin(["ESR1", "PGR", "ERBB2"])]
de_coding = de[~de.index.str.match(r"^(A[CLP]\d{6}|LINC|RNU|MIR|SNOR|.*-AS\d|.*-DT)")]
report["top_up_TNBC_by_padj"] = de_coding[de_coding.log2FC > 0].sort_values("padj_lfc1").head(15).index.tolist()
report["top_down_TNBC_by_padj"] = de_coding[de_coding.log2FC < 0].sort_values("padj_lfc1").head(15).index.tolist()

# ------------------------------------------------ cross-assay agreement with PanCan PAM50 (comparator)
def two_by_two(a, b, name):
    a, b = a.astype(bool), b.astype(bool)
    t = pd.crosstab(a, b).reindex(index=[True, False], columns=[True, False], fill_value=0)
    (tp, fn), (fp, tn_) = t.values
    orr = (tp + .5) * (tn_ + .5) / ((fn + .5) * (fp + .5))
    se = np.sqrt(1 / (tp + .5) + 1 / (tn_ + .5) + 1 / (fn + .5) + 1 / (fp + .5))
    from sklearn.metrics import cohen_kappa_score
    return {"comparison": name, "both": int(tp), "only_clinical": int(fn), "only_PAM50": int(fp), "neither": int(tn_),
            "pct_clinical_called_PAM50": round(100 * tp / (tp + fn), 1), "pct_PAM50_with_clinical": round(100 * tp / (tp + fp), 1),
            "odds_ratio": round(orr, 1), "OR_CI": [round(np.exp(np.log(orr) - 1.96 * se), 1), round(np.exp(np.log(orr) + 1.96 * se), 1)],
            "kappa": round(cohen_kappa_score(a, b), 3), "fisher_p": stats.fisher_exact(t.values)[1]}


cg = lab[lab.group.isin(["HR+/HER2-", "HER2+", "TNBC"]) & lab.PAM50_pancan.notna()]
agree = [two_by_two(cg.group.eq("TNBC"), cg.PAM50_pancan.eq("Basal"), "TNBC vs PAM50 Basal"),
         two_by_two(cg.group.eq("HER2+"), cg.PAM50_pancan.eq("Her2"), "clinical HER2+ vs PAM50 HER2-enriched")]
er = lab[lab.ER.isin(["Positive", "Negative"]) & lab.PAM50_pancan.notna()]
erp = er[er.ER == "Positive"]
agree.append({"comparison": "ER+ tumours called non-luminal (Basal or Her2)", "n_ER_pos": len(erp),
              "non_luminal": int(erp.PAM50_pancan.isin(["Basal", "Her2"]).sum()),
              "pct": round(100 * erp.PAM50_pancan.isin(["Basal", "Her2"]).mean(), 1),
              "pct_incl_Normal_like": round(100 * erp.PAM50_pancan.isin(["Basal", "Her2", "Normal"]).mean(), 1)})
report["cross_assay_agreement"] = agree
pd.DataFrame(agree).to_csv(OUT / "cross_assay_agreement.tsv", sep="\t", index=False)

dump(report, OUT / "enrichment_report.json")
write_params3(OUT, {"gene_sets": "MSigDB Hallmark v2024.1 (HGNC-mapped, step 2)", "ranking": "Wald statistic (lfc_null=0); "
                    "duplicate symbols -> largest |stat|", "permutations": N_PERM, "min_size": 15, "max_size": 500,
                    "ora": "Fisher one-sided (greater), BH over 50 sets per direction, universe = tested genes with symbol",
                    "phenotype_check": "gseapy.gsea, signal_to_noise, phenotype permutation, all TNBC vs all HR+/HER2- HiSeqV2, "
                                       "genes with mean log2 >= 1", "camera": "Wu & Smyth 2012, estimated inter-gene correlation, design ~ group + histology"})
print(json.dumps({k: report[k] for k in ["gsea_main", "ora_main", "de_robustness_vs_main"]}, indent=1, default=str))
