"""Step 3.1: differential expression with PyDESeq2 on STAR counts (plan section 3C).

Main contrast: TNBC vs HR+/HER2- with design ~ histology + group (histology = IDC / ILC / Other).
A gene is DE when padj < 0.05 for H0: |log2FC| <= 1 (lfc_null = 1, alt = greaterAbs), i.e. the fold-change
threshold is tested, not applied after the test. The plain Wald test (lfc_null = 0) is also run: its statistic ranks
genes for GSEA, and its DE count shows why the threshold test is needed. Shrunk LFCs (apeGLM-type prior, as
DESeq2's lfcShrink) are used for the volcano plot and the rankings.

Runs (each is one PyDESeq2 fit, written to results/de/<run>_results.tsv):
  main                TNBC vs HRpos, ~histology + group, all label calls
  idc_only            same, IDC tumours only, ~group (is lobular biology being reported as "HR+ biology"?)
  age                 ~age_z + histology + group
  strict_conf         drop every HER2 call that rests on the status text alone (step 2 stricter sensitivity)
  matched_erher2      ER-/HER2- vs ER+/HER2- (PR ignored), ~histology + group; matched to the METABRIC contrast
  her2                HER2+ vs HR+/HER2-, ~histology + group (secondary)
  her2_drop_low_pos   same, low-confidence HER2+ calls dropped (step 2 main sensitivity)
Usage: python scripts/10_de.py [run ...]     (no argument = all runs; ~1-2 h on 4 cores)
"""
import sys, time
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.common3 import PROC, RES3, SEED, load_labels, gene_symbols_star, write_params3, dump  # noqa: E402
from pydeseq2.dds import DeseqDataSet  # noqa: E402
from pydeseq2.ds import DeseqStats  # noqa: E402
from pydeseq2.default_inference import DefaultInference  # noqa: E402

OUT = RES3 / "de"
OUT.mkdir(parents=True, exist_ok=True)
N_CPUS = 4
ALPHA, LFC_NULL = 0.05, 1.0
sym, ann = gene_symbols_star()
counts_all = pd.read_parquet(PROC / "tcga_star_counts_tumour_filtered.parquet")
assert counts_all.index.is_unique


def design_for(run):
    """(metadata DataFrame with `group` levels A (test) / B (reference), design formula, description)."""
    if run in ("main", "idc_only", "age"):
        lab = load_labels()
        m = lab[lab.group.isin(["TNBC", "HR+/HER2-"])].copy()
        m["group"] = m.group.map({"TNBC": "TNBC", "HR+/HER2-": "HRpos"})
        test, ref = "TNBC", "HRpos"
    elif run == "strict_conf":
        lab = load_labels(her2_filter="drop_low_all")
        m = lab[lab.group.isin(["TNBC", "HR+/HER2-"])].copy()
        m["group"] = m.group.map({"TNBC": "TNBC", "HR+/HER2-": "HRpos"})
        test, ref = "TNBC", "HRpos"
    elif run == "matched_erher2":
        lab = load_labels()
        m = lab[lab.erher2.notna()].copy()
        m["group"] = m.erher2.map({"ERneg_HER2neg": "ERneg", "ERpos_HER2neg": "ERpos"})
        test, ref = "ERneg", "ERpos"
    elif run in ("her2", "her2_drop_low_pos"):
        lab = load_labels(her2_filter="drop_low_pos" if run == "her2_drop_low_pos" else None)
        m = lab[lab.group.isin(["HER2+", "HR+/HER2-"])].copy()
        m["group"] = m.group.map({"HER2+": "HER2pos", "HR+/HER2-": "HRpos"})
        test, ref = "HER2pos", "HRpos"
    else:
        raise ValueError(run)
    m = m[m.index.isin(counts_all.index)]
    if run == "idc_only":
        m = m[m.histology == "IDC"]
        formula = "~group"
    elif run == "age":
        m = m[m.age.notna()]
        m["age_z"] = (m.age - m.age.mean()) / m.age.std()
        formula = "~age_z + histology + group"
    else:
        formula = "~histology + group"
    # reference levels first alphabetically is not guaranteed for every contrast, so order categories explicitly
    m["group"] = pd.Categorical(m.group, categories=[ref, test])
    m["histology"] = pd.Categorical(m.histology, categories=["IDC", "ILC", "Other"])
    keep = ["group", "histology"] + (["age_z"] if "age_z" in m else [])
    return m[keep], formula, test, ref


def run_one(run):
    t0 = time.time()
    md, formula, test, ref = design_for(run)
    if run == "idc_only":
        md = md.drop(columns="histology")
    X = counts_all.loc[md.index]
    # keep the step 2 gene filter (>= 10 counts in >= smallest-group n tumours) but re-apply it to this subset
    nmin = int(md.group.value_counts().min())
    X = X.loc[:, (X >= 10).sum(0) >= nmin]
    inf = DefaultInference(n_cpus=N_CPUS)
    dds = DeseqDataSet(counts=X, metadata=md, design=formula, refit_cooks=True, inference=inf, quiet=True)
    dds.deseq2()
    coef = [c for c in dds.obsm["design_matrix"].columns if c.startswith("group")]
    assert len(coef) == 1, coef
    # threshold test (the DE call)
    st1 = DeseqStats(dds, contrast=["group", test, ref], alpha=ALPHA, lfc_null=LFC_NULL, alt_hypothesis="greaterAbs",
                     inference=inf, quiet=True)
    st1.summary()
    r1 = st1.results_df
    # plain Wald test (ranking statistic, and the count that shows why the threshold matters)
    st0 = DeseqStats(dds, contrast=["group", test, ref], alpha=ALPHA, inference=inf, quiet=True)
    st0.summary()
    r0 = st0.results_df.copy()
    st0.lfc_shrink(coeff=coef[0])
    rs = st0.results_df
    res = pd.DataFrame({
        "symbol": sym.reindex(X.columns).values, "chrom": ann.chrom.reindex(X.columns).values,
        "baseMean": r0.baseMean, "log2FC": r0.log2FoldChange, "lfcSE": r0.lfcSE,
        "log2FC_shrunk": rs.log2FoldChange, "wald_stat": r0.stat, "pvalue_wald": r0.pvalue, "padj_wald": r0.padj,
        "stat_lfc1": r1.stat, "pvalue_lfc1": r1.pvalue, "padj_lfc1": r1.padj}, index=X.columns)
    res.index.name = "ensembl_id"
    res["DE"] = (res.padj_lfc1 < ALPHA)
    res["direction"] = np.where(res.DE & (res.log2FC > 0), "up", np.where(res.DE & (res.log2FC < 0), "down", ""))
    res.sort_values("pvalue_lfc1").to_csv(OUT / f"{run}_results.tsv", sep="\t")
    summ = {
        "run": run, "formula": formula, "contrast": f"{test} vs {ref}", "n_test": int((md.group == test).sum()),
        "n_ref": int((md.group == ref).sum()), "genes_tested": int(X.shape[1]), "gene_filter_min_samples": nmin,
        "DE_lfc1_padj05": int(res.DE.sum()), "DE_up": int((res.direction == "up").sum()),
        "DE_down": int((res.direction == "down").sum()),
        "wald_padj05_no_threshold": int((res.padj_wald < ALPHA).sum()),
        "wald_padj05_and_absLFC_gt1_post_hoc": int(((res.padj_wald < ALPHA) & (res.log2FC.abs() > 1)).sum()),
        "histology_counts": md.get("histology", pd.Series(dtype=str)).astype(str).value_counts().to_dict()
        if "histology" in md else {"IDC": len(md)},
        "cooks_refit_genes": int(dds.varm["refitted"].sum()) if "refitted" in dds.varm else None,
        "seconds": round(time.time() - t0)}
    dump(summ, OUT / f"{run}_summary.json")
    print(run, {k: summ[k] for k in ["n_test", "n_ref", "genes_tested", "DE_lfc1_padj05", "DE_up", "DE_down",
                                     "wald_padj05_no_threshold", "seconds"]}, flush=True)


RUNS = ["main", "idc_only", "age", "strict_conf", "matched_erher2", "her2", "her2_drop_low_pos"]
if __name__ == "__main__":
    np.random.seed(SEED)
    runs = sys.argv[1:] or RUNS
    write_params3(OUT, {"alpha": ALPHA, "lfc_null": LFC_NULL, "alt_hypothesis": "greaterAbs", "refit_cooks": True,
                        "shrinkage": "pydeseq2 lfc_shrink (apeGLM-like Cauchy prior) on the group coefficient",
                        "gene_filter": ">=10 counts in >= n tumours of the smaller group, recomputed per run",
                        "runs": RUNS, "n_cpus": N_CPUS}, name="params.json")
    for r in runs:
        run_one(r)
