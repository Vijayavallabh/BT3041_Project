"""Step 1: build the frozen, analysis-ready data set (data freeze v1).

Inputs : data/raw/** (from scripts/00_download.sh)
Outputs: data/processed/*  (data dictionary in README.md) and results/preprocess/*

Every parameter that changes the data is in PARAMS and is written to
results/preprocess/params.json, so the report can quote it verbatim. All group counts
in preprocess_report.json ending in "_final" are after every exclusion and are the
numbers to quote; keys ending in "_supervised" are the female supervised set (DE groups and
classifier targets). Cross-tabulations are stored row-first: the key name says which variable is the row.
Run from the project root:  python scripts/01_preprocess.py
"""
import json, sys
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import cohen_kappa_score

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.common import RAW, PROC, RES, write_params

IDENTITY_RHO_MIN = 0.80      # absolute floor for same-tumour HiSeqV2-vs-STAR Spearman
IDENTITY_NULL_PCTL = 99      # self-rho must also beat this percentile of between-patient rho
IDENTITY_MIN_MEAN_LOG2 = 2   # identity QC uses genes on both platforms with HiSeqV2 mean log2 above this
STAR_MIN_COUNT = 10  # minimum-sample threshold is set from the data below (smallest DE group)
SEX_EXPR_CUT = 6     # log2 expression cut for XIST / RPS4Y1 in the expression-based sex call
LOW_LIB_FRACTION = 0.25

PARAMS = {
    "tumour_sample_code": "01 = primary solid tumour (used); 11 = adjacent normal (kept separately); 06 = metastatic (dropped)",
    "vial_dedup_rule": "if a patient has >1 primary-tumour vial in STAR counts, keep vial 01A (checked: 01A is the better match to HiSeqV2 in 11/11; 4 of these 11 patients still fail the identity QC and are excluded)",
    "ER_PR_source": "EMR fields breast_carcinoma_{estrogen,progesterone}_receptor_status; Indeterminate/blank = missing",
    "HER2_rule": ("v1 (used): (1) ISH/FISH Positive -> pos, FISH Negative -> neg (FISH has precedence). "
                  "(2) Without a conclusive FISH, IHC score her2_immunohistochemistry_level_result: 3+ -> pos, 0/1+ -> neg, 2+ -> missing (equivocal). "
                  "(3) If the score is blank, fall back to the IHC status text (Positive/Negative). "
                  "(4) If status text and score contradict (status Positive with score 0/1+, or Negative with 3+) and there is no FISH -> missing. "
                  "(5) IHC 3+ with a negative FISH is discordant (IHC 3+ alone is positive under ASCO/CAP) -> missing. "
                  "HER2_confidence = low for calls that rest on the IHC status text alone (no score, no FISH). "
                  "v0 (team plan, kept as HER2_v0 for sensitivity): FISH overrides IHC status text; IHC equivocal without FISH = missing."),
    "clinical_group_rule": "needs ER, PR and HER2 all known; HER2+ if HER2 pos; else HR+/HER2- if ER or PR pos; else TNBC",
    "ERneg_HER2neg_rule": "ER neg and HER2 neg (PR ignored, to match METABRIC)",
    "exclusions_all_analyses": (f"(a) sample identity not confirmed: same-tumour Spearman between HiSeqV2 and STAR log2-CPM < {IDENTITY_RHO_MIN}, "
                                f"or below the {IDENTITY_NULL_PCTL}th percentile of between-patient Spearman, or the tumour's own STAR profile is not its best match in both directions; "
                                "(b) case flagged 'Redacted' in TCGA-CDR (Liu 2018); the adjacent normals of redacted cases are dropped too"),
    "identity_qc_genes": (f"HGNC genes present on both platforms with HiSeqV2 mean log2 > {IDENTITY_MIN_MEAN_LOG2}; STAR as log2(CPM + 1), "
                          "Ensembl IDs collapsed per HGNC symbol by maximum; Spearman via ranks"),
    "expression_sex": (f"Male if RPS4Y1 > {SEX_EXPR_CUT} and XIST < {SEX_EXPR_CUT}; Female if XIST >= {SEX_EXPR_CUT}; otherwise unclear "
                       "(used only when the recorded sex is blank)"),
    "supervised_set": "tumours not excluded AND female (recorded sex, or XIST/RPS4Y1 expression when the record is blank); the 12 male patients are kept for unsupervised analyses only (male breast cancer differs in biology; too few to model)",
    "histology_collapse": "IDC = Infiltrating Ductal; ILC = Infiltrating Lobular; Other = everything else incl. mixed and blank",
    "stage_collapse": "Stage I/IA/IB -> I; II/IIA/IIB -> II; III*/IV -> III-IV; X, blank, [Discrepancy] -> missing",
    "pam50_main": "PanCanAtlas TCGASubtype.20170308 Subtype_mRNA (Xena PAM50Call_RNAseq is a subset of the same calls)",
    "hiseq_gene_filter_at_freeze": "drop genes with zero variance across the analysed tumours (no other filter at freeze)",
    "star_count_conversion": ("counts = round(2**x - 1) from Xena log2(count+1), np.rint (half to even). Non-integers are exact .5 values "
                              "confined to a few 01A columns that the hub stored as the mean of two aliquots (flagged star_averaged_aliquots)"),
    "star_gene_filter": f">= {STAR_MIN_COUNT} counts in >= n_min analysed primary tumours, n_min = size of the smallest clinical group "
                        "in the supervised set with STAR counts (computed at run time, logged as star_min_samples)",
    "gene_harmonisation": "HGNC complete set: approved symbol > Entrez ID (incl. IDs embedded in LOC<id> / ?|<id> names) > Ensembl ID > unique previous symbol > unique alias; otherwise keep the original name (method 'unmapped')",
    "duplicate_gene_rule": "when several rows map to one HGNC symbol keep the row with the highest mean expression",
    "sex_specific_genes": "chrY genes (GENCODE v36) plus XIST and TSIX; listed in sex_specific_genes.txt and removed from unsupervised gene sets",
    "erbb2_amplicon_windows_bp": [1_000_000, 1_500_000, 2_000_000],
    "rppa": "primary tumours only; phospho/cleaved antibodies flagged and excluded from mRNA-protein concordance (not collapsed); pan-isoform antibodies flagged isoform_ambiguous",
    "metabric_er_label": "ER_IHC ('Positve' typo fixed); blank = missing",
    "metabric_her2_label": "HER2_SNP6 GAIN = gained; NEUTRAL/LOSS = not gained; UNDEF = missing",
    "metabric_missing_values": "genes with >5% missing dropped; remaining missing values imputed by the gene median",
    "low_library_rule": f"flag (not drop) STAR samples with total counts below {LOW_LIB_FRACTION:.0%} of the median",
}

OUT = RES / "preprocess"
for d in (PROC, OUT):
    d.mkdir(parents=True, exist_ok=True)
report = {}


def log(k, v):
    report[k] = v
    print(f"{k}: {v}", flush=True)


def to_parquet(df, path):
    df = df.copy()
    df.index.name = "sample"
    df.columns = df.columns.astype(str)
    df.to_parquet(path)


# ---------------------------------------------------------------- HGNC map
hg = pd.read_csv(RAW / "ref/hgnc_complete_set.txt", sep="\t", low_memory=False,
                 usecols=["symbol", "prev_symbol", "alias_symbol", "entrez_id", "ensembl_gene_id"])
APPROVED = set(hg.symbol)
prev, alias = defaultdict(set), defaultdict(set)
for s, p, a in zip(hg.symbol, hg.prev_symbol, hg.alias_symbol):
    for x in str(p).split("|") if isinstance(p, str) else []:
        prev[x].add(s)
    for x in str(a).split("|") if isinstance(a, str) else []:
        alias[x].add(s)
ENTREZ = {int(e): s for e, s in zip(hg.entrez_id, hg.symbol) if pd.notna(e)}
ENSG = {e: s for e, s in zip(hg.ensembl_gene_id, hg.symbol) if isinstance(e, str)}


def to_hgnc(symbol=None, entrez=None, ensg=None):
    if isinstance(symbol, str) and symbol in APPROVED:
        return symbol, "approved"
    if entrez is not None and pd.notna(entrez) and int(entrez) in ENTREZ:
        return ENTREZ[int(entrez)], "entrez"
    if isinstance(ensg, str) and ensg.split(".")[0] in ENSG:
        return ENSG[ensg.split(".")[0]], "ensembl"
    if isinstance(symbol, str):
        if len(prev.get(symbol, ())) == 1:
            return next(iter(prev[symbol])), "previous"
        if len(alias.get(symbol, ())) == 1:
            return next(iter(alias[symbol])), "alias"
        return symbol, "unmapped"
    return None, "unmapped"


gene_map_rows = []


def harmonise(df, dataset, symbols=None, entrez=None, ensg=None):
    """Return df re-indexed by HGNC symbol (duplicates collapsed by highest mean) and log the map."""
    n = len(df)
    symbols = symbols if symbols is not None else [None] * n
    entrez = entrez if entrez is not None else [None] * n
    ensg = ensg if ensg is not None else [None] * n
    mapped = [to_hgnc(s, e, g) for s, e, g in zip(symbols, entrez, ensg)]
    new = pd.Series([m[0] for m in mapped], index=df.index)
    how = pd.Series([m[1] for m in mapped], index=df.index)
    for i, s, e, g in zip(df.index, symbols, entrez, ensg):
        gene_map_rows.append((dataset, i, s, e, g, new[i], how[i]))
    keep = new.notna()
    d = df[keep.values].copy()
    d["_sym"], d["_mean"] = new[keep].values, d.mean(axis=1, numeric_only=True).values
    d = d.sort_values("_mean", ascending=False).drop_duplicates("_sym")
    d = d.set_index("_sym").drop(columns="_mean")
    d.index.name = "gene"
    log(f"{dataset}_harmonisation", {"rows_in": n, **how.value_counts().to_dict(), "unique_hgnc_out": len(d)})
    return d


def kappa(a, b):
    m = a.notna() & b.notna()
    return {"n": int(m.sum()), "kappa": round(float(cohen_kappa_score(a[m], b[m])), 3),
            "agreement": round(float((a[m] == b[m]).mean()), 3)}


# ---------------------------------------------------------------- TCGA clinical labels
clin = pd.read_csv(RAW / "tcga/BRCA_clinicalMatrix.tsv", sep="\t", low_memory=False).set_index("sampleID")
hiseq = pd.read_csv(RAW / "tcga/HiSeqV2.gz", sep="\t", index_col=0)
log("hiseq_raw_shape_genes_x_samples", list(hiseq.shape))
log("hiseq_sample_types", hiseq.columns.str[13:15].value_counts().to_dict())
tum_ids = sorted(c for c in hiseq.columns if c[13:15] == "01")
nrm_ids = sorted(c for c in hiseq.columns if c[13:15] == "11")
assert len(set(c[:12] for c in tum_ids)) == len(tum_ids), "one tumour per patient expected"

c = clin.reindex(tum_ids)
pos_neg = lambda s: s.where(s.isin(["Positive", "Negative"]))
ihc = c["lab_proc_her2_neu_immunohistochemistry_receptor_status"]
score = c["her2_immunohistochemistry_level_result"]
fish = c["lab_procedure_her2_neu_in_situ_hybrid_outcome_type"]

# HER2 v0 (team plan rule, kept for sensitivity)
her2_v0 = pd.Series(np.nan, index=c.index, dtype=object)
her2_v0[ihc == "Negative"] = "Negative"
her2_v0[fish == "Negative"] = "Negative"
her2_v0[(ihc == "Positive") & (fish != "Negative")] = "Positive"
her2_v0[fish == "Positive"] = "Positive"

# HER2 v1 (used): FISH first, then IHC score, then status text; contradictions -> missing
her2 = pd.Series(np.nan, index=c.index, dtype=object)
src = pd.Series(np.nan, index=c.index, dtype=object)
fpos, fneg = fish == "Positive", fish == "Negative"
nofish = ~(fpos | fneg)
conflict = nofish & (((ihc == "Positive") & score.isin(["0", "1+"])) | ((ihc == "Negative") & (score == "3+")))
ihc3_fishneg = fneg & (score == "3+")
rules = [(fpos, "Positive", "FISH"), (fneg & ~ihc3_fishneg, "Negative", "FISH"),
         (nofish & ~conflict & (score == "3+"), "Positive", "IHC score"),
         (nofish & ~conflict & score.isin(["0", "1+"]), "Negative", "IHC score"),
         (nofish & score.isna() & (ihc == "Positive"), "Positive", "IHC status (no score)"),
         (nofish & score.isna() & (ihc == "Negative"), "Negative", "IHC status (no score)")]
for m, v, sname in rules:
    her2[m], src[m] = v, sname
src[conflict] = "conflict (status vs score) -> missing"
src[ihc3_fishneg] = "IHC 3+ vs FISH negative -> missing"
src[nofish & (score == "2+")] = src[nofish & (score == "2+")].fillna("IHC 2+ without FISH -> missing")
log("HER2_label_source", src.value_counts(dropna=False).to_dict())
log("HER2_v0_rows_v1_cols", pd.crosstab(her2_v0.fillna("missing"), her2.fillna("missing")).to_dict(orient="index"))

lab = pd.DataFrame(index=pd.Index(tum_ids, name="sample"))
lab["patient"] = lab.index.str[:12]
lab["tss"] = lab.index.str[5:7]
lab["sex"] = c["gender"].str.title()
lab["age"] = pd.to_numeric(c["age_at_initial_pathologic_diagnosis"], errors="coerce")
lab["ER"] = pos_neg(c["breast_carcinoma_estrogen_receptor_status"])
lab["PR"] = pos_neg(c["breast_carcinoma_progesterone_receptor_status"])
lab["HER2"] = her2
lab["HER2_source"] = src
lab["HER2_confidence"] = np.where(her2.isna(), None, np.where(src == "IHC status (no score)", "low", "high"))
lab["HER2_v0"] = her2_v0
lab["HER2_IHC_status"], lab["HER2_IHC_score"], lab["HER2_FISH"] = ihc, score, fish


def groups(her2col):
    complete = lab[["ER", "PR"]].notna().all(axis=1) & her2col.notna()
    g = pd.Series("Incomplete", index=lab.index)
    g[complete & (her2col == "Positive")] = "HER2+"
    g[complete & (her2col == "Negative") & ((lab.ER == "Positive") | (lab.PR == "Positive"))] = "HR+/HER2-"
    g[complete & (her2col == "Negative") & (lab.ER == "Negative") & (lab.PR == "Negative")] = "TNBC"
    return g


lab["clinical_group"] = groups(lab.HER2)
lab["clinical_group_v0"] = groups(lab.HER2_v0)
eh = lab[["ER", "HER2"]].notna().all(axis=1)
lab["ERneg_HER2neg"] = np.where(eh, ((lab.ER == "Negative") & (lab.HER2 == "Negative")).astype(int), np.nan)
for k in ["ER_Status_nature2012", "PR_Status_nature2012", "HER2_Final_Status_nature2012"]:
    lab[k] = pos_neg(c[k])
h = c["histological_type"]
lab["histology"] = np.select([h == "Infiltrating Ductal Carcinoma", h == "Infiltrating Lobular Carcinoma"],
                             ["IDC", "ILC"], "Other")
st = c["pathologic_stage"].fillna("")
lab["stage"] = np.select([st.str.fullmatch(r"Stage I[AB]?"), st.str.fullmatch(r"Stage II[AB]?"),
                          st.str.match(r"Stage (III|IV)")], ["I", "II", "III-IV"], None)
sub = pd.read_csv(RAW / "tcga/TCGASubtype.20170308.tsv.gz", sep="\t", index_col=0)
sub = sub[sub.Subtype_Selected.astype(str).str.startswith("BRCA")]
lab["PAM50_pancan"] = sub["Subtype_mRNA"].reindex(lab.index)
lab["PAM50_xena"] = c["PAM50Call_RNAseq"]
lab["PAM50_nature2012"] = c["PAM50_mRNA_nature2012"].replace({"HER2-enriched": "Her2", "Luminal A": "LumA", "Luminal B": "LumB",
                                                             "Basal-like": "Basal", "Normal-like": "Normal"})
lab["SigClust_2012"] = c["SigClust_Unsupervised_mRNA_nature2012"]
surv = pd.read_csv(RAW / "tcga/BRCA_survival.tsv", sep="\t", index_col=0)
for k in ["OS", "OS.time", "PFI", "PFI.time", "DSS", "DSS.time"]:
    lab[k] = surv[k].reindex(lab.index)
lab["cdr_redacted"] = surv["Redaction"].reindex(lab.index).eq("Redacted")

log("n_primary_tumours_hiseq", len(lab))
log("n_adjacent_normals_hiseq", len(nrm_ids))
log("normals_paired_with_tumour", int(sum(n[:12] in set(lab.patient) for n in nrm_ids)))
log("kappa_HER2_v1_vs_nature2012", kappa(lab.HER2, lab.HER2_Final_Status_nature2012))
log("kappa_HER2_v0_vs_nature2012", kappa(lab.HER2_v0, lab.HER2_Final_Status_nature2012))
both = lab.HER2.notna() & lab.HER2_v0.notna() & lab.HER2_Final_Status_nature2012.notna()
log("kappa_HER2_common_set", {"note": "tumours labelled by v1, v0 and TCGA 2012 alike",
                              "v1": kappa(lab.HER2[both], lab.HER2_Final_Status_nature2012[both]),
                              "v0": kappa(lab.HER2_v0[both], lab.HER2_Final_Status_nature2012[both]),
                              "v0_v1_opposite_calls": int(((lab.HER2 == "Positive") & (lab.HER2_v0 == "Negative")).sum()
                                                          + ((lab.HER2 == "Negative") & (lab.HER2_v0 == "Positive")).sum())})
log("note_label_sources", "ER/PR nature2012 = EMR fields where both exist (kappa 1.0, same source) and Xena PAM50Call_RNAseq = PanCan Subtype_mRNA "
                          "(kappa 1.0); these are not independent checks. Only HER2 (rule vs nature2012) is informative.")
log("kappa_ER_EMR_vs_nature2012", kappa(lab.ER, lab.ER_Status_nature2012))
log("kappa_PAM50_pancan_vs_xena", kappa(lab.PAM50_pancan, lab.PAM50_xena))
log("xena_PAM50_calls_on_adjacent_normals", clin[clin.index.str[13:15] == "11"]["PAM50Call_RNAseq"].value_counts().to_dict())

# ---------------------------------------------------------------- HiSeqV2 (harmonise; matrices written after QC)
hs = hiseq[tum_ids + nrm_ids]
emb = pd.Series(hs.index).str.extract(r"^(?:LOC|\?\|)(\d+)$")[0].astype(float)
hs_h = harmonise(hs, "HiSeqV2", symbols=list(hs.index), entrez=list(emb))
tum_all = hs_h[tum_ids]

# ---------------------------------------------------------------- STAR counts
pm = pd.read_csv(RAW / "ref/gencode.v36.probemap.tsv", sep="\t", index_col=0)
star = pd.read_csv(RAW / "tcga/TCGA-BRCA.star_counts.tsv.gz", sep="\t", index_col=0)
log("star_raw_shape_genes_x_samples", list(star.shape))
log("star_sample_types", star.columns.str[13:15].value_counts().to_dict())
raw = np.power(2.0, star.values) - 1
cnt = np.rint(raw)
nonint = pd.Series((np.abs(raw - cnt) > 1e-3).sum(0), index=star.columns)
log("star_fraction_non_integer_values_rounded", round(float(nonint.sum() / raw.size), 5))
log("star_columns_with_non_integer_values", nonint[nonint > 0].to_dict())
cnt = pd.DataFrame(cnt.astype(np.int64), index=star.index, columns=star.columns)
cols = pd.Series(star.columns, index=star.columns)
keep_cols = cols.groupby(cols.str[:15]).apply(lambda s: sorted(s)[0])
dups = {k: sorted(v) for k, v in cols.groupby(cols.str[:15]) if len(v) > 1}
log("star_vial_duplicates", dups)
averaged = sorted(v[:15] for v in keep_cols.values if nonint[v] > 0)
cnt = cnt[keep_cols.values]
cnt.columns = cnt.columns.str[:15]
assert cnt.columns.is_unique
cnt = cnt[[c for c in cnt.columns if c[13:15] in ("01", "11")]]
star_tum = sorted(c for c in cnt.columns if c[13:15] == "01")
star_nrm = sorted(c for c in cnt.columns if c[13:15] == "11")
log("star_primary_tumours_after_dedup", len(star_tum))
log("star_normals", len(star_nrm))
log("star_tumours_not_in_hiseq", sorted(set(star_tum) - set(tum_ids)))
log("hiseq_tumours_not_in_star", sorted(set(tum_ids) - set(star_tum)))
ann = pm.reindex(cnt.index)[["gene", "chrom", "chromStart", "chromEnd"]]
mapped = [to_hgnc(s, None, e) for s, e in zip(ann.gene, ann.index)]
ann["hgnc"], ann["hgnc_method"] = [m[0] for m in mapped], [m[1] for m in mapped]
for i, s, g, m in zip(ann.index, ann.gene, ann.hgnc, ann.hgnc_method):
    gene_map_rows.append(("GENCODEv36_STAR", i, s, None, i, g, m))
log("GENCODEv36_STAR_harmonisation", ann.hgnc_method.value_counts().to_dict())
lab["star_available"] = lab.index.isin(star_tum)
lab["star_averaged_aliquots"] = lab.index.isin(averaged)
lib = cnt.sum(axis=0)
lab["star_library_size"] = lib.reindex(lab.index)
low_lib = lib[lib < LOW_LIB_FRACTION * lib[star_tum].median()]
log("star_low_library_flagged", sorted(low_lib.index))

# ---------------------------------------------------------------- sample identity QC (HiSeqV2 vs STAR)
star_log = np.log2(cnt[star_tum].div(cnt[star_tum].sum(0), axis=1) * 1e6 + 1)
star_sym = star_log.groupby(ann.hgnc.reindex(star_log.index)).max()
common_g = tum_all.index.intersection(star_sym.index)
common_g = common_g[tum_all.loc[common_g].mean(axis=1) > IDENTITY_MIN_MEAN_LOG2]
cs = [s for s in tum_ids if s in star_sym.columns]


def zrank(M):
    r = M.rank().values
    return (r - r.mean(0)) / r.std(0)


A, B = zrank(tum_all.loc[common_g, cs]), zrank(star_sym.loc[common_g, cs])
R = A.T @ B / len(common_g)                       # rows = HiSeqV2 tumour, cols = STAR tumour
self_rho = pd.Series(np.diag(R), index=cs)
off = R[~np.eye(len(cs), dtype=bool)]
null_cut = float(np.percentile(off, IDENTITY_NULL_PCTL))
best_row = pd.Series(np.array(cs)[R.argmax(1)], index=cs)
best_col = pd.Series(np.array(cs)[R.argmax(0)], index=cs)
mutual = (best_row == best_row.index) & (best_col == best_col.index)
ident_fail = (self_rho < IDENTITY_RHO_MIN) | (self_rho < null_cut) | ~mutual
log("identity_qc", {"genes": len(common_g), "self_rho_median": round(float(self_rho.median()), 3),
                    "between_patient_rho_median": round(float(np.median(off)), 3),
                    f"between_patient_rho_p{IDENTITY_NULL_PCTL}": round(null_cut, 3),
                    "failed": {s: {"self_rho": round(float(self_rho[s]), 3), "best_STAR_match": best_row[s],
                                   "best_match_rho": round(float(R[cs.index(s)].max()), 3)} for s in self_rho.index[ident_fail]}})
lab["hiseq_star_spearman"] = self_rho.reindex(lab.index)
lab["identity_unconfirmed"] = ident_fail.reindex(lab.index, fill_value=False).astype(bool)
pd.DataFrame({"self_rho": self_rho, "best_between_patient_rho": pd.Series(np.where(np.eye(len(cs), dtype=bool), -1, R).max(1), index=cs)}).to_csv(
    OUT / "identity_qc.tsv", sep="\t")

# vial choice check: does vial 01A or the other vial match the HiSeqV2 profile better?
vc = []
for k, vials in dups.items():
    x = np.rint(np.power(2.0, star[vials].values) - 1)
    x = pd.DataFrame(np.log2(x / x.sum(0) * 1e6 + 1), index=star.index, columns=vials)
    x = x.groupby(ann.hgnc.reindex(x.index)).max().reindex(common_g)
    for vv in vials:
        vc.append((k, vv, round(float(tum_all.loc[common_g, k].corr(x[vv], method="spearman")), 3)))
vc = pd.DataFrame(vc, columns=["sample", "vial", "spearman_vs_hiseq"])
vc.to_csv(OUT / "vial_choice_check.tsv", sep="\t", index=False)
best = vc.loc[vc.groupby("sample").spearman_vs_hiseq.idxmax(), "vial"]
log("vial_01A_matches_hiseq_best", f"{int(best.str.endswith('A').sum())}/{len(best)}")

# ---------------------------------------------------------------- sex check (XIST vs RPS4Y1)
lab["XIST"], lab["RPS4Y1"] = tum_all.loc["XIST"], tum_all.loc["RPS4Y1"]
lab["expr_sex"] = np.where((lab.RPS4Y1 > SEX_EXPR_CUT) & (lab.XIST < SEX_EXPR_CUT), "Male",
                           np.where(lab.XIST >= SEX_EXPR_CUT, "Female", "unclear"))
log("sex_recorded_rows_expression_cols", pd.crosstab(lab.sex.fillna("NA"), lab.expr_sex).to_dict(orient="index"))
sexgenes = sorted(set(ann.loc[ann.chrom == "chrY", "hgnc"].dropna()) | {"XIST", "TSIX"})
pd.Series(sexgenes).to_csv(PROC / "sex_specific_genes.txt", index=False, header=False)

# ---------------------------------------------------------------- final flags
lab["exclude"] = lab.identity_unconfirmed | lab.cdr_redacted
lab["exclude_reason"] = np.select([lab.identity_unconfirmed & lab.cdr_redacted, lab.identity_unconfirmed, lab.cdr_redacted],
                                  ["identity unconfirmed; CDR redacted", "identity unconfirmed (HiSeqV2 vs STAR)", "TCGA-CDR redacted"], "")
lab["supervised_set"] = ~lab.exclude & ((lab.sex == "Female") | (lab.sex.isna() & (lab.expr_sex == "Female")))
log("excluded_tumours", lab.loc[lab.exclude, ["exclude_reason", "clinical_group"]].to_dict(orient="index"))
ok = lab[~lab.exclude]
sv = lab[lab.supervised_set]
log("n_tumours_final", len(ok))
log("n_supervised_set_final", len(sv))
log("ER_counts_final", ok.ER.value_counts(dropna=False).to_dict())
log("PR_counts_final", ok.PR.value_counts(dropna=False).to_dict())
log("HER2_counts_final", ok.HER2.value_counts(dropna=False).to_dict())
log("clinical_group_final_all", ok.clinical_group.value_counts().to_dict())
log("clinical_group_final_supervised", sv.clinical_group.value_counts().to_dict())
log("ER_counts_final_supervised", sv.ER.value_counts(dropna=False).to_dict())
log("HER2_counts_final_supervised", sv.HER2.value_counts(dropna=False).to_dict())
log("HER2_source_of_HER2pos_supervised", sv[sv.clinical_group == "HER2+"].HER2_source.value_counts().to_dict())
log("clinical_group_v0_final_supervised", sv.clinical_group_v0.value_counts().to_dict())
log("clinical_group_final_supervised_star", sv[sv.star_available].clinical_group.value_counts().to_dict())
log("ERneg_HER2neg_vs_ERpos_HER2neg_supervised_star", {
    "ER-/HER2-": int(((sv.star_available) & (sv.ER == "Negative") & (sv.HER2 == "Negative")).sum()),
    "ER+/HER2-": int(((sv.star_available) & (sv.ER == "Positive") & (sv.HER2 == "Negative")).sum())})
log("histology_by_group_supervised_star", pd.crosstab(sv[sv.star_available].clinical_group, sv[sv.star_available].histology).to_dict(orient="index"))
log("ER_rows_HER2_cols_final", pd.crosstab(ok.ER, ok.HER2).to_dict(orient="index"))
log("sex_final", ok.sex.value_counts(dropna=False).to_dict())
log("stage_final", ok.stage.value_counts(dropna=False).to_dict())
log("PAM50_pancan_final", ok.PAM50_pancan.value_counts(dropna=False).to_dict())
log("survival_events_final", {"OS": int(ok.OS.sum()), "PFI": int(ok.PFI.sum()), "DSS": int(ok.DSS.sum())})
cox = sv[sv.clinical_group.isin(["HR+/HER2-", "HER2+", "TNBC"]) & sv.age.notna() & sv.stage.notna()]
log("survival_events_cox_set", {"definition": "supervised set, complete clinical group, age and stage known", "n": len(cox),
                                "OS": int(cox.OS.sum()), "PFI": int(cox.PFI.sum()), "DSS": int(cox.DSS.sum())})

# ---------------------------------------------------------------- write matrices (excluded tumours removed)
keep_t = list(ok.index)
tum = tum_all[keep_t]
const = tum.var(axis=1) == 0
log("hiseq_constant_genes_dropped", int(const.sum()))
tum = tum[~const]
to_parquet(tum.T.astype("float32"), PROC / "tcga_hiseq_tumour_log2.parquet")
red_cases = set(surv.index[surv.Redaction.eq("Redacted")].str[:12])
nrm_keep = [n for n in nrm_ids if n[:12] not in red_cases]
log("normals_excluded_cdr_redacted", sorted(set(nrm_ids) - set(nrm_keep)))
to_parquet(hs_h.loc[tum.index, nrm_keep].T.astype("float32"), PROC / "tcga_hiseq_normal_log2.parquet")
log("hiseq_tumour_matrix_samples_x_genes", [tum.shape[1], tum.shape[0]])
star_nrm_drop = [c for c in star_nrm if c[:12] in red_cases]
log("star_normals_excluded_cdr_redacted", star_nrm_drop)
to_parquet(cnt.drop(columns=star_nrm_drop).T.astype(np.int32), PROC / "tcga_star_counts_all.parquet")
ct = cnt[[s for s in star_tum if s in set(keep_t)]]
STAR_MIN_SAMPLES = int(sv[sv.star_available].clinical_group.value_counts().reindex(["HR+/HER2-", "HER2+", "TNBC"]).min())
log("star_min_samples", STAR_MIN_SAMPLES)
expressed = (ct >= STAR_MIN_COUNT).sum(axis=1) >= STAR_MIN_SAMPLES
log("star_genes_passing_filter", int(expressed.sum()))
log("star_filtered_matrix_samples_x_genes", [ct.shape[1], int(expressed.sum())])
to_parquet(ct[expressed].T.astype(np.int32), PROC / "tcga_star_counts_tumour_filtered.parquet")
dupsym = ann.loc[expressed[expressed].index, "hgnc"].duplicated(keep=False)
log("star_filtered_genes_sharing_a_symbol", int(dupsym.sum()))
ann.to_csv(PROC / "gencode_v36_annotation.tsv", sep="\t")

# ---------------------------------------------------------------- 17q12 amplicon and PAM50 genes
e = ann[ann.gene == "ERBB2"].iloc[0]
log("ERBB2_coordinates_GRCh38", f"{e.chrom}:{int(e.chromStart)}-{int(e.chromEnd)}")
amp = ann[(ann.chrom == e.chrom)].copy()
amp["dist_to_ERBB2_start"] = (amp.chromStart - e.chromStart).abs()
amp = amp[amp.dist_to_ERBB2_start <= max(PARAMS["erbb2_amplicon_windows_bp"])]
for w in PARAMS["erbb2_amplicon_windows_bp"]:
    amp[f"within_{w // 1000}kb"] = amp.dist_to_ERBB2_start <= w
    log(f"amplicon_genes_within_{w // 1000}kb", {"gencode": int(amp[f"within_{w // 1000}kb"].sum()),
                                                  "in_hiseq": int(amp[amp[f'within_{w // 1000}kb']].hgnc.isin(tum.index).sum())})
amp.sort_values("chromStart").to_csv(PROC / "erbb2_amplicon_genes.tsv", sep="\t")
p50 = pd.read_csv(RAW / "ref/pam50_centroids.txt", sep="\t", index_col=0)
p50_h = harmonise(p50, "PAM50_centroids", symbols=list(p50.index))
p50_h.to_csv(PROC / "pam50_centroids_hgnc.tsv", sep="\t")
log("pam50_genes_in_hiseq", int(p50_h.index.isin(tum.index).sum()))
assert p50_h.index.isin(tum.index).all()

# ---------------------------------------------------------------- Hallmark GMT
gmt = {}
for line in open(RAW / "ref/h.all.v2024.1.Hs.symbols.gmt"):
    f = line.rstrip("\n").split("\t")
    gmt[f[0]] = f[2:]
cov = []
with open(PROC / "h.all.v2024.1.hgnc.gmt", "w") as fh:
    for k, genes in gmt.items():
        hg_genes = sorted({to_hgnc(g)[0] for g in genes})
        fh.write("\t".join([k, "hgnc-harmonised"] + hg_genes) + "\n")
        cov.append((k, len(genes), int(pd.Index(hg_genes).isin(tum.index).sum()),
                    int(pd.Index(hg_genes).isin(set(ann.hgnc)).sum())))
cov = pd.DataFrame(cov, columns=["set", "n_genes", "in_hiseq", "in_star"])
all_h = {g for v in gmt.values() for g in v}
log("hallmark_unique_genes_matched_hiseq_raw_vs_harmonised",
    [int(len(all_h & set(hiseq.index))), int(len({to_hgnc(g)[0] for g in all_h} & set(tum.index))), len(all_h)])

# ---------------------------------------------------------------- RPPA
rppa = pd.read_csv(RAW / "tcga/RPPA_RBN.gz", sep="\t", index_col=0)
log("rppa_shape_antibodies_x_samples", list(rppa.shape))
log("rppa_sample_types", rppa.columns.str[13:15].value_counts().to_dict())
rppa = rppa[[s for s in rppa.columns if s[13:15] == "01"]]
AB_GENE = {"ERALPHA": "ESR1", "ERALPHAPS118": "ESR1", "PR": "PGR", "HER2": "ERBB2", "HER2PY1248": "ERBB2",
           "GATA3": "GATA3", "AR": "AR", "BCL2": "BCL2", "CYCLINB1": "CCNB1", "ECADHERIN": "CDH1",
           "EGFR": "EGFR", "EGFRPY1068": "EGFR", "EGFRPY1173": "EGFR", "HER3": "ERBB3", "HER3PY1298": "ERBB3",
           "CYCLIND1": "CCND1", "CYCLINE1": "CCNE1", "P53": "TP53", "PTEN": "PTEN", "INPP4B": "INPP4B",
           "CAVEOLIN1": "CAV1", "CLAUDIN7": "CLDN7", "PCADHERIN": "CDH3", "NCADHERIN": "CDH2", "CMYC": "MYC",
           "IGFBP2": "IGFBP2", "STATHMIN": "STMN1", "CDK1": "CDK1", "PCNA": "PCNA", "FIBRONECTIN": "FN1",
           "BETACATENIN": "CTNNB1", "YB1": "YBX1", "SYK": "SYK", "LCK": "LCK", "CKIT": "KIT", "ASNS": "ASNS",
           "PAI1": "SERPINE1", "COLLAGENVI": "COL6A1", "MIG6": "ERRFI1", "RAD51": "RAD51", "XRCC1": "XRCC1",
           "MRE11": "MRE11", "ATM": "ATM", "BAX": "BAX", "BAK": "BAK1", "BID": "BID", "BIM": "BCL2L11",
           "BCLXL": "BCL2L1", "SRC": "SRC", "MTOR": "MTOR", "AKT": "AKT1", "TUBERIN": "TSC2", "VEGFR2": "KDR",
           "NOTCH1": "NOTCH1", "NF2": "NF2", "SMAD3": "SMAD3", "SMAD4": "SMAD4", "SMAD1": "SMAD1", "YAP": "YAP1",
           "1433EPSILON": "YWHAE", "EIF4E": "EIF4E", "EEF2": "EEF2", "EEF2K": "EEF2K", "RAD50": "RAD50",
           "53BP1": "TP53BP1", "KU80": "XRCC5", "CHK1": "CHEK1", "CHK2": "CHEK2", "LKB1": "STK11",
           "PKCALPHA": "PRKCA", "PEA15": "PEA15", "PAXILLIN": "PXN", "IRS1": "IRS1", "HSP70": "HSPA1A",
           "GAB2": "GAB2", "DJ1": "PARK7", "DVL3": "DVL3", "ERK2": "MAPK1", "MEK1": "MAP2K1", "CRAF": "RAF1",
           "JNK2": "MAPK9", "P38MAPK": "MAPK14", "P27": "CDKN1B", "P70S6K": "RPS6KB1", "S6": "RPS6",
           "FOXO3A": "FOXO3", "CD31": "PECAM1", "CD49B": "ITGA2", "BECLIN": "BECN1", "CIAP": "BIRC2",
           "STAT5ALPHA": "STAT5A", "4EBP1": "EIF4EBP1", "ACC1": "ACACA", "AMPKALPHA": "PRKAA1",
           "PI3KP110ALPHA": "PIK3CA", "GSK3ALPHABETA": "GSK3B", "CASPASE7CLEAVEDD198": "CASP7"}
ISOFORMS = {"AKT": "AKT1|AKT2|AKT3", "GSK3ALPHABETA": "GSK3A|GSK3B", "AMPKALPHA": "PRKAA1|PRKAA2",
            "HSP70": "HSPA1A|HSPA1B", "COLLAGENVI": "COL6A1|COL6A2|COL6A3", "CIAP": "BIRC2|BIRC3",
            "JNK2": "MAPK9", "ERK2": "MAPK1", "S6": "RPS6"}
abm = pd.DataFrame({"antibody": rppa.index})
abm["gene"] = abm.antibody.map(AB_GENE)
abm["isoform_genes"] = abm.antibody.map(ISOFORMS).fillna(abm.gene)
abm["isoform_ambiguous"] = abm.isoform_genes.fillna("").str.contains(r"\|")
abm["phospho_or_cleaved"] = abm.antibody.str.contains(r"P[STY]\d|PT\d|CLEAVED") & ~abm.antibody.isin(
    ["PTEN", "PAI1", "P27", "P53", "PCNA", "PAXILLIN", "PCADHERIN", "PEA15", "PKCALPHA", "PI3KP110ALPHA", "P38MAPK", "P70S6K"])
abm["total_protein"] = ~abm.phospho_or_cleaved
abm["in_hiseq"] = abm.gene.isin(tum.index)
abm.to_csv(PROC / "rppa_antibody_map.tsv", sep="\t", index=False)
to_parquet(rppa.T.astype("float32"), PROC / "tcga_rppa.parquet")
lab["rppa_available"] = lab.index.isin(rppa.columns)
log("rppa_tumours_by_group_supervised", lab[lab.rppa_available & lab.supervised_set].clinical_group.value_counts().to_dict())
log("rppa_antibodies_total_with_gene", int((abm.total_protein & abm.gene.notna()).sum()))
log("rppa_antibodies_isoform_ambiguous", abm.loc[abm.isoform_ambiguous, "antibody"].tolist())

lab.to_csv(PROC / "tcga_labels.tsv", sep="\t")

# ---------------------------------------------------------------- METABRIC
mpat = pd.read_csv(RAW / "metabric/data_clinical_patient.txt", sep="\t", comment="#", index_col=0)
msam = pd.read_csv(RAW / "metabric/data_clinical_sample.txt", sep="\t", comment="#", index_col=1)
mx = pd.read_csv(RAW / "metabric/data_mrna_illumina_microarray.txt", sep="\t")
log("metabric_expression_raw_rows_x_samples", [mx.shape[0], mx.shape[1] - 2])
assert mx.shape[1] - 2 == 1980, "METABRIC expression file incomplete?"
mx.index = [f"row{i}" for i in range(len(mx))]
sym, ent = mx.pop("Hugo_Symbol"), mx.pop("Entrez_Gene_Id")
na_frac = mx.isna().mean(axis=1)
log("metabric_missing_value_fraction_overall", round(float(mx.isna().values.mean()), 5))
keep = na_frac <= 0.05
log("metabric_rows_dropped_gt5pct_missing", int((~keep).sum()))
mx, sym, ent = mx[keep], sym[keep], ent[keep]
mx = mx.apply(lambda r: r.fillna(r.median()), axis=1)
mh = harmonise(mx, "METABRIC", symbols=list(sym), entrez=list(ent))
to_parquet(mh.T.astype("float32"), PROC / "metabric_log2.parquet")
miss50 = sorted(set(p50_h.index) - set(mh.index))
log("metabric_pam50_coverage", {"present": 50 - len(miss50), "missing": miss50})
assert not miss50
log("genes_shared_tcga_hiseq_and_metabric", int(len(tum.index.intersection(mh.index))))

ml = pd.DataFrame(index=pd.Index(mh.columns, name="sample"))
p, s = mpat.reindex(ml.index), msam.reindex(ml.index)
ml["ER_IHC"] = p.ER_IHC.replace({"Positve": "Positive"}).where(lambda x: x.isin(["Positive", "Negative"]))
ml["HER2_SNP6"] = p.HER2_SNP6
ml["HER2_gain"] = p.HER2_SNP6.map({"GAIN": 1, "NEUTRAL": 0, "LOSS": 0})
both = ml.ER_IHC.notna() & ml.HER2_gain.notna()
g = pd.Series(np.nan, index=ml.index, dtype=object)
g[both & (ml.ER_IHC == "Negative") & (ml.HER2_gain == 0)] = "ER-/HER2-notgained"
g[both & (ml.ER_IHC == "Positive") & (ml.HER2_gain == 0)] = "ER+/HER2-notgained"
g[both & (ml.HER2_gain == 1)] = "HER2-gained"
ml["contrast_group"] = g
for k in ["AGE_AT_DIAGNOSIS", "CLAUDIN_SUBTYPE", "THREEGENE", "INTCLUST", "COHORT", "OS_MONTHS", "OS_STATUS",
          "VITAL_STATUS", "RFS_MONTHS", "RFS_STATUS", "HISTOLOGICAL_SUBTYPE", "CHEMOTHERAPY", "HORMONE_THERAPY"]:
    ml[k] = p[k]
for k in ["ER_STATUS", "PR_STATUS", "HER2_STATUS", "GRADE", "TUMOR_SIZE", "TUMOR_STAGE"]:
    ml[k] = s[k]
ml["OS_event"] = ml.OS_STATUS.str.startswith("1").astype(float).where(ml.OS_STATUS.notna())
ml["BC_death"] = (ml.VITAL_STATUS == "Died of Disease").astype(float).where(ml.VITAL_STATUS.notna())
ml.to_csv(PROC / "metabric_labels.tsv", sep="\t")
log("metabric_samples_with_expression", len(ml))
log("metabric_ER_IHC", ml.ER_IHC.value_counts(dropna=False).to_dict())
log("metabric_HER2_SNP6", ml.HER2_SNP6.value_counts(dropna=False).to_dict())
log("metabric_both_labels_known", int(both.sum()))
log("metabric_contrast_groups", ml.contrast_group.value_counts().to_dict())
log("metabric_claudin_subtype", ml.CLAUDIN_SUBTYPE.value_counts(dropna=False).to_dict())
log("metabric_events", {"OS": int(ml.OS_event.sum()), "BC_death": int(ml.BC_death.sum())})

# ---------------------------------------------------------------- gene map and coverage
gm = pd.DataFrame(gene_map_rows, columns=["dataset", "row_id", "original_symbol", "entrez", "ensembl", "hgnc_symbol", "method"])
gm.to_csv(PROC / "gene_map.tsv", sep="\t", index=False)
gm.groupby(["dataset", "method"]).size().unstack(fill_value=0).to_csv(OUT / "gene_harmonisation_coverage.tsv", sep="\t")
cov.to_csv(OUT / "hallmark_coverage.tsv", sep="\t", index=False)

write_params(OUT, PARAMS)
json.dump(report, open(OUT / "preprocess_report.json", "w"), indent=2, default=str)
print("done")
