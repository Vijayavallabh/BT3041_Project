"""Step 3.6: survival, exploratory (plan sections 3F.4 and 4.3).

TCGA (TCGA-CDR endpoints, Liu 2018: PFI primary, OS secondary)
  - Kaplan-Meier + log-rank by clinical group and by the label-free k = 2 cluster (basal-like vs rest).
  - Cox PFI ~ group + age (per 10 y) + stage (I / II / III-IV): 5 parameters for 100 events (20 events per
    parameter), sized to the step 2 Cox set rather than the 173 events the plan assumed. Schoenfeld test
    (lifelines proportional_hazard_test); a covariate that fails is moved to strata and the model refitted.
METABRIC (long follow-up)
  - Endpoints: OS and breast-cancer death as a cause-specific hazard (other-cause deaths censored), follow-up
    truncated at 10 years (hazards of ER- and ER+ disease cross later); a separate first-5-years model.
  - KM + log-rank by IHC/SNP6 group, by the frozen TCGA classifier's predicted ER-/HER2- group, and by the TCGA k = 2
    centroid cluster.
  - Cox (breast-cancer death, 10 y) ~ group + age + grade + tumour size, with the Schoenfeld check; if PH fails,
    the effect is split at 5 years (time-varying coefficient via episode splitting).
  - Direction check against Nielsen 2010 (ER+, tamoxifen-treated, 5-year disease-specific): ER-IHC+ tumours that got
    hormone therapy, breast-cancer death within 5 years, HRs vs LumA by METABRIC's published PAM50 call (comparator).
"""
import sys, warnings
from pathlib import Path
import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import multivariate_logrank_test, proportional_hazard_test
from lifelines.utils import to_episodic_format

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.common3 import PROC, RES3, load_labels, write_params3, dump  # noqa: E402

warnings.filterwarnings("ignore")
OUT = RES3 / "survival"
OUT.mkdir(parents=True, exist_ok=True)
report = {}


def km_table(df, t, e, g, label, horizons):
    rows = []
    for lvl, d in df.groupby(g):
        k = KaplanMeierFitter().fit(d[t], d[e])
        r = {"analysis": label, "group": lvl, "n": len(d), "events": int(d[e].sum())}
        for h in horizons:
            r[f"S({h})"] = float(k.survival_function_at_times(h).iloc[0])
        rows.append(r)
    lr = multivariate_logrank_test(df[t], df[g], df[e])
    return rows, {"chi2": lr.test_statistic, "df": int(df[g].nunique() - 1), "p": lr.p_value}


def cox_with_ph(df, duration, event, strata=None, label=""):
    cph = CoxPHFitter()
    cph.fit(df, duration, event, strata=strata)
    ph = proportional_hazard_test(cph, df, time_transform="rank")
    phs = ph.summary[["test_statistic", "p"]]
    out = {"label": label, "n": len(df), "events": int(df[event].sum()),
           "params": int(len(cph.params_)), "events_per_parameter": round(df[event].sum() / len(cph.params_), 1),
           "concordance": cph.concordance_index_,
           "HR": cph.summary[["exp(coef)", "exp(coef) lower 95%", "exp(coef) upper 95%", "p"]].round(4).to_dict("index"),
           "schoenfeld": phs.round(4).to_dict("index"), "strata": strata}
    return cph, out


# ======================================================== TCGA
lab = load_labels()
cl = pd.read_csv(PROC / "tcga_clusters.tsv", sep="\t", index_col=0)
lab["k2"] = cl.kmeans_k2.reindex(lab.index).map({0: "luminal-like (C0)", 1: "basal-like (C1)"})
tc_rows, tests = [], {}
for ep in ["PFI", "OS"]:
    d = lab[lab.group.isin(["HR+/HER2-", "HER2+", "TNBC"])].dropna(subset=[ep, f"{ep}.time"])
    d = d.assign(years=d[f"{ep}.time"] / 365.25)
    r, t = km_table(d, "years", ep, "group", f"TCGA {ep} by clinical group", [3, 5])
    tc_rows += r; tests[f"TCGA_{ep}_group"] = t
    d2 = lab.dropna(subset=[ep, f"{ep}.time", "k2"]).assign(years=lambda x: x[f"{ep}.time"] / 365.25)
    r, t = km_table(d2, "years", ep, "k2", f"TCGA {ep} by k=2 cluster", [3, 5])
    tc_rows += r; tests[f"TCGA_{ep}_k2"] = t
    # reverse KM median follow-up
    k = KaplanMeierFitter().fit(d.years, 1 - d[ep])
    tests[f"TCGA_{ep}_median_followup_years"] = float(k.median_survival_time_)
report["tcga_logrank"] = tests

cox_rows = []
cs = lab[lab.group.isin(["HR+/HER2-", "HER2+", "TNBC"]) & lab.age.notna() & lab.stage.notna()].copy()
X = pd.DataFrame({"HER2pos": (cs.group == "HER2+").astype(int), "TNBC": (cs.group == "TNBC").astype(int),
                  "age_per10y": (cs.age - cs.age.mean()) / 10, "stage_II": (cs.stage == "II").astype(int),
                  "stage_III_IV": (cs.stage == "III-IV").astype(int)}, index=cs.index)
for ep in ["PFI", "OS"]:
    d = X.join(cs[[ep, f"{ep}.time"]]).dropna()
    d = d[d[f"{ep}.time"] > 0]
    d["years"] = d.pop(f"{ep}.time") / 365.25
    cph, out = cox_with_ph(d, "years", ep, label=f"TCGA Cox {ep}")
    failing = [c for c, v in out["schoenfeld"].items() if v["p"] < 0.05]
    if failing:
        out["ph_fails"] = failing
        if any(c.startswith("stage") for c in failing):
            d2 = d.drop(columns=["stage_II", "stage_III_IV"]).assign(stage=cs.loc[d.index, "stage"])
            _, out2 = cox_with_ph(d2, "years", ep, strata=["stage"], label=f"TCGA Cox {ep} stratified by stage")
            cox_rows.append(out2)
    cox_rows.append(out)
    print(out["label"], {k: v["exp(coef)"] for k, v in out["HR"].items()}, "PH:", {k: v["p"] for k, v in out["schoenfeld"].items()}, flush=True)
report["tcga_cox"] = cox_rows

# ======================================================== METABRIC
ml = pd.read_csv(PROC / "metabric_labels.tsv", sep="\t", index_col=0)
mc = pd.read_csv(PROC / "metabric_clusters.tsv", sep="\t", index_col=0)
pr = pd.read_csv(RES3 / "transfer" / "metabric_predictions.tsv", sep="\t", index_col=0)
m = ml.join(mc.tcga_centroid_k2.map({0: "luminal-like (C0)", 1: "basal-like (C1)"}).rename("k2")).join(pr[["pred_ERneg_HER2neg", "pred_PAM50_from_TCGA"]])
m["group"] = m.contrast_group.map({"ER+/HER2-notgained": "ER+/HER2 not gained", "HER2-gained": "HER2 gained",
                                   "ER-/HER2-notgained": "ER-/HER2 not gained"})
m["pred_group"] = m.pred_ERneg_HER2neg.map({1: "predicted ER-/HER2-", 0: "predicted other"})
m["died_other"] = (m.VITAL_STATUS == "Died of Other Causes").astype(int)
m = m[m.OS_MONTHS.notna() & m.BC_death.notna()]


def truncate(df, years):
    t = df.OS_MONTHS / 12
    out = df.assign(years=np.minimum(t, years))
    out["OS_t"] = ((df.OS_event == 1) & (t <= years)).astype(int)
    out["BCD_t"] = ((df.BC_death == 1) & (t <= years)).astype(int)
    return out


mb_rows, mtests = [], {}
for yrs in [10, 5]:
    d = truncate(m, yrs)
    for g in ["group", "pred_group", "k2"]:
        dd = d.dropna(subset=[g])
        for ep in ["BCD_t", "OS_t"]:
            r, t = km_table(dd, "years", ep, g, f"METABRIC {ep[:-2]} {yrs}y by {g}", [5, 10] if yrs == 10 else [5])
            mb_rows += r; mtests[f"MB_{ep[:-2]}_{yrs}y_{g}"] = t
report["metabric_logrank"] = mtests
km_all = pd.DataFrame(tc_rows + mb_rows)
km_all.to_csv(OUT / "km_summary.tsv", sep="\t", index=False)

mcox = []
d = truncate(m, 10).dropna(subset=["group", "AGE_AT_DIAGNOSIS", "GRADE", "TUMOR_SIZE"])
Xm = pd.DataFrame({"HER2_gained": (d.group == "HER2 gained").astype(int), "ERneg_HER2notgained": (d.group == "ER-/HER2 not gained").astype(int),
                   "age_per10y": (d.AGE_AT_DIAGNOSIS - d.AGE_AT_DIAGNOSIS.mean()) / 10, "grade3": (d.GRADE == 3).astype(int),
                   "size_per10mm": (d.TUMOR_SIZE - d.TUMOR_SIZE.mean()) / 10, "years": d.years, "BCD": d.BCD_t}, index=d.index)
Xm = Xm[Xm.years > 0]
cph, out = cox_with_ph(Xm, "years", "BCD", label="METABRIC Cox breast-cancer death, 10 y")
mcox.append(out)
failing = [c for c, v in out["schoenfeld"].items() if v["p"] < 0.05]
out["ph_fails"] = failing
if failing:
    # split follow-up at 5 years: separate group effects for 0-5 y and 5-10 y
    ep = to_episodic_format(Xm.reset_index(drop=True).assign(id=range(len(Xm))), duration_col="years", event_col="BCD", id_col="id", time_gaps=5.0)
    late = (ep.start >= 5).astype(int)
    for c in ["HER2_gained", "ERneg_HER2notgained"]:
        ep[f"{c}_x_after5y"] = ep[c] * late
    ctv = CoxPHFitter().fit(ep.drop(columns=["id"]), duration_col="stop", event_col="BCD", entry_col="start")
    s = ctv.summary
    early = {c: float(np.exp(s.loc[c, "coef"])) for c in ["HER2_gained", "ERneg_HER2notgained"]}
    latehr = {c: float(np.exp(s.loc[c, "coef"] + s.loc[f"{c}_x_after5y", "coef"])) for c in ["HER2_gained", "ERneg_HER2notgained"]}
    mcox.append({"label": "METABRIC Cox BCD, group effect split at 5 y", "HR_0_5y": early, "HR_5_10y": latehr,
                 "interaction_p": {c: float(s.loc[f"{c}_x_after5y", "p"]) for c in ["HER2_gained", "ERneg_HER2notgained"]}})
d5 = truncate(m, 5).dropna(subset=["group", "AGE_AT_DIAGNOSIS", "GRADE", "TUMOR_SIZE"])
X5 = Xm.drop(columns=["years", "BCD"]).join(d5[["years", "BCD_t"]].rename(columns={"BCD_t": "BCD"}), how="inner")
X5 = X5[X5.years > 0]
_, out5 = cox_with_ph(X5, "years", "BCD", label="METABRIC Cox breast-cancer death, first 5 y")
mcox.append(out5)
# predicted ER-/HER2- group (frozen TCGA model), same adjustment
dp = truncate(m, 10).dropna(subset=["pred_ERneg_HER2neg", "AGE_AT_DIAGNOSIS", "GRADE", "TUMOR_SIZE"])
Xp = pd.DataFrame({"pred_ERneg_HER2neg": dp.pred_ERneg_HER2neg.astype(int), "age_per10y": (dp.AGE_AT_DIAGNOSIS - dp.AGE_AT_DIAGNOSIS.mean()) / 10,
                   "grade3": (dp.GRADE == 3).astype(int), "size_per10mm": (dp.TUMOR_SIZE - dp.TUMOR_SIZE.mean()) / 10,
                   "years": dp.years, "BCD": dp.BCD_t})
Xp = Xp[Xp.years > 0]
_, outp = cox_with_ph(Xp, "years", "BCD", label="METABRIC Cox BCD 10 y, predicted ER-/HER2- (frozen TCGA model)")
mcox.append(outp)

# Nielsen 2010 direction check (comparator; published PAM50 calls of METABRIC)
dn = truncate(m[(m.ER_IHC == "Positive") & (m.HORMONE_THERAPY == "YES") & m.CLAUDIN_SUBTYPE.isin(["LumA", "LumB", "Her2", "Basal"])], 5)
Xn = pd.get_dummies(dn.CLAUDIN_SUBTYPE)[["LumB", "Her2", "Basal"]].astype(int).assign(years=dn.years, BCD=dn.BCD_t)
Xn = Xn[Xn.years > 0]
cn = CoxPHFitter().fit(Xn, "years", "BCD")
nielsen = {"n": len(Xn), "events": int(Xn.BCD.sum()), "subtype_counts": dn.CLAUDIN_SUBTYPE.value_counts().to_dict(),
           "HR_vs_LumA": cn.summary[["exp(coef)", "exp(coef) lower 95%", "exp(coef) upper 95%", "p"]].round(3).to_dict("index"),
           "published_Nielsen2010_HR_vs_LumA": {"LumB": 1.99, "Her2": 3.65, "Basal": 17.71}}
report["metabric_cox"] = mcox
report["nielsen_direction_check"] = nielsen
print("METABRIC cox:", {k: v["exp(coef)"] for k, v in mcox[0]["HR"].items()}, "PH fails:", mcox[0]["ph_fails"], flush=True)
print("Nielsen check:", {k: v["exp(coef)"] for k, v in nielsen["HR_vs_LumA"].items()}, flush=True)
dump(report, OUT / "survival_report.json")
write_params3(OUT, {"tcga_endpoints": "PFI primary, OS secondary (TCGA-CDR)", "tcga_cox": "group + age/10y + stage (I ref)",
                    "metabric_endpoints": "OS; breast-cancer death cause-specific (other deaths censored)",
                    "truncation_years": [10, 5], "metabric_cox": "group + age/10y + grade 3 + size/10mm",
                    "ph_test": "lifelines proportional_hazard_test, rank time transform; fail -> stratify (TCGA stage) or split at 5 y"})
