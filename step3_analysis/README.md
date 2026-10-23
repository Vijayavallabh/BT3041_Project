# BT3041 term project, step 3: core analysis and validation

This folder holds step 3 of the team plan (`../project_plan.md`, sections 3C–3F
and 4): differential expression with FDR, pathway enrichment, classification of the clinical receptor status, and
validation against protein (RPPA), an independent cohort (METABRIC), survival and published numbers. It reads step 2's
`data/processed/` and never writes to it. Slides 6–9 (and a survival backup slide) follow step 2's slides 3–5.

**Labels are frozen (23 Oct 2026).** Step 2's HER2 rule v1, the 9 exclusions and the groups 591/160/157 were reviewed
against ASCO/CAP 2018 HER2 guidance (Wolff et al. 2018, doi 10.1200/JCO.2018.77.8738) and kept unchanged; the reasons are
in "Decisions and deviations" item 0. Labels stay parameterised for the sensitivity runs: every script gets its labels from `src/common3.py::load_labels`, driven by `config.json`
(`her2_rule`: `v1` step 2 rule | `v0` plan rule; `her2_filter`: `none` | `drop_low_pos` | `drop_low_all`) or by the
environment variables `BT3041_HER2_RULE` / `BT3041_HER2_FILTER`. If step 2 changes `tcga_labels.tsv`, `bash run_all.sh`
reruns everything; figure titles and headlines are computed from the result files, so nothing typed by hand goes stale.

## How to reproduce

```bash
python3.12 -m venv .venv && . .venv/bin/activate && pip install -r requirements.txt
bash run_all.sh            # about 1.5 h on 4 cores (classification ~45 min); needs ../step2_preprocessing/data/processed
```

| Script | Does | Main outputs |
|---|---|---|
| `10_de.py` | PyDESeq2, 7 runs (below) | `results/de/` |
| `11_enrichment.py` | GSEA prerank, Fisher ORA, CAMERA, label-permutation GSEA, DE robustness, PAM50 agreement | `results/enrichment/` |
| `12_classify.py` | Nested CV, 4 models × ablations, baselines, corrected t-tests, controls, HER2 label sensitivity, PAM50 benchmark | `results/classify/` |
| `13_freeze_transfer.py` | Freezes the METABRIC models on TCGA (reads only METABRIC gene names) and writes `CHANGELOG.md` `metabric-freeze` | `results/transfer/frozen/` |
| `14_metabric.py` | Transfer to METABRIC, DE replication, Hallmark agreement | `results/transfer/` |
| `15_rppa.py` | RPPA tests, mRNA–protein concordance with permutation null, ER error analysis | `results/rppa/` |
| `15b_classifier_features.py` | Stable classifier features and their Hallmark ORA | `results/classify/stable_*` |
| `16_survival.py` | KM, log-rank, Cox with Schoenfeld checks (TCGA, METABRIC), Nielsen 2010 direction check | `results/survival/` |
| `17_figures_de.py`, `18_figures_validation.py` | Report figures and 16:9 slides | `figures/` |
| `19_literature_table.py` | Our numbers next to published ones, with DOIs | `results/literature_table.{tsv,md}` |

Every `results/<step>/` has a `params.json` (all parameters, package versions, label configuration) and a
`*_report.json` holding every number quoted below. All randomness uses seed 42. Logs of the last run are in `logs/`.

## Figures

| File | For | Shows |
|---|---|---|
| `slide6_de` | slide 6 | Volcano (TNBC vs HR+/HER2−) and DEG counts per sensitivity run |
| `slide7_pathways` | slide 7 | Hallmark NES with correlation-aware checks; TNBC vs PAM50 Basal 2×2; METABRIC NES agreement |
| `slide8_classification` | slide 8 | ROC for ER and HER2 (gene alone, best model, ablated) and AUC of all 4 models × gene sets |
| `slide9_validation` | slide 9 | DE vs RPPA concordance; ER error analysis vs ERα protein; METABRIC transfer |
| `slideS1_survival` | backup | KM: TCGA PFI by group; METABRIC breast-cancer death by group and by the frozen TCGA model |
| `fig11`–`fig16` | report | Report-sized versions (DE, Hallmarks, classification, RPPA, METABRIC, survival) |
| `figS10_her2_volcano`, `figS11_nes_runs` | supplement | HER2+ vs HR+/HER2− volcano (17q12 dominates); Hallmark NES across all DE runs |

Slides are exactly 13.33 × 7.5 in (16:9) with text ≥ 10.5 pt, in step 2's colours (HR+/HER2− blue, HER2+ orange, TNBC green).

## Main results

Counts are for the female supervised set under the current labels (HR+/HER2− 591, HER2+ 160, TNBC 157).

**1. Differential expression (slide 6, fig11).** PyDESeq2 on STAR counts, `~ histology + group`, testing
|log2FC| > 1 directly (`lfc_null = 1`, `greaterAbs`). TNBC (157) vs HR+/HER2− (591): **2,741 DEGs** at FDR 5%
(1,484 higher, 1,257 lower in TNBC) of 23,439 genes; the plain Wald test calls 18,604 (79%), which is why we test the
threshold itself. Positive controls lead (ESR1 log2FC −5.1, PGR −4.4); beyond them the ER regulon is lower in TNBC (GATA3 −3.3,
AGR3 −4.8, FOXA1 −2.9, XBP1 −2.4) and basal transcription factors are higher (FOXC1 +3.1, EN1, SOX10; KRT5/KRT17 +2.4).
MKI67 and BIRC5 are higher (log2FC ≈ 1.1–1.2) but not significantly beyond 2-fold.
Robustness (share of main DEGs still DE with the same sign): + age 94%, IDC only 91% (so lobular biology is not being
reported as "HR+ biology"), strict labels (every low-confidence HER2 call dropped, 111 vs 431) 87%, and the
METABRIC-matched contrast ER−/HER2− vs ER+/HER2− 96%. HER2+ vs HR+/HER2−: 148 DEGs, 33 within 2 Mb of ERBB2
(GRB7, STARD3, PGAP3, MIEN1, PNMT, TCAP, PPP1R1B, MED1, CDK12…); dropping low-confidence HER2+ calls gives 189 DEGs.

**2. Pathways (slide 7, fig12).** GSEA prerank on the Wald statistic: 28 Hallmarks at FDR < 0.05. Up in TNBC:
E2F targets (NES +3.11), G2M (+2.88), allograft rejection (+2.74), MYC V1 (+2.74), MYC V2 (+2.38), IFN-γ (+2.18),
mTORC1 (+2.16). Down: estrogen response early (−3.08) and late (−2.40). Fisher ORA on the DEG lists agrees
(estrogen early/late down; G2M, TNFα/NF-κB up at q < 0.05). **Caveat:** two tests that
respect co-expression, CAMERA (Wu & Smyth 2012; median inter-gene correlation 0.06, up to 0.25 for E2F) and
phenotype-permutation GSEA, agree with the prerank direction for all 28 sets (NES ρ 0.92) but only **estrogen
response early** stays at q < 0.05 (CAMERA). So the proliferation/MYC/immune programmes are shifted consistently and
with large effect sizes, but "more shifted than the average gene" is not established once co-expression is accounted for; report NES
and leading edges, not q ≈ 0. Clinical TNBC vs PanCan PAM50 Basal: κ 0.78, 83% of TNBC are Basal (TCGA 2012: 86%),
81% of Basal are TNBC (80%); clinical HER2+ called HER2-enriched 33% (κ 0.40).

**3. Classification (slide 8, fig13).** Nested CV (outer 5-fold × 3, inner 5-fold), all filtering/scaling/selection
inside the training folds. **No model beats the receptor gene alone.** ER: ESR1 AUC 0.959 vs best model (Random
Forest) 0.954 (Nadeau-Bengio corrected t p = 0.69). HER2: ERBB2 0.845 vs Random Forest 0.851 (p = 0.66). With
ESR1/PGR/ERBB2, every gene within 1.5 Mb of ERBB2 and the 50 PAM50 genes removed, **ER is unchanged** (0.948–0.953;
paired tests vs all genes p ≥ 0.18), so the ER signal comes from a regulon, not one gene (stable features AGR3, DACH1).
**HER2 drops** by 0.02–0.06 AUC (0.784–0.820; significant for L1-logistic and RF), identical at 2 Mb, so ~0.80 AUC of
HER2 signal lies outside the amplicon (largely the ER-low/proliferative character of HER2+ tumours; call it
"copy-number-independent signal" only with that caveat). Shuffled labels give MCC −0.01 (ER) and +0.04 (HER2).
HER2 label sensitivity: dropping low-confidence HER2+ calls raises ERBB2 alone to 0.894 and models to 0.88; dropping
all low-confidence calls gives 0.886. The noisy status-text calls cost ~0.05 AUC, as step 2 predicted.
Secondary benchmark (circular, agreement with a centroid rule): 5-class PAM50 87.9% accuracy (87.2% without PAM50
genes), 4-class 91.0%.

**4. Protein (slide 9, fig14).** 92 of 131 antibodies differ between TNBC (109) and HR+/HER2− (395) at q < 0.05.
DE log2FC vs RPPA median difference over 85 genes with a total-protein antibody: **Spearman ρ 0.80, permutation p =
0.001** (null 95th percentile |ρ| 0.49); all 15 DEGs with an antibody agree in sign. HER2+ vs HR+/HER2−: ρ 0.58
(p = 0.001). **Error analysis:** 15 IHC-ER− tumours are called ER+ out of fold; 8 of the 10 with RPPA have ERα above
the 95th percentile of correctly called ER− tumours (median 1.87 vs −1.10) and most are PAM50 luminal, so they are
likely IHC false negatives. The 68 IHC-ER+ tumours called ER− are mostly PAM50 Basal/Her2 with low ERα (26/43 below the 5th
percentile of correct ER+): ER-low tumours or IHC false positives.

**5. METABRIC (slide 9, fig15).** Models frozen on TCGA before any METABRIC label was read (`CHANGELOG.md`).
Transfer AUC, METABRIC vs TCGA CV: ER (IHC) **0.959** (95% CI 0.947–0.970) vs 0.950; ER ablated 0.952; HER2 (SNP6
gain) **0.906** (0.886–0.925) vs 0.838; HER2 ablated 0.763; ER−/HER2− vs rest 0.935 vs 0.964. No AUC drop for ER and
a gain for HER2, unlike the 13–22% drop reported by Shah & Patel 2026 for TCGA→CPTAC: METABRIC's copy-number HER2
label tracks mRNA more closely than TCGA's IHC/FISH label. Quantile mapping instead of z-scoring changes AUC by
< 0.012. The frozen models still do not beat the gene alone in METABRIC (ESR1 0.960, ERBB2 0.900).
DE replication (OLS with age, 286 vs 1,220): **86% of the 1,857 TCGA DEGs on the array replicate** with the same
sign at q < 0.05 (sign agreement 95%; effect ρ 0.79); Hallmark NES ρ 0.94, 26/26 TCGA-significant sets same sign.

**6. Survival (backup slide, fig16).** TCGA (median follow-up ~2.4 y): PFI differs by group (log-rank p = 0.007);
Cox PFI ~ group + age + stage (884 tumours, 100 events, 5 parameters): **TNBC HR 2.69 (1.70–4.26)**, HER2+ 0.98.
METABRIC, breast-cancer death within 10 years (cause-specific, 1,846 tumours, 492 deaths; age, grade, size): HER2
gained HR 1.90, ER−/HER2 not gained 1.60, but PH fails for ER−/HER2−: split at 5 years the HR is **2.46 in years 0–5
and 0.51 in years 5–10**, the known early-relapse pattern of ER-negative disease. The frozen TCGA classifier's
predicted ER−/HER2− group also separates METABRIC survival (log-rank p = 3e-8). Nielsen 2010 direction check (ER+
endocrine-treated, 5-year breast-cancer death, METABRIC's published PAM50): HR vs LumA LumB 2.57, HER2E 3.37, Basal
8.17 (published 1.99, 3.65, 17.71): same order.

**7. Literature table:** `results/literature_table.md` (21 rows, each with a DOI checked against PubMed on 21 Oct).

## Decisions and deviations from the plan

0. **Label sign-off (23 Oct 2026): current labels confirmed, no change.**
   - *HER2 rule v1 follows ASCO/CAP 2018:* IHC 3+ positive, 0/1+ negative, 2+ equivocal and decided by ISH (FISH);
     the lab's final status is used only when no score or FISH is recorded. The one departure is that IHC 3+ with a
     negative FISH (5 tumours) is set to missing rather than positive: ASCO/CAP asks for such discordant cases to be
     re-reviewed, which a public dataset cannot do, so leaving them out is the conservative reading.
   - *Low-confidence HER2 calls are kept in the main analysis* and dropped only in the sensitivity runs. Dropping them
     from the main set would mean removing labels because they disagree with ERBB2 expression, which biases the HER2
     classifier upwards (the outcome we test). Reporting both (ERBB2 AUC 0.845 main, 0.894 without them) avoids
     that bias and shows its size.
   - *The 9 exclusions* are required: 5 sample-identity failures (the expression may belong to another patient) and
     redacted TCGA-CDR cases, which TCGA itself says not to analyse.
   - *Groups* use the standard clinical surrogate definitions (HR+ = ER+ or PR+; TNBC = ER−, PR−, HER2−) on female
     patients, so they match the published comparators in the literature table.
1. **Counts differ from the plan** (step 2's final labels): DE 157 vs 591, not 160 vs 605; Cox set 100 events, not 173,
   so the TCGA Cox model has 5 parameters (20 events each).
2. **Two correlation-aware enrichment checks.** The plan asked for one phenotype-permutation check; CAMERA was added
   because the check came back with no set at FDR < 0.05. The slide reports this rather than q ≈ 0.
3. **PLS-DA threshold** is the midpoint of the two class means of the training-fold scores, not 0.5 (0.5 favours the
   majority class when classes are unbalanced).
4. **Random Forest** tunes only k (number of genes), not tree parameters, to keep run time near 45 min.
5. **Stable features** use L1-logistic weights in ≥ 12/15 folds. With strong regularisation chosen by the inner CV,
   only 2–4 genes are stable per configuration, so their Hallmark ORA finds nothing; the slides do not rely on it.
6. **Transfer threshold** frozen at probability 0.5; prevalence-matched metrics are also in `transfer_metrics.tsv`.
7. **TCGA Cox:** age violates proportional hazards in the PFI model (Schoenfeld p = 0.03); kept as a covariate
   (the group effect, which is what we report, passes PH, p = 0.49).
8. **METABRIC survival:** the plan's 10-year truncation plus a 5-year model; the group effect is additionally split at
   5 years because PH fails for ER−/HER2−.

## Limitations

- TNBC vs HR+/HER2− differ in proliferation and immune infiltration, which are partly tumour-content effects; DE was
  not adjusted for purity.
- PAM50 is used only as a comparator; the TNBC/Basal agreement is a concordance, not a validation of either.
- METABRIC ER−/HER2− uses PR-agnostic groups and copy-number HER2, so it is not called TNBC.
- TCGA survival is exploratory: short follow-up, treatment not modelled.
- Classifier comparisons rest on 15 overlapping folds; the corrected t-test is conservative and the paired bootstrap
  is the second check (both in `summary.tsv`).

## AI use
The AI prompts and the AI working log for this step are submitted with the report as supplementary material (not in
this repository).
