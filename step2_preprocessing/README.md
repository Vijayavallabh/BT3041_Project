# BT3041 term project, step 2: TCGA-BRCA data preparation and exploratory analysis

This folder holds the code, figures, results and logs for step 2 of the team plan
(`../project_plan.md`): download, clean, normalise and filter the data, then explore it with PCA, UMAP,
violin plots, heatmaps and clustering, and check whether the clusters match the clinical receptor groups.
Steps 3 (DE, enrichment, classifiers, validation) and 4 (slides, report) start from `data/processed/`.

**Data freeze v1:** see `CHANGELOG.md`. After 23 Oct, every change to `data/processed/` gets a new CHANGELOG entry.

## Contents

| Path | What it is |
|---|---|
| `scripts/` | Every step, numbered in run order (table below) |
| `src/common.py` | Paths, seed, the fixed colour palette, plotting style, `drop_excluded()`, `write_params()` |
| `data/processed/` | Analysis-ready matrices and tables (see the data dictionary below and `DATA_DICTIONARY.tsv`) (not committed; `run_all.sh` creates it) |
| `data/SHA256SUMS.freeze` | Checksums of every raw download at the data freeze |
| `data/raw_small/` | Copies of the small raw files (clinical tables, references, HGNC set) so the freeze can be rebuilt exactly |
| `figures/` | Every figure as PNG (300 dpi) and PDF. `fig*` = main (report-sized; `fig1_sample_flow` is 7.5 in wide with 9 pt text), `figS*` = supplementary, `slide*` = 16:9 slide composites with text of 11 pt or more (`slide3_sample_flow_notes.txt` holds slide 3's footnotes for the speaker notes) |
| `results/<step>/` | Tables behind every figure, `params.json` (every parameter plus package versions; `results/cluster/params_figures.json` for the clustering figures) and a `*_report.json` of every number quoted |
| `logs/` | Console output of the last full run |
| `requirements.txt`, `environment.yml` | Pinned Python 3.12 environment |

## How to reproduce

```bash
python3.12 -m venv .venv && . .venv/bin/activate && pip install -r requirements.txt
#   or: conda env create -f environment.yml && conda activate bt3041
bash run_all.sh          # about 25 min on 4 cores, plus the downloads (METABRIC expression is 689 MB)
```

Windows: run the commands from Git Bash and activate with `source .venv/Scripts/activate`, or skip activation and
call `PY=.venv/Scripts/python.exe bash run_all.sh`. `run_all.sh` passes the same interpreter to every step.

`run_all.sh` first copies the small frozen raw files from `data/raw_small/` into `data/raw/` (never overwriting),
downloads anything still missing, and checks every raw file against `data/SHA256SUMS.freeze`. It stops on any
mismatch, including the HGNC set, which is re-released monthly at the same URL (`scripts/verify_sha256.py
--allow-volatile` accepts a newer HGNC copy with a warning, but the freeze is then no longer exact). It then runs every
script and writes its log to `logs/`. The large raw files are not in the shared folder; they are re-downloaded
(HiSeqV2, STAR counts and METABRIC expression, about 850 MB in all). All randomness uses seed 42.

| Script | Does | Main outputs | Time |
|---|---|---|---|
| `00_download.sh` | Downloads all sources, verifies checksums | `data/raw/` | download |
| `01_preprocess.py` | Labels, QC, exclusions, gene harmonisation, matrices | `data/processed/*`, `results/preprocess/` | ~3 min |
| `02_explore.py` | PCA, η² batch check, UMAP, marker tests, PAM50 heatmap, RNA–protein | fig2–6, figS1–S5 | ~2 min |
| `03_cluster.py` | k-means scan, consensus clustering, hierarchical, agreement with labels, sensitivity | `results/cluster/`, `tcga_clusters.tsv`, centroids | ~15 min |
| `03b_cluster_figures.py` | Clustering figures, ARI ceilings, purity/site check, Netanely comparison | fig7–9, figS7–S8 | <1 min |
| `04_metabric_explore.py` | METABRIC PCA and batch check, TCGA centroids applied to METABRIC | fig10, figS9 | ~1 min |
| `05_sample_flow.py` | Sample-flow diagram (report and 16:9 slide versions), library-size and identity QC | fig1, slide3, figS0 | seconds |
| `06_slide_figures.py` | Two 16:9 slide composites with large type | slide4, slide5 | seconds |
| `07_data_dictionary.py` | Shapes and checksums of the processed files | `DATA_DICTIONARY.tsv`, `SHA256SUMS` | <1 min |

## Data sources

| Source | File | Used for |
|---|---|---|
| UCSC Xena TCGA hub | `TCGA.BRCA.sampleMap/HiSeqV2` (log2 RSEM normalised + 1), `BRCA_clinicalMatrix`, `RPPA_RBN` | expression, IHC/FISH labels, protein |
| UCSC Xena GDC hub | `TCGA-BRCA.star_counts` (log2 count + 1), GENCODE v36 probemap | integer counts for DE |
| Xena survival / PanCanAtlas | `BRCA_survival.txt` (TCGA-CDR, Liu 2018), `TCGASubtype.20170308` | survival, PAM50 comparator |
| cBioPortal datahub (Git LFS) | `brca_metabric`: Illumina HT-12 expression, patient and sample tables | external validation |
| HGNC | `hgnc_complete_set.txt` (official Google Cloud bucket) | gene-symbol harmonisation |
| MSigDB (mirror) | Hallmark v2024.1 GMT | enrichment (step 3) |
| Parker 2009 / UNC | PAM50 centroids | PAM50 gene list and comparator |

## Data dictionary (`data/processed/`)

Sample IDs are 15-character TCGA barcodes (`TCGA-XX-XXXX-01` tumour, `-11` adjacent normal); a patient is the first
12 characters. Every matrix is **samples × genes**, with the index named `sample`. Matrices are **unscaled**:
classifiers must fit any filter, scaling or feature selection inside each training fold (plan section 2.2).
Exact shapes, the row key of each file and checksums are in `DATA_DICTIONARY.tsv` and `SHA256SUMS`. There, `n_cols`
counts every column of a TSV (its key column included) but only the data columns of a parquet matrix (genes or
antibodies; the `sample` index is not counted), as in the table below.

| File | Rows × columns | Contents | Used by |
|---|---|---|---|
| `tcga_labels.tsv` | 1,097 × 45 | One row per primary tumour: ER, PR, HER2 (v1 rule) with `HER2_source` and `HER2_confidence`, `HER2_v0`, `clinical_group`, `clinical_group_v0`, `ERneg_HER2neg`, histology, stage, PAM50, SigClust 2012, OS/PFI/DSS, QC flags, `exclude`, `exclude_reason`, `supervised_set`, `star_available`, `rppa_available` | everyone |
| `tcga_hiseq_tumour_log2.parquet` | 1,088 × 20,201 | log2(RSEM + 1), HGNC symbols, excluded tumours removed, only constant genes dropped | M2, M4 |
| `tcga_hiseq_normal_log2.parquet` | 113 × 20,201 | adjacent normals (the normal of the redacted case dropped), same genes | M2, M3 |
| `tcga_star_counts_tumour_filtered.parquet` | 1,087 × 24,660 | integer counts, versioned Ensembl IDs, genes with ≥10 counts in ≥ n tumours (n = smallest clinical group) | M3 (PyDESeq2) |
| `tcga_star_counts_all.parquet` | 1,207 × 60,660 | all STAR tumours after vial dedup (excluded tumours included: filter with `tcga_labels.exclude`) and 112 normals (the redacted case's normal removed), no gene filter | M3 |
| `tcga_rppa.parquet` | 742 × 131 | RPPA (RBN) for primary tumours; filter with `tcga_labels.exclude` | M5 |
| `rppa_antibody_map.tsv` | 131 antibodies | gene, `total_protein`, `phospho_or_cleaved`, `isoform_ambiguous` | M5 |
| `metabric_log2.parquet` | 1,980 × 20,193 | METABRIC log2 intensity, HGNC symbols | M2, M3, M4 |
| `metabric_labels.tsv` | 1,980 rows | `ER_IHC`, `HER2_SNP6`, `HER2_gain`, `contrast_group`, survival, PAM50 + claudin-low (comparison only) | M2, M3, M4 |
| `tcga_clusters.tsv` | 1,088 rows | cluster IDs for every method and k (C0 = largest cluster) | M2, M5 |
| `tcga_kmeans_k{2,4,5}_centroids_z.tsv` | 2,000 genes × k | centroids in z-units; assign another cohort by Euclidean nearest centroid after z-scoring within it | M2, M4 |
| `metabric_clusters.tsv` | 1,980 rows | nearest TCGA centroid per k (Pearson columns are a labelled sensitivity only) | M2 |
| `gene_map.tsv` | one row per input gene per platform | HGNC approved symbol and how it was matched | everyone |
| `gencode_v36_annotation.tsv` | 60,660 Ensembl IDs | GENCODE gene name, chromosome, start, end, HGNC symbol | M3 |
| `erbb2_amplicon_genes.tsv` | 215 genes | genes within 1 / 1.5 / 2 Mb of ERBB2 (for the amplicon ablation) | M4 |
| `pam50_centroids_hgnc.tsv` | 50 genes | PAM50 centroids with HGNC symbols | M2, M4 |
| `h.all.v2024.1.hgnc.gmt` | 50 sets | Hallmark sets mapped to HGNC symbols | M3 |
| `sex_specific_genes.txt` | list | chrY genes, XIST, TSIX (left out of unsupervised gene sets) | M2 |

Names in the plan map to these files as follows: `labels.tsv` → `tcga_labels.tsv`, `tumour_log2.parquet` →
`tcga_hiseq_tumour_log2.parquet`, `counts_int.parquet` → `tcga_star_counts_tumour_filtered.parquet`,
`rppa.parquet` → `tcga_rppa.parquet`, `clusters.tsv` → `tcga_clusters.tsv`.

## Final sample counts

All numbers below come from `results/preprocess/preprocess_report.json` and are drawn in `figures/slide3_sample_flow` (all four data sets) and `figures/fig1_sample_flow` (report version: TCGA RNA as boxes, RPPA and METABRIC as text).

| Set | n | HR+/HER2− | HER2+ | TNBC | Labels incomplete |
|---|---|---|---|---|---|
| HiSeqV2 primary tumours downloaded | 1,097 | | | | |
| Analysed tumours (9 excluded), all unsupervised work | 1,088 | 596 | 166 | 157 | 169 |
| Supervised set (female), HiSeqV2: classifiers | 1,076 | 591 | 160 | 157 | 168 |
| Supervised set with STAR counts: DE | 1,075 | 591 | 160 | 157 | 167 |
| Supervised set with RPPA: protein validation | 730 | 395 | 114 | 109 | 112 |
| TCGA contrast matched to METABRIC (STAR, PR ignored) | | ER+/HER2− 582 | | ER−/HER2− 170 | |
| METABRIC with expression | 1,980 | ER+/HER2 not gained 1,220 | HER2 gained 426 | ER−/HER2 not gained 286 | 48 |

Under the plan's HER2 rule (`clinical_group_v0`), the supervised set would be 595 / 176 / 159. In the supervised set,
ER is positive in 789, negative in 235 and missing in 52; HER2 is positive in 161, negative in 752 and missing in 163
(`ER_counts_final_supervised`, `HER2_counts_final_supervised`). Of the 160 HER2+ tumours, 72 were called by FISH,
58 by IHC score and 30 by status text alone (`HER2_confidence == "low"`).

Survival events in the 1,088 analysed tumours: OS 151, PFI 144, DSS 83 (METABRIC: OS 1,143, breast-cancer deaths 646).
A Cox model on the supervised set with a complete clinical group and known age and stage has 895 tumours with
100 OS and 100 PFI events (`survival_events_cox_set`), not the 173 PFI events the plan assumed, so keep a PFI
model to about 7–10 parameters.


## Decisions and deviations from the plan (for the team to confirm before the freeze)

Each item changes a number in the plan or adds a rule the plan did not have. The plan's preliminary counts were
computed before QC, so every count downstream should come from `results/preprocess/preprocess_report.json`
(keys ending in `_final`), not from the plan.

1. **HER2 rule (v1, used).** (1) FISH positive or negative decides when present, except that IHC 3+ with a negative
   FISH is treated as missing (the two tests disagree; 5 tumours). (2) Without a conclusive FISH, the IHC score: 3+
   positive, 0 or 1+ negative, 2+ missing (equivocal). (3) The IHC status text is used only when the score is blank.
   (4) If the status text and the score contradict each other and there is no FISH, HER2 is missing. The plan's rule
   (FISH, then the status text only) is kept as `HER2_v0` / `clinical_group_v0` for a sensitivity analysis.
   *Why:* under v0, 17 of the 1,097 tumours were HER2-positive on a status text of "Positive" without FISH while
   their recorded IHC score was 1+ (11) or 2+ (6), and 2 tumours scored 3+ without FISH were negative.
   *What it does and does not change:* v1 turns no v0 positive into a negative or the reverse; it only sets
   conflicting records to missing (24 v0 calls become missing, 2 v0-missing tumours get a call). On the 661 tumours
   that v1, v0 and TCGA 2012 all label, both rules agree with the TCGA 2012 HER2 calls at the same κ (0.88), so v1 is
   cleaner, not more accurate on the shared tumours (`kappa_HER2_common_set` in the report).
   *Confidence flag:* `HER2_confidence` is "low" for calls that rest on the status text alone (no FISH, no score).
   In the supervised set that is 239 calls: 30 of the 160 HER2+, 160 of the 591 HR+/HER2−, 46 of the 157 TNBC and
   3 with incomplete labels. ERBB2 mRNA separates HER2+ from HER2− with AUC 0.70 for these calls, against 0.97 for
   IHC-score calls and 0.81 for FISH calls (FISH is mostly done on equivocal 2+ cases) (`ERBB2_auc_by_HER2_source`).
   The noise sits in the low-confidence positives: above the 99th percentile of ERBB2 in high-confidence negatives
   (log2 14.15) lie 89 of 131 high-confidence positives but only 11 of 30 low-confidence positives, and 0 of 209
   low-confidence negatives. *For step 3:* the main sensitivity analysis drops `HER2 == "Positive"` with
   `HER2_confidence == "low"` (HER2+ 160 → 130 in the supervised groups); a stricter one drops every low-confidence
   call (TNBC 111, HR+/HER2− 431, HER2+ 130). The main TNBC vs HR+/HER2− contrast keeps all calls.
   *Effect:* HER2+ falls from 176 to 160 in the supervised set; TNBC 159 → 157.
2. **Nine tumours excluded from every analysis** (`exclude == True`, reason in `exclude_reason`):
   - 5 fail the sample-identity check between HiSeqV2 and STAR (Spearman on the 15,862 shared genes with HiSeqV2
     mean log2 > 2; STAR as log2(CPM + 1)): a tumour passes only
     if its own two profiles correlate at ≥ 0.80, above the 99th percentile of correlations between different patients
     (0.87), and are each other's best match. Three match another patient better than themselves (TCGA-AC-A3QQ-01:
     ρ 0.54 with itself), so at least one of the two files may not describe the same sample.
   - 5 cases are flagged "Redacted" in TCGA-CDR (Liu 2018), the survival resource the plan uses.
   - One tumour (TCGA-A7-A0DC-01) is in both lists. Its adjacent normal (TCGA-A7-A0DC-11) is dropped as well, which
     leaves 113 HiSeqV2 and 112 STAR normals.
   - 4 of the 5 identity failures are among the 11 patients with two primary-tumour vials in STAR (item 6), which
     suggests vial or aliquot mix-ups at the source rather than a problem with the check.
3. **Male patients.** The 12 male patients stay in the unsupervised analyses but are left out of the supervised set
   (`supervised_set == False`): male breast cancer differs in biology and 12 is too few to model. Two tumours with a
   blank sex field are female by XIST/RPS4Y1 expression and stay in. No tumour recorded as female looks male by
   expression or vice versa; 3 recorded females have low XIST and no Y-gene expression, which is common in tumours that
   lose the inactive X (`sex_recorded_rows_expression_cols` in the report).
4. **Sex-specific genes** (chrY, XIST, TSIX) are removed before the variance ranking for PCA, UMAP and clustering,
   so that sex cannot drive a component.
5. **STAR gene filter** keeps genes with ≥ 10 counts in at least as many tumours as the smallest clinical group in the
   supervised set (computed at run time, `star_min_samples` in the report; 157 now, 160 in the plan).
6. **Vial rule.** For the 11 patients with two primary-tumour vials in STAR, vial 01A is kept, as the plan says. Checked:
   in 11 of 11 patients, 01A is the vial that matches HiSeqV2 best. Even so, 4 of these 11 patients fail the identity
   check and are excluded (item 2).
7. **Half-integer STAR values.** After `round(2^x − 1)`, 0.14% of values were exact .5 values. They sit in 5 columns
   that the data hub stored as the mean of two aliquots (`star_averaged_aliquots`); they are rounded half to even.
8. **Gene symbols** are mapped to current HGNC approved symbols in every dataset (order: approved symbol, Entrez ID,
   Ensembl ID, unique previous symbol, unique alias). Duplicates keep the row with the highest mean. Coverage per
   dataset is in `results/preprocess/gene_harmonisation_coverage.tsv`; Hallmark genes matched rise from 4,191 to 4,378
   of 4,384.
9. **Choice of k** was fixed before looking at agreement with labels: the k with the lowest PAC (consensus clustering,
   100 resamples of 80%). k = 4 and k = 5 are also reported (k = 5 to compare with Netanely 2016, who used
   correlation distance and included adjacent normals; here k-means is Euclidean on tumours only).
10. **PAM50 genes** are left out of the clustering genes. Only 20 of the 50 pass the top-2,000 variance cut, so the
    sensitivity runs with them included (and with all 50 forced in) give nearly the same clusters (ARI ≥ 0.92 at k = 2, 4 and 5).
11. **Assigning METABRIC to TCGA clusters** uses the k-means rule itself: Euclidean distance to the TCGA centroids in
    z-units, after z-scoring METABRIC within its own cohort. The script checks that the rule gives back the TCGA
    k-means labels (ARI 1.0). An earlier draft used Pearson correlation, which inflated the small cluster (667 vs 307
    tumours at k = 2); it is kept only as a labelled sensitivity column. Z-scoring each cohort separately assumes a
    similar subtype mix; this holds for ER (22.7% ER-negative in both cohorts), and a quantile-mapping sensitivity
    gives the same assignment (ARI 0.96 at k = 2; `quantile_mapping_sensitivity` in `metabric_report.json`).
12. **CPTAC 2020** (the plan's fallback cohort) is descoped: METABRIC downloaded and passed every check.
13. **Label sources that are not independent.** The ER/PR "nature2012" fields equal the EMR fields wherever both exist
    (κ 1.0), and Xena `PAM50Call_RNAseq` equals PanCan `Subtype_mRNA` (κ 1.0): these are the same data, not
    confirmations. Only the HER2 comparison against TCGA 2012 is informative.


## Quality control

- **Library sizes** (STAR): median and threshold in `figS0_qc_library_identity`; no tumour falls below 25% of the median.
- **Sample identity** (HiSeqV2 vs STAR): `results/preprocess/identity_qc.tsv` and figS0; 5 failures excluded (above).
- **PCA outliers:** robust z > 5 on any of PC1–5 would be flagged; none were.
- **Batch effects:** η² of PC1–10 on tissue source site, with a permutation chance level (99th percentile of 1,000
  shuffles), is in `figS1_pc_variance_by_factor`. Site explains little of PC1–4, which follow ER/PAM50, histology and
  immune signal; PC5 is the most site-associated component (η² 0.17). `figS1b_pca_histology_site` shows the sites
  that contribute most to that η² (`site_contribution_to_eta2_on_site_PC` in the report), each with its mean line.
- **Clusters and confounders:** `results/cluster/cluster_confounder_check.json` gives Cramér's V of the clusters
  against site and against PAM50, and marker-gene proxies for adipocyte, stromal, immune and proliferation content
  per cluster (`cluster_purity_site_k*.tsv`).
- **STAR vs HiSeqV2:** the same PCA run on STAR log2-CPM gives PC scores that correlate |r| ≥ 0.88 with HiSeqV2 for PC1–5.
- **METABRIC:** processing cohort explains little of PC1–3 (η² in fig10); PC1 is the ER axis.
- **Bootstrap CIs for ARI and NMI** resample tumours with the fitted partition held fixed, so they show uncertainty
  from the choice of tumours, not from re-running the clustering; clustering stability is measured separately by PAC
  and the gene-set sensitivity runs.


## Main results of step 2

These are first-pass results for slides 3–5; step 3 builds on them. Every number is in a `*_report.json`.

- **The main axis of variation is ER.** PC1 (14.4% of variance) separates ER-positive from ER-negative tumours
  (AUC 0.95 for IHC ER; η² 0.59 for clinical group, 0.82 for PAM50). Tissue source site explains η² ≤ 0.08 of PC1–4
  (chance level 0.03) and 0.17 of PC5 (`figS1`, `figS1b`).
- **A basal-like island holds 82% of TNBC.** In the UMAP, a separate island of 189 tumours holds 129 of the 157 TNBC
  and 183 of the 189 tumours there are PAM50 Basal. It also holds 22 HR+/HER2−, 9 HER2+ and 29 tumours with incomplete
  labels (`fig3`, `slide4_exploration`). HER2+ tumours do not form their own island: they mostly overlap HR+/HER2−
  tumours, and their clearest difference is ERBB2 itself.
- **Marker genes behave as expected** in the supervised set (Mann-Whitney, BH over 16 tests, `fig4`): ESR1, PGR,
  GATA3 and FOXA1 are lower in TNBC than in HR+/HER2− (rank-biserial r −0.83 to −0.95); FOXC1 is higher (r +0.83);
  ERBB2 is higher in HER2+ (r +0.65).
- **mRNA predicts the receptor tests well but not perfectly:** AUC 0.96 for ESR1 vs IHC ER, 0.92 for PGR vs PR,
  0.85 for ERBB2 vs HER2 (95% bootstrap CIs in `results/explore/explore_report.json`). For HER2 the AUC depends on how
  the label was made: 0.97 for IHC-score calls, 0.81 for FISH calls (mostly equivocal 2+ cases) and 0.70 for calls
  from the status text alone (`HER2_confidence == "low"`).
- **Clinical groups vs PAM50 (comparator, `fig6`):** 83% of TNBC are PAM50 Basal (130/157; TCGA 2012: 65/76, 86%).
  32% of clinical HER2+ are HER2-enriched (53/166), against "~50%" in TCGA 2012. Like for like, the TCGA 2012 HER2
  labels and PAM50 calls on these tumours give 39/75 (52%), which matches the paper. The gap to 32% has three roughly
  equal parts: the larger cohort (39% on the tumours that have a 2012 call), the PAM50 call (PanCan RNA-seq vs 2012
  microarray: 45%), and the HER2 label (52% with the 2012 label). Among ER+ tumours, 6.2% are Basal or HER2-enriched
  (9.1% with Normal-like), against 8.7% (10.7%) in Nielsen 2010, a different cohort of 786 ER+ tamoxifen-treated tumours
  typed by qRT-PCR.
- **mRNA vs protein (RPPA):** median Spearman 0.39 over 85 total-protein antibodies in 738 tumours; 0.89 for ERα,
  0.77 for PR, 0.64 for HER2 (`figS5`).
- **Clustering without labels and without the PAM50 genes** is most stable at k = 2 (PAC 0.02; silhouette also
  highest): a basal-like group of 217 tumours splits off, holding 188 of 189 PAM50 Basal, 138 of 157 TNBC and 179 of
  235 ER-negative tumours. It agrees with IHC ER at ARI 0.67 (95% CI 0.61–0.73), with the clinical groups at 0.42 (0.37–0.47) and
  with PAM50 at 0.41. With only 2 clusters the best reachable ARIs are 1.0, 0.88 and 0.73; PAM50 itself, collapsed to
  2 groups (luminal vs Basal + Her2-enriched) for a class-matched comparison, reaches 0.74 against ER and 0.53 against
  the clinical groups (`fig7`–`fig9`, `slide5_clustering`). k-means, consensus, Ward and average-linkage agree at
  k = 2 (ARI 0.82–1.0), the split survives other gene sets (ARI ≥ 0.95), and each of the two clusters has a
  best-match Jaccard of at least 0.82 in every alternative run (`cluster_stability_jaccard_k2.tsv`).
- **k = 5** (for comparison with Netanely 2016) gives a clean basal cluster (97% Basal), a HER2-enriched cluster
  (55% PAM50 Her2-enriched, 55% clinical HER2+), an ILC-rich LumA cluster and two mixed LumA/LumB clusters; ARI vs PAM50 is
  0.35. The ILC-rich cluster (47% ILC, 10% Normal-like) is high on the adipocyte, stromal and immune proxies and low on
  proliferation. It matches Netanely's homogeneous LumA cluster, which holds their lobular, immune-high LumA-R2 group;
  high values on every non-tumour proxy are also consistent with lower tumour content. Netanely clustered 1,035
  tumours plus 113 adjacent normals with correlation distance, so their normal cluster has no counterpart here, and
  the cluster names on slide 5 come from a heuristic (majority PAM50 class), not their published assignment.
  k = 5 is much less stable (PAC 0.14; ARI 0.67 against a 5,000-gene rerun; 0.57 against Ward). Per cluster
  (best-match Jaccard with each alternative run, `cluster_stability_jaccard_k5.tsv`), the basal cluster is
  reproducible (0.97 against both Ward and the 5,000-gene rerun) and the ILC-rich LumA cluster fairly so (0.78 and
  0.86). The HER2-enriched cluster is the least stable (0.40 against Ward, 0.54 against the 5,000-gene rerun), so it
  should be described as a tendency, not a robust group. Clusters are associated with PAM50 far more than with site
  (Cramér's V 0.63 vs 0.20).
- **METABRIC.** On 1,919 of the 2,000 clustering genes, PC1 is again the ER axis (AUC 0.94 for IHC ER) and the
  processing cohort explains η² ≤ 0.06 of PC1–3 (0.23 of PC6) (`fig10`). The TCGA k = 2 centroids put 15.5% of
  METABRIC tumours (TCGA: 19.9%) in the basal-like cluster: 307 tumours, of which 176 are PAM50 Basal, 100 claudin-low,
  25 Her2-enriched and 6 Normal-like; 91% are ER-negative by IHC (TCGA: 87%) (`figS9`). METABRIC is z-scored within its
  own cohort, which assumes a similar subtype mix; the ER-negative share is 22.7% in both cohorts, and mapping each
  METABRIC gene onto the TCGA distribution instead gives the same result (16.1% of tumours, 90% ER-negative, ARI 0.96 vs
  the main assignment). Agreement is lower than in TCGA: ARI 0.58 against IHC ER (TCGA 0.67) and 0.20 against
  4-class PAM50 (TCGA 0.41); the cluster holds 176 of the 209 PAM50 Basal tumours (TCGA 188 of 189), and a third of
  it is claudin-low. This supports the basal-like / ER-negative split in an independent cohort, platform and label
  source, more loosely than in TCGA.

**Literature numbers** were checked against the papers: TCGA 2012 (Nature 490:61, doi 10.1038/nature11412: 76 TNBC,
81 basal-like, 65 both; "~50%" of clinical HER2+ are HER2E), Nielsen 2010 (Clin Cancer Res 16:5222,
doi 10.1158/1078-0432.CCR-10-1282: 786 ER+ tumours, HER2E 8.1%, Basal 0.6%, Normal-like 2.0%) and Netanely 2016
(Breast Cancer Res 18:74, doi 10.1186/s13058-016-0724-2). TCGA 2012 analysed a subset of these same tumours, so it
is a consistency check, not an independent cohort.
