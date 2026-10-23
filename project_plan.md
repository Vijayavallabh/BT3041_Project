# BT3041 term project plan (final): TCGA-BRCA, 5-person team, Python

**Summary**
- **Question:** in the clinic, breast tumours are grouped by ER, PR and HER2 tests on tissue (IHC and FISH). What does tumour gene expression show about these groups? Does the result hold at the protein level and in an independent cohort?
- **Labels:** the differential expression (DE) groups and the classifier targets come from pathology tests, not from RNA-seq. This avoids circularity. PAM50 is used only to compare with published work.
- **Course methods covered:**
  - Exploration: PCA, UMAP, box/violin plots and heatmaps.
  - Clustering: k-means, hierarchical and consensus clustering.
  - Hypothesis tests: PyDESeq2 with a fold-change threshold and BH-FDR, then GSEA and Fisher's exact tests.
  - Classification: L1 logistic regression, linear SVM, Random Forest and PLS-DA.
  - Validation: against published numbers, against protein data (RPPA, same patients) and against METABRIC.
- **Extra:** external validation in METABRIC, using receptor labels from other assays (ER by IHC, HER2 by SNP6 copy number). It covers classifier transfer, DE replication and survival with long follow-up. **Stretch:** heterogeneity within HR+ tumours, only if METABRIC is finished by Wed 4 Nov.
- **Compute:** no GPUs are needed. Everything runs on laptops with Python 3.12.

---

## 1. Title, research questions, and how circularity is handled

**Title:** *Beyond the receptor test: what tumour transcriptomes reveal about clinically defined breast cancer groups in TCGA-BRCA, confirmed at the protein level and in an independent cohort.*

**Primary research question.** Which genes and pathways separate the clinical groups in TCGA-BRCA RNA-seq: triple-negative (TNBC) vs HR+/HER2−, and ER− vs ER+? How well does expression predict the clinical ER and HER2 status, once the receptor genes and their neighbouring genes are removed?

**Sub-questions**
- **SQ1 (unsupervised).** Cluster the tumours without labels and without the 50 PAM50 genes. How much of the clinical grouping and of the PAM50 subtypes do the clusters recover (ARI and NMI, with bootstrap CIs)? Where do they disagree?
- **SQ2 (hypothesis-driven).** In TNBC vs HR+/HER2− tumours, at BH-FDR < 0.05, adjusted for histology:
  - Are the proliferation hallmarks (E2F targets, G2M checkpoint, MYC targets) enriched?
  - Are the estrogen-response hallmarks depleted?
  - Do the DE effects agree with protein differences in the same patients (RPPA)?
- **SQ3 (prediction and transfer).** Does a classifier trained on TCGA expression predict IHC ER and IHC/FISH HER2 status better than a single-gene baseline (ESR1, ERBB2)? Use a paired statistical test, not eyeballed SDs. Does it still work without the receptor genes, the 17q12 amplicon and the PAM50 genes? In METABRIC, does it predict ER-IHC and HER2 copy-number gain, and do its predicted groups separate survival?

**Sample sizes for planning.** Counts come from BRCA_clinicalMatrix and were checked in the cloud container. M1 re-checks them on day 2 and finalises the HER2 rule.

| Item | n |
|---|---|
| HiSeqV2 columns | 1,218: 1,097 primary tumours (one per patient), 114 adjacent normals (all paired with a tumour), 7 metastatic |
| GDC STAR counts columns (used for DE) | 1,226: 1,106 primary-tumour vials. 11 patients have two vials (01A plus 01B, or 01A plus 01C for TCGA-A7-A0DB), so after keeping vial 01A there are **1,095 tumours**. Also 113 normals (all paired) and 7 metastatic. All 1,095 are in HiSeqV2; 2 HiSeqV2 tumours are not in STAR. |
| ER (EMR field `breast_carcinoma_estrogen_receptor_status`) | 807 positive / 238 negative / 52 missing or indeterminate |
| ER_Status_nature2012 (cross-check only) | 601 / 179 / 317 missing |
| PR (EMR field `breast_carcinoma_progesterone_receptor_status`) | 698 / 344 |
| HER2 (rule: FISH overrides IHC; IHC equivocal with no FISH = missing) | 185 positive / 770 negative / 142 missing |
| ER × HER2, both known (n = 954) | ER+/HER2− 596, ER+/HER2+ 142, ER−/HER2+ 42, **ER−/HER2− 174** |
| Clinical group, HiSeqV2 | HR+/HER2− 605, HER2+ 184, TNBC 161, incomplete 147 |
| Clinical group, STAR set (DE input) | HR+/HER2− 605, HER2+ 184, **TNBC 160**, incomplete 146 |
| Histology in the DE groups (STAR) | HR+/HER2−: 402 IDC, 146 ILC, 57 other or mixed. TNBC: 136 IDC, 6 ILC, 18 other. |
| Histology, all 1,097 HiSeqV2 tumours | 784 IDC, 203 ILC, 110 other |
| PAM50, PanCanAtlas `Subtype_mRNA` (**main comparator**) | Covers all 1,097: LumA 566, LumB 217, Basal 192, Her2 82, Normal 40 |
| PAM50Call_RNAseq (Xena; kappa cross-check) | LumA 421, LumB 192, Basal 141, Her2 67, Normal 23, blank 253 |
| TCGA 2012 unsupervised clusters (`SigClust_Unsupervised_mRNA_nature2012`) | 522 tumours |

**How each label is used**

| Use | Label | Circular? | Handling |
|---|---|---|---|
| DE groups and classifier targets | IHC/FISH ER, PR, HER2 and the TNBC / HR+/HER2− / HER2+ groups built from them | No. Pathology measures these on tissue with a different assay. | ESR1, PGR and ERBB2 will top every list. Treat them as positive controls, interpret the genes below them, and run an ablation that removes them. |
| Comparator | PanCan `Subtype_mRNA` (main) and Xena PAM50Call_RNAseq (cross-check). Both are computed from expression with centroids. | Yes, if predicted from expression | Used only (a) to compare with the unsupervised clusters, (b) in a secondary benchmark labelled "agreement with a centroid rule", run with and without the PAM50 genes, and (c) for concordance with the IHC groups. |
| Clustering input | Genes | It would be partly tautological if the PAM50 genes were included | The main run excludes the 50 PAM50 genes. A sensitivity run includes them. |
| METABRIC labels | ER_IHC (IHC), HER2_SNP6 (copy number) | No | Use these two only. METABRIC's ER_STATUS, PR_STATUS, HER2_STATUS, THREEGENE and CLAUDIN_SUBTYPE come from the same HT-12 arrays, so they are reported only as "agreement with an expression-derived call" (section 4). |
| Independent evidence | RPPA protein (another assay), METABRIC (other patients, another platform, ER/HER2 from other assays), survival | No | Section 3F and section 4 |
| Leakage | — | — | Every scaling, filtering, selection and tuning step happens inside the cross-validation (CV) folds (section 2.2, step 3). |

**Label quality.**
- 95 of the 119 PAM50Call_RNAseq "Normal" calls are adjacent normal tissue. Only 23 are primary tumours, and 12 adjacent normals are even called LumA. Netanely 2016 saw the same problem.
- So tumour analyses use -01 samples only, and Normal-like is dropped from 4-class subtype tasks.

---

## 2. Data

### 2.1 Files

**Core files.** All were downloaded and checked from the cloud container, about 205 MB in total. On your own machines, download them again from these URLs.

| # | File | URL | Size | Used for |
|---|---|---|---|---|
| 1 | HiSeqV2, log2(RSEM norm + 1) | https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.BRCA.sampleMap%2FHiSeqV2.gz | 64 MB | Exploration, clustering, classifiers |
| 2 | BRCA_clinicalMatrix | https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.BRCA.sampleMap%2FBRCA_clinicalMatrix | 1,247 rows | IHC labels, PAM50, TCGA 2012 clusters, stage, age, histology |
| 3 | GDC STAR counts, stored as log2(count + 1) | https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-BRCA.star_counts.tsv.gz | 138.5 MB | PyDESeq2 |
| 3b | GENCODE v36 probemap (Ensembl ID, symbol, chromosome, start, end) | https://gdc-hub.s3.us-east-1.amazonaws.com/download/gencode.v36.annotation.gtf.gene.probemap | 3.2 MB | ID mapping and amplicon coordinates |
| 4 | Survival, TCGA-CDR endpoints OS/DSS/DFI/PFI (Liu 2018) | https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/survival%2FBRCA_survival.txt | 68 KB | TCGA survival (exploratory) |
| 5 | PanCanAtlas TCGASubtype (main PAM50 source) | https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/TCGASubtype.20170308.tsv.gz | 45 KB | PAM50 comparator |
| 6 | RPPA, replicate-based normalisation: 131 antibodies × 747 tumours, 0% missing | https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.BRCA.sampleMap%2FRPPA_RBN.gz | 0.74 MB | Protein validation |
| 7 | MSigDB Hallmark v2024.1 GMT (third-party GitHub mirror) | https://raw.githubusercontent.com/jokergoo/GSEAtopics/master/h.all.v2024.1.Hs.symbols.gmt | 49 KB | GSEA and ORA. Cite MSigDB and Liberzon 2015, and compare once with the official file from gsea-msigdb.org on your own machine. |
| 8 | PAM50 genes and centroids (Parker 2009) | https://raw.githubusercontent.com/perou-lab/PAM50_classifier/main/PAM50_R/pam50_centroids.txt | 5 KB | PAM50 gene list |
| 9 | HGNC complete gene set (approved, previous and alias symbols; Entrez and Ensembl IDs) | From genenames.org. **Download on your own machine** (not reachable or checked from the container). | — | One gene-harmonisation table for every dataset and the GMT |

**External cohorts**

| # | File | Status | Used for |
|---|---|---|---|
| 10 | METABRIC `data_clinical_patient.txt` and `data_clinical_sample.txt` (cBioPortal datahub `brca_metabric`) | Downloaded through the datahub's Git LFS server; sha256 checked | ER_IHC, HER2_SNP6, survival, covariates |
| 11 | METABRIC `data_mrna_illumina_microarray.txt` (1,980 samples, log2 HT-12 intensity, with Hugo_Symbol and Entrez_Gene_Id columns) | 689 MB. Reachable from the container but not fully downloaded. **Download on your own machine.** PAM50 gene coverage in this file has not been checked yet. | Section 4 |
| 12 | CPTAC 2020 (Krug), cBioPortal `brca_cptac_2020`: `data_mrna_seq_fpkm.txt` plus clinical files | Downloaded: 18 MB; 122 tumours; all 50 PAM50 genes; ER_UPDATED_CLINICAL_STATUS 81/39; ERBB2_GENE_AMPLIFIED 17; no survival | Fallback cohort for classifier transfer |

**METABRIC download (on your own machine).**
- **Primary route:** the stdlib downloader. It was tested against the same LFS server and checks sha256.
  ```
  python3 -I fetch_cbio_lfs.py brca_metabric data/raw/metabric data_clinical_patient.txt data_clinical_sample.txt data_mrna_illumina_microarray.txt
  ```
  The script is in the shared files (see the appendix) and goes into the repo on day 1.
- **Fallback:** git-lfs, using a repo-local config so your global git settings stay unchanged. This route is untested. The files on raw.githubusercontent.com are only LFS pointers.
  ```
  GIT_LFS_SKIP_SMUDGE=1 git clone --depth 1 https://github.com/cBioPortal/datahub.git
  cd datahub
  git lfs install --local
  git -c lfs.fetchexclude='' lfs pull -I public/brca_metabric
  ```
- **Licence:** METABRIC data are under ODbL. Put this attribution in the report's data-availability line.

**Deliberately left out**
- SCAN-B: file names unverified, and GEO was unreachable from the container.
- GTEx (Toil): 1.3 GB, and mixing pipelines is risky.
- The Mertins-77 proteome: only an unofficial mirror is available.
- The Krug 2020 proteome: blocked from the container.
- PanCan and GDC RPPA: alternatives only if more antibodies are needed.

**Environment**
- Python 3.12 or newer.
- Pin pandas 2.3.x, because lifelines 0.30.3 needs pandas < 3.
- Packages: scikit-learn 1.9.1, pydeseq2 0.5.4, gseapy 1.3.1, lifelines 0.30.3, statsmodels 0.15.0, scipy, umap-learn 0.5.12, seaborn 0.13.2.
- pydeseq2, gseapy and lifelines have been installed and run on real data. umap-learn is untested, so everyone installs on day 1.

### 2.2 Preprocessing (owner: M1)

1. **Sample IDs.**
   - Use 15-character barcodes. STAR uses 16 characters, so truncate them only **after** removing duplicate vials.
   - **Vial dedup rule:** for the 11 STAR patients with two tumour vials (TCGA-A7-A26J, -A26E, -A13E, -A26F, -A13D, -A26I, -A0DB, -A13G, -A0DC; TCGA-AC-A2QH, -A3OD), keep vial 01A. All 11 have one. Log this rule in the sample-flow figure.
   - Patient = the first 12 characters.
   - Use -01 for tumour analyses. Drop the 7 metastatic (-06) samples.
2. **Labels table** (one row per tumour):
   - ER and PR from the EMR fields, and HER2 by the documented rule.
   - The clinical group, plus the ER−/HER2− flag (ignores PR, so it matches METABRIC).
   - The nature2012 receptor labels, with Cohen's kappa against the EMR fields.
   - PanCan `Subtype_mRNA` and Xena PAM50Call_RNAseq, with the kappa between them.
   - Histology collapsed to IDC / ILC / other, age, and stage collapsed to I / II / III–IV (X or blank = missing).
   - PFI and OS from file 4.
3. **HiSeqV2 matrix: the freeze ships the unscaled log2 matrix only.**
   - At the freeze, drop only genes that are constant across all tumours.
   - **Classifiers:** the expression filter (for example, expressed in at least 20% of samples), z-scoring, variance filtering and top-k selection all sit inside the sklearn `Pipeline` and are refitted in each training fold. No object may be fitted on data that includes the test fold.
   - **Exploration and clustering:** these use no labels and have no test fold, so M2 may filter and z-score the whole tumour matrix in their own notebook.
4. **STAR counts.**
   - Convert with `counts = round(2**x - 1)`.
   - Keep the Ensembl ID as the key.
   - Keep genes with at least 10 counts in at least 160 samples (the smallest group).
   - Re-tabulate the group sizes on the STAR set and use those numbers, not the HiSeqV2 ones. Use 113, not 114, for normals.
5. **Gene harmonisation table.**
   - Map every symbol in HiSeqV2 (2012-era symbols), GENCODE v36 (Ensembl), METABRIC (Entrez), CPTAC and the Hallmark GMT to the current HGNC approved symbol, using previous and alias symbols.
   - The gap is real. Without the table, only 4,191 of 4,384 Hallmark genes and 17,756 of 23,109 CPTAC symbols match HiSeqV2.
   - The table also covers the PAM50 fixes: CDCA1→NUF2, KNTC2→NDC80, ORC6L→ORC6.
   - Merge duplicate probes or symbols with one documented rule: keep the probe with the highest mean expression.
   - Supplementary table: genes matched per dataset, and Hallmark coverage per gene set.
6. **Quality control.**
   - Check library sizes and PCA outliers.
   - Colour the PCA by tissue source site (barcode characters 6–7) to look for batch effects.
   - Log every exclusion.
7. **RPPA.**
   - Join on the 15-character ID.
   - Map antibodies to genes: ERALPHA→ESR1, PR→PGR, HER2→ERBB2, GATA3, AR, BCL2, CYCLINB1→CCNB1, ECADHERIN→CDH1, EGFR.
   - Flag phospho antibodies (for example ERALPHAPS118, HER2PY1248, EGFRPY1068) so they can be collapsed or excluded in concordance tests.
8. **METABRIC.**
   - Keep the 1,980 samples that have expression.
   - Labels come **only** from other assays:
     - ER_IHC from the patient table. Fix the 'Positve' spelling: 1,498 positive / 439 negative / 43 blank.
     - HER2_SNP6 GAIN (438) vs NEUTRAL/LOSS (1,537); UNDEF (5) = missing. This is a copy-number gain call, not an amplification call, so say so in the report.
   - Both labels known: 1,932.
   - Contrast groups: **ER-IHC−/HER2-not-gained (n = 286)** vs ER-IHC+/HER2-not-gained (n = 1,220). Do not call the first group "TNBC": METABRIC has no IHC PR.
   - Apply the harmonisation table using Entrez_Gene_Id.
   - Check that all 50 PAM50 genes are present before the go/no-go.
   - For PAM50 comparisons, drop claudin-low (218) and NC (6).
9. **Data freeze v1 on Fri 23 Oct.** After that, every change gets a new version and a CHANGELOG line.

---

## 3. Analyses mapped to the course methods

### A. Exploration (M2, slide 4)
- **What:**
  - PCA and UMAP of tumours plus normals, coloured by clinical group, PAM50, tissue type, histology and tissue source site.
  - Violin plots of ESR1, PGR, ERBB2, GATA3, FOXA1, FOXC1, KRT5 and MKI67 by group.
  - Heatmap of the PAM50 genes, ordered by PAM50 and annotated with the IHC group.
- **Settings:** top-variance genes, z-scored. UMAP with n_neighbors 15–30 and seed 42. Mann-Whitney tests for the markers, with BH across the 8 genes.
- **Published comparator:**
  - The main axes should be tumour vs normal and basal/TNBC vs luminal (Perou 2000, PMID 10963602; Sørlie 2001, PMID 11553815; Sørlie 2003, PMID 12829800; TCGA 2012, PMID 23000897).
  - ESR1 and ERBB2 mRNA vs clinical status: Spearman 0.62–0.77 (Gong 2007, PMID 17329190).

### B. Clustering (M2, slide 5)
- **What:**
  - k-means for k = 2–8 with 50 restarts.
  - Hierarchical clustering: Ward, plus average linkage on 1 − Pearson.
  - Consensus clustering: 80% subsampling, 100 repeats.
- **Settings:**
  - Choose k with PAC and silhouette (Wilkerson 2010, PMID 20427518; Senbabaoglu 2014, PMID 25158761).
  - ARI and NMI with 1,000-sample bootstrap 95% CIs against:
    - PanCan PAM50 (n = 1,057 tumours without Normal-like),
    - the clinical group,
    - the TCGA 2012 SigClust clusters (n = 522).
  - Contingency heatmaps.
  - Sensitivity runs: with the PAM50 genes, and with 5,000 genes.
- **Published comparator:**
  - Netanely 2016 (PMID 27386846) used the same Xena HiSeqV2 data with k-means K = 5 on the top 2,000 genes. It found clean basal, HER2E and normal clusters, while LumA split into a pure LumA cluster and a mixed LumA/LumB cluster. Netanely reported no ARI, so ours is new.
  - Huo & Tseng 2017 (PMID 28959370) report ARI of about 0.27–0.38. The values come from a search snippet, so M5 checks Table 2 before quoting them.
  - Krug 2020 (PMID 33212010) also found a mixed LumA/B cluster.
  - Expect moderate ARI. A very high ARI would be suspicious.

### C. Differential expression with BH-FDR (M3, slide 6)
- **Main contrast:** TNBC (160) vs HR+/HER2− (605) on STAR counts, with `design="~ histology + group"` (histology = IDC / ILC / other). Use plain level names such as `TNBC` and `HRpos`; avoid "+" and "/".
- **Test against the threshold, not after it.** With 765 samples almost every expressed gene passes BH (the 30 vs 30 smoke test already gave 11,560 genes). So test |log2FC| > 1 directly:
  ```python
  ds = DeseqStats(dds, contrast=["group", "TNBC", "HRpos"], lfc_null=1, alt_hypothesis="greaterAbs")
  ```
  Call a gene DE when padj < 0.05 under this test.
  - Also report the plain Wald count (lfc_null = 0) as one sentence, to show why the threshold test is needed.
  - Use `ds.lfc_shrink()` on the group coefficient for the volcano plot and the gene rankings.
- **Sensitivity runs:**
  - IDC-only (402 vs 136). This is the main check that lobular biology is not being reported as "HR+ biology".
  - Age added to the design.
  - Report how many DEGs survive each run.
- **Matched contrast for METABRIC:** ER−/HER2− vs ER+/HER2− on STAR, ignoring PR, with the same design. On HiSeqV2 the groups are 174 vs 596; M1 re-tabulates them on STAR.
- **Dropped:** the paired tumour-vs-normal PyDESeq2 run (113 patient levels, little payoff).
- **Expected direction:**
  - Lower in TNBC: ESR1 and PGR (positive controls), then FOXA1, GATA3, XBP1, AGR3, MLPH, SLC39A6, BCL2.
  - Higher in TNBC: FOXC1 (Ray 2010, PMID 20406990), KRT5/6/17, MKI67, BIRC5, CCNB1, UBE2C (TCGA 2012; Perou 2000).
  - The pipeline smoke test (30 Basal vs 30 LumA; PAM50-defined, so only a check) gave ESR1 log2FC −6.0, PGR −7.0, FOXA1 −3.9, FOXC1 +3.6, KRT5 +2.2, MKI67 +2.1.

### D. Enrichment and categorical tests (M3, slide 7)
- **GSEA prerank** (gseapy) on the 50 Hallmarks.
  - Rank on the Wald statistic from the lfc_null = 0 test; 1,000 permutations; seed 42.
  - Present NES and leading-edge genes, not "q ≈ 0": gene-permutation q-values are too optimistic when genes are correlated.
  - One check run: phenotype-permutation `gseapy.gsea` on the log2 HiSeqV2 data. A balanced subsample is fine for speed.
- **ORA with Fisher's exact test** on the up and down DE lists, BH over the 50 sets.
  - **Universe for DE lists:** all genes tested after the count filter, mapped to HGNC symbols.
  - **Universe for M4's classifier features:** the genes eligible after the in-fold variance filter. The features themselves are those selected in at least 80% of outer folds.
- **Cross-assay agreement.** Report the concordance, not just a p-value: a Fisher p-value here will be tiny and says nothing new.
  - TNBC vs PanCan Basal: odds ratio with 95% CI, Cohen's kappa, and % overlap. TCGA 2012 found 65 of 76 TNBC were basal-like and 65 of 81 basal-like were TNBC.
  - Clinical HER2+ vs HER2E: TCGA 2012 found only about 50%.
  - ER+ tumours called non-luminal: Nielsen 2010 found 10% (PMID 20837693).
- **Published comparator:**
  - E2F_TARGETS and G2M_CHECKPOINT up in basal (Krug 2020).
  - MYC targets up (Chandriani 2009, PMID 19690609).
  - ESTROGEN_RESPONSE_EARLY/LATE down (TCGA 2012 FOXA1/ER hub; Liberzon 2015, PMID 26771021; GSEA method: Subramanian 2005, PMID 16199517).
  - Pipeline check: E2F NES +3.36, G2M +3.09, MYC_V1 +2.71, ESTROGEN_EARLY −2.99, LATE −2.21.

### E. Classification (M4, slide 8)
- **Targets:** IHC ER (807/238) and IHC/FISH HER2 (185/770).
  - Plus, for METABRIC transfer only: ER−/HER2− vs the rest (174 vs 780), with L1 logistic regression.
- **Models (4):** L1 logistic regression, linear SVM, Random Forest (500 trees), and PLS-DA (`PLSRegression` with n_components tuned, threshold at 0.5). PLS-DA covers the course "PLS" item.
- **Ablations (2):**
  - (a) All genes.
  - (b) Without the receptor genes (ESR1, PGR, ERBB2), without the whole 17q12 amplicon (every gene whose start is within ±1.5 Mb of ERBB2 at chr17:39,687,914, GRCh38, from the GENCODE v36 probemap) and without the 50 PAM50 genes.
  - The ±1.5 Mb window holds about 160 GENCODE genes, about 87 of them in HiSeqV2. The ±1 Mb window already includes co-amplified genes such as MED1, CDK12, NEUROD2, PPP1R1B, TCAP, PNMT, IKZF3, ORMDL3, GSDMB, PSMD3, MED24, THRA, NR1D1, TOP2A and the PAM50 gene CDC6.
  - Check (b) at ±2 Mb (about 110 HiSeqV2 genes).
  - Call the remaining HER2 AUC "copy-number-independent signal" only if it holds after (b).
  - For ER, list the top features left after (b) and interpret them as the ER regulon (GATA3, AGR3, XBP1 and so on), not as new markers.
- **Grid:** 2 targets × 4 models × 2 ablations = 16 configurations.
- **Settings:**
  - Nested CV: outer stratified 5-fold × 3 repeats; inner 5-fold for tuning.
  - Inside each training fold: expression filter, `StandardScaler`, variance filter, then `SelectKBest(f_classif)`. k is tuned in the inner CV from {50, 200, 1,000}.
  - `class_weight="balanced"`.
  - **No DE list from M3 is used as a CV feature set.** Those lists were computed on all TCGA tumours, test folds included.
- **Metrics and tests:**
  - ROC-AUC with 95% CI (bootstrap on pooled out-of-fold predictions), plus balanced accuracy, MCC and F1.
  - Compare each model with the single-gene baseline (ESR1 or ERBB2; threshold chosen inside the folds by Youden's J) using the **Nadeau-Bengio corrected resampled t-test** on the paired fold scores. Use DeLong or a paired bootstrap on the pooled predictions as a second check.
  - Controls: shuffled labels (MCC should be about 0) and the majority class.
- **Secondary benchmark (labelled circular):** PanCan PAM50, L1 logistic regression only, all genes vs non-PAM50 genes.
  - Run both 4-class (Normal-like dropped) and 5-class (Normal = 40), so the published 5-class numbers compare like for like.
- **Leakage audit (M4 QA checklist, applied to every workstream):** "no scaler, filter, selector or model has been fitted on data that includes the test fold."
- **Published comparator:**
  - Single gene: ER 88–96%, HER2 89–93% accuracy (Gong 2007).
  - ER test MCC 0.80–0.83 (Chierici 2020, PMID 32714870; arrays, multi-omics).
  - Multimodal upper bound: ER AUC 0.99, HER2 0.98 (Shah & Patel 2026, PMID 42407417).
  - 5-class PAM50: CV 88%, test 85–88% (Cascianelli 2020, PMID 32826944); CNN 88.42% (Mostavi 2020, PMID 32241303); RF on 11 genes 88.2% (Seo 2020, PMID 33255759); SVM 0.866 accuracy / 0.861 F1 / 0.808 MCC (Choi & Chae 2023, PMID 37101124).
  - Ceiling: published PAM50 calls shift about 15% with the reference set (85.5% concordance, Cascianelli 2020), and AIMS agrees with PAM50 only 76–77% (Paquet & Hallett 2015, PMID 25479802).

### F. Validation (M5 for literature, RPPA and error analysis; M2 for survival; slide 9)
1. **Literature table:** every headline number next to its published counterpart, with PMID.
2. **RPPA (747 tumours; M1 tabulates how many fall in each group):**
   - mRNA–protein Spearman for about 10 matched genes.
   - Mann-Whitney TNBC vs HR+/HER2− for all 131 antibodies, BH-corrected.
   - **Concordance:** Spearman between the shrunk DE log2FC and the RPPA effect (difference in medians) over the matched genes, using total-protein antibodies only (phospho forms collapsed). Get the null by permuting TNBC/HR+ labels 1,000 times within the RPPA set. This replaces the binomial sign test, which wrongly assumed the 131 antibodies are independent.
3. **Error analysis (cheap, original):** for IHC-ER+ tumours that the classifier calls ER− out of fold, and the reverse, show RPPA ERα and PAM50. This separates likely IHC errors from real ER-low biology. One figure panel.
4. **TCGA survival (M2):**
   - Kaplan-Meier and log-rank by clinical group and by cluster.
   - Cox model: group + age + stage (I / II / III–IV).
   - With 173 PFI events, keep to about one parameter per 10–15 events.
   - Check proportional hazards with lifelines `check_assumptions` (Schoenfeld residuals).
   - PFI is primary and OS secondary, because Liu 2018 (PMID 29625055) recommends PFI/DFI for BRCA and caution for OS. Present survival as exploratory.
- **Published comparator:**
  - Preliminary RPPA Spearman: ESR1~ERα 0.89, BCL2 0.84, GATA3 0.79, PGR 0.76, CCNB1 0.73, AR 0.67, ERBB2 0.64, CDH1 0.54, EGFR 0.49.
  - Mertins 2016 (PMID 27251275) and Krug 2020 show basal vs luminal separation at the proteome level.
  - TCGA survival should be weak: smoke test Basal vs LumA OS p = 0.55, PFI p = 0.10; median follow-up 17 months (TCGA 2012).

---

## 4. Differentiating extra and optional stretch

### Extra: external validation in METABRIC
METABRIC has 1,980 tumours on Illumina arrays, with median OS_MONTHS of 116.5 and 646 deaths from breast cancer. It gives independent patients, a different platform, and ER/HER2 labels from different assays (IHC and SNP6). Its PR, HER2_STATUS and PAM50 columns are array-derived, so they are not independent.

1. **Classifier transfer (M4).**
   - **Freeze first.** Before opening any METABRIC label, freeze the gene list, hyperparameters and decision threshold on TCGA, and record them in a CHANGELOG entry tagged `metabric-freeze`.
   - Restrict to the harmonised shared genes.
   - Standardise each gene within each cohort.
   - Refit once on TCGA with the frozen hyperparameters (the data are already standardised, so drop the in-pipeline scaler).
   - Predict METABRIC:
     - ER vs ER_IHC (1,498 / 439).
     - HER2 vs HER2_SNP6 GAIN (438 / 1,537). TCGA's label is IHC/FISH positivity, METABRIC's is copy-number gain, and prevalence is 19% vs 22%.
     - ER−/HER2− vs the 286 group.
   - **AUC with bootstrap CI is the primary metric.** Give balanced accuracy and MCC at both the frozen threshold and a prevalence-matched threshold.
   - Robustness check: map each METABRIC gene onto TCGA's quantiles instead of z-scoring.
   - Also test M3's full-TCGA DE list as a fixed feature set. That is legitimate here because METABRIC is a true external test.
   - Report agreement with METABRIC ER_STATUS and HER2_STATUS only as "agreement with an expression-threshold call", next to the ESR1/ERBB2 single-gene baseline.
   - Secondary: 4-class PAM50 vs CLAUDIN_SUBTYPE (LumA 700, LumB 475, Her2 224, Basal 209), reported as macro-F1 and Cohen's kappa. This is agreement with array-derived labels, not ground truth.
   - Shah & Patel 2026 saw AUC fall 13–22% across cohorts. We report our own drop.
2. **DE replication (M3).**
   - ER-IHC−/HER2-not-gained (286) vs ER-IHC+/HER2-not-gained (1,220).
   - OLS on log2 intensity with age (statsmodels), BH-corrected.
   - Compare with the matched TCGA contrast ER−/HER2− vs ER+/HER2−. Report:
     - the share of TCGA DE genes that replicate with the same sign,
     - the Spearman correlation of effect sizes,
     - agreement of hallmark NES (prerank on the t-statistic).
3. **Survival (M2).**
   - Endpoints: OS, and breast-cancer death. Model breast-cancer death as a cause-specific hazard, with the 496 other-cause deaths censored.
   - Truncate follow-up at 10 years, because ER− hazards are high early and ER+ recurrences come late, so hazards cross over 355 months.
   - Fit a separate first-5-years model to compare with Nielsen 2010.
   - Kaplan-Meier and log-rank by IHC group, by predicted ER−/HER2− group, and by TCGA cluster (assigned by nearest centroid).
   - Cox model adjusted for age, grade and tumour size, with a Schoenfeld check. If PH fails, stratify or split time.
   - Compare the ordering with Sørlie 2001 (LumA best, basal worst).
   - Compare directions with Nielsen 2010 (HRs vs LumA: LumB 1.99, HER2E 3.65, basal 17.71). That cohort was ER+ and tamoxifen-treated, with disease-specific survival over 5 years, so compare direction and order, not size.

**Why it is low-risk**
- The clinical tables are downloaded and checked; the expression file is confirmed reachable.
- The code reuses the frozen TCGA pipeline.
- The go/no-go is Sat 31 Oct.
- Fallback: CPTAC 2020 for classifier transfer only (122 tumours, RNA-seq, clinical ER 81/39, ERBB2 amplification 17; harmonise genes first).

### Optional stretch: heterogeneity within HR+ tumours
**Cancelled unless METABRIC is finished by Wed 4 Nov.** Report only.
- Re-cluster the HR+/HER2− tumours into 2 groups.
- Characterise them with GSEA. Netanely 2016 found a proliferative, ductal group R1 and an immune, lobular group R2 with lower 5-year recurrence (p = 0.0076). Use the histology labels here.
- Project the split onto METABRIC and test relapse-free survival.

**Not recommended**
- Deep learning (e.g. a CNN). Mostavi's CNN reached 88.42%, about the same as Cascianelli's regularised logistic regression (88%). It is not in the course methods and it costs time.
- GTEx and SCAN-B.

---

## 5. Team split

| Member | Workstream | Outputs (files) | Handoffs | Shared role | Presents |
|---|---|---|---|---|---|
| **M1** | Data and integration | `labels.tsv`, `tumour_log2.parquet` (unscaled), `counts_int.parquet` (vial-deduplicated), `rppa.parquet`, `gene_map.tsv` (HGNC harmonisation), QC report, sample-flow figure (v1 by Fri 23 Oct); METABRIC matrices and labels plus the PAM50-coverage check (by Thu 29 Oct) | Data v1 to everyone (23 Oct); METABRIC data to M2, M3 and M4 (29 Oct) | **Code integrator**: GitHub repo, pull-request reviews, `run_all.py`, final clean-clone run. Compiles the AI-use logs. | Slides 1–3 |
| **M2** | Exploration, clustering, **survival** (TCGA and METABRIC) | Exploration figures; `clusters.tsv` and centroids; ARI/NMI with CIs; PAC plot; KM and Cox results with PH checks | Centroids to M4 for METABRIC cluster assignment (28 Oct) | **Slide-deck owner**: template, figure style (one fixed colour per clinical group), deck assembly | Slides 4–5 |
| **M3** | DE, enrichment, **METABRIC DE replication** | `de_tnbc_vs_hr.tsv` (threshold test, shrunk LFC), IDC-only sensitivity, matched ER−/HER2− contrast, GSEA (with leading edges), ORA, concordance tables, METABRIC replication | Top genes and shrunk LFC to M5 for RPPA (Wed 28 Oct); receives M4's stable features for ORA | **Report editor**: outline, merging, page limit | Slides 6–7 |
| **M4** | Classification and transfer | `cv_metrics.tsv`, corrected t-tests, ablation table, out-of-fold predictions, stable features, frozen models (`models/*.joblib`), `metabric-freeze` CHANGELOG entry, METABRIC/CPTAC predictions | OOF predictions to M5 for error analysis (Thu 29 Oct); METABRIC predictions to M2 (Sat 31 Oct) | **Statistics QA**: leakage audit, BH usage, seeds. Rehearsal timekeeper. Leads methods questions in Q&A. | Slide 8 |
| **M5** | Validation and literature | Benchmark table, `rppa_tests.tsv`, mRNA–protein and concordance figures, error-analysis panel | Benchmark numbers to every owner for their slides | **References** (PMIDs, DOIs; checks the partly verified numbers), intro, abstract, **Q&A bank** | Slides 9–10 |

**Load balance**
- M4's grid drops from 162 configurations to 16, plus 5 small extra runs.
- Survival moves to M2, because it uses M2's clusters. METABRIC DE replication moves to M3, because it uses the same code. The Q&A bank moves to M5.
- M5 starts on day 1 with RPPA and the literature table, which depend on no one.

**Backups:** M2 covers M1, M3 covers M2, M4 covers M3, M5 covers M4, and M1 covers M5.

**Working rules**
- One notebook per workstream (`01_data` … `05_validation`), plus a shared `src/`.
- Only M1 writes to `data/processed/`.
- Seed = 42. Figures saved as `figN_name.png` at 300 dpi.
- Daily 10-minute stand-up in the group chat.
- From day 1, everyone saves AI prompts and chat exports to `ai_logs/<member>/<date>_<topic>.md`.

---

## 6. Timeline

**Calendar check (at kick-off).**
- Diwali is around **Sun 8 Nov 2026**.
- IIT Madras July–Nov end-semester exams probably overlap the report weeks. Check the academic calendar on Wed 21 Oct.
- If exams start before 20 Nov, move every Phase 2 date earlier by the overlap.

### Phase 1: presentation (worst case Tue 3 Nov)

| Date | Milestone | Who |
|---|---|---|
| Wed 21 Oct | Kick-off (1 h): agree this plan, the label rules, the figure style, and the calendar check. Set up the repo and `requirements.txt`. Commit the helper scripts from the appendix. Everyone installs and downloads the core files and the HGNC set. M1 starts the METABRIC download with `fetch_cbio_lfs.py`. M5 starts the benchmark table and the Q&A file. | All |
| Thu 22 Oct | Labels table v0: ER/PR/HER2 rules, vial dedup, histology, stage collapse, kappa between label sources, STAR group re-tabulation. M2 and M4 prototype on raw HiSeqV2. M3 tests the count conversion. M5 joins the RPPA data. | M1, all |
| Fri 23 Oct | **Data freeze v1** (unscaled), gene map v1, sample-flow figure, batch PCA. | M1 |
| Sat 24 – Mon 26 | Round 1: PCA/UMAP and clustering (M2). PyDESeq2 threshold test plus IDC-only run (M3). Nested CV with baselines (M4). RPPA tests (M5). TCGA KM curves (M2). METABRIC harmonisation (M1). | All |
| **Tue 27 Oct** | **Checkpoint 1** (1 h): each member shows one key figure with its number next to the literature value. M4 runs the leakage audit on every notebook. **Scope freeze for the talk.** | All |
| Wed 28 – Thu 29 | Round 2: GSEA, ORA and concordance (M3). Ablations, corrected t-tests and PAM50 benchmark (M4). Consensus clustering and PAC (M2). RPPA concordance and error analysis (M5). METABRIC data and PAM50 coverage ready (M1, Thu). Everyone drafts slides and speaker notes. | All |
| Fri 30 Oct | Deck v1 assembled. **Dry run 1** (timed at 8 min). `metabric-freeze` entry written, then the first METABRIC transfer. | M2, M4 |
| **Sat 31 Oct** | **METABRIC go/no-go**: if ready, slide 9 shows the transfer AUC; if not, slide 9 shows RPPA and lists METABRIC under "for the report". Figures fixed. Q&A bank of at least 20 questions. **Report sections assigned.** | All |
| Sun 1 Nov | **Dry run 2** with 2 minutes of mock Q&A. Deck frozen and exported as PDF and PPTX to Drive and two devices. | All |
| Mon 2 Nov | Buffer and final run-through. METABRIC work continues alongside. | All |
| **Tue 3 Nov** | **Presentation** (the real slot is between 3 and 6 Nov). | All |

**Ready at talk time** (this goes on slide 10 as a "done / for the final report" line)
- **Done:**
  - Data preparation with every n.
  - Exploration and clustering with ARI.
  - DE with the threshold test and the histology check.
  - GSEA, ORA and concordance.
  - Nested-CV classification with baselines, ablation and corrected tests.
  - RPPA validation.
  - METABRIC transfer AUC, if the 31 Oct go/no-go passes.
- **Planned methods done:** about 70%, comfortably above the 40% rule.
- **For the report:** METABRIC DE replication and survival, TCGA Cox, the error analysis, and the stretch goal if it runs.

### Phase 2: report (due Wed 25 Nov)

| Date | Milestone |
|---|---|
| Wed 4 Nov | Stretch go/no-go: only if METABRIC is finished. Otherwise the stretch is cancelled. |
| Wed 4 – Sat 7 Nov | Act on Q&A feedback. Finish METABRIC (transfer, DE replication, survival), TCGA Cox and ORA on the classifier features. |
| **Sun 8 Nov** | **Results freeze** (Diwali; finish on Sat 7 if people are travelling). |
| Mon 9 – Tue 10 | M1 reruns everything from a clean clone. Final set: at most 6 figures and 2 tables. |
| Mon 9 – Thu 12 | Each owner writes their Methods and Results subsection (about 1 page), starting from their speaker notes. M5 writes the intro and literature; M3 writes the discussion skeleton. |
| **Fri 13 Nov** | **Report v1** merged by M3. |
| Sat 14 – Wed 18 | Cross-review (each member reviews another's section). Trim to 8 pages, Arial 11, single spacing. Check every number against the code outputs, and check the references. Ask the instructor whether references count toward the 8 pages. |
| Thu 19 – Fri 20 | Supplementary material: tagged GitHub release (code), README with reproduction steps, and the AI logs via a private Drive link (not in the public repo). Add the AI-use statement and the METABRIC ODbL attribution. |
| **Sat 21 – Sun 22 Nov** | **Submit.** Earlier if exams start sooner. |
| Wed 25 Nov | Deadline. |

**Page budget**

| Section | Pages |
|---|---|
| Title and abstract | 0.3 |
| Introduction | 0.7 |
| Data and methods | 1.5 |
| Results (5–6 figures) | 3.5 |
| Discussion (literature validation, biology, limitations) | 1.2 |
| References, data availability (ODbL), AI-use statement | 0.8 |

---

## 7. Ten-slide outline (8 minutes; total 7:40)

Speakers present in contiguous blocks, which gives 4 handovers.

| # | Slide | Content | Main visual | Time | Presenter |
|---|---|---|---|---|---|
| 1 | Title and question | Title, primary question, team | — | 0:15 | M1 |
| 2 | Why it matters and how we avoid circularity | Receptor tests vs PAM50. PAM50 labels are circular on this RNA-seq. Our design: IHC/FISH labels, protein data, and METABRIC with ER by IHC and HER2 by copy number. | 3-box design diagram | 0:40 | M1 |
| 3 | Data and preparation | Label rules; vial dedup; Normal-like calls are mostly adjacent normals; QC and batch check; gene harmonisation | Flow: 1,218 HiSeqV2 → 1,097 tumours (1,095 in STAR) → 949 complete → 605 / 184 / 160; RPPA 747; METABRIC 1,980 (1,932 with both labels) | 0:50 | M1 |
| 4 | Exploration | Main axes of variation; markers by clinical group | PCA/UMAP coloured by group, plus 4 violins | 0:45 | M2 |
| 5 | Clustering vs labels | k chosen by PAC; ARI/NMI with CIs vs PAM50, IHC and TCGA 2012; where LumA and LumB mix | Contingency heatmap, ARI next to Netanely and Huo | 0:50 | M2 |
| 6 | Differential expression | DEGs with \|LFC\| > 1 at FDR 5%, adjusted for histology; how many survive IDC-only; top genes other than the receptors vs literature | Volcano plot (shrunk LFC) | 0:50 | M3 |
| 7 | Pathways and cross-assay agreement | Hallmark NES with leading-edge genes; TNBC~Basal kappa and overlap vs TCGA 2012 | NES bar chart and a 2×2 table | 0:45 | M3 |
| 8 | Predicting clinical status | One panel with three numbers: single-gene AUC, best-model AUC (corrected-test p), AUC after removing receptor, amplicon and PAM50 genes. Shuffled-label MCC ≈ 0. | ROC panel with the 3 AUCs | 1:05 | M4 |
| 9 | Validation | RPPA concordance (permutation p); error-analysis example; METABRIC transfer AUC if go | mRNA vs protein scatter, plus AUC | 1:00 | M5 |
| 10 | Conclusions | Three take-home findings, limitations, "done / for the final report" line | 3 bullets | 0:40 | M5 |

Methods and results take slides 4–9, about 65% of the talk. No backup slides unless the instructor allows them; keep backup answers in M5's Q&A bank. Spread Q&A answers across members to even out credit. M4 takes methods questions.

**Likely questions** (M5 collects them; owners write the answers)
- Why not just use PAM50?
- How was HER2 defined?
- Why are the METABRIC groups not called TNBC?
- Why test |LFC| > 1 instead of filtering afterwards?
- Why BH rather than Bonferroni?
- How did you prevent leakage?
- Why does ER prediction still work without ESR1?
- Why remove the whole amplicon?
- Why is TCGA survival not significant?
- How did you handle the RNA-seq vs microarray difference, and the change in prevalence?

---

## 8. Top risks and mitigations

| Risk | Mitigation | Owner |
|---|---|---|
| Examiner calls the analysis "circular" | IHC/FISH labels drive DE and classifiers; PAM50 is only a comparator; clustering and ablations exclude the PAM50 genes; METABRIC uses only ER_IHC and HER2_SNP6. Shown on slide 2. | M5, M4 |
| METABRIC labels are array-derived (ER_STATUS, PR_STATUS, HER2_STATUS, PAM50) | Use only ER_IHC and HER2_SNP6. Name the contrast ER-IHC−/HER2-not-gained (286). Report the array-derived columns only as agreement with an expression call. | M1, M4 |
| Leakage from DE lists or global scaling | Unscaled freeze; all preprocessing inside the sklearn Pipeline; DE lists used only for interpretation and external tests; leakage audit on 27 Oct. | M4 |
| Duplicate STAR vials crash PyDESeq2 or create pseudo-replicates | Keep vial 01A for the 11 patients; log it; assert unique IDs in code. | M1 |
| Gene ID mismatches silently drop genes | HGNC harmonisation table for every dataset and the GMT; supplementary coverage table; PAM50 check in METABRIC before go/no-go. | M1 |
| HER2 classifier just learns co-amplified genes | Coordinate-based amplicon removal (±1.5 Mb, checked at ±2 Mb). | M4 |
| Histology confounds the TNBC vs HR+ contrast | Histology in the design; IDC-only sensitivity run. | M3 |
| Missing HER2 data or a debatable label rule | Written rule; sensitivity analysis with the nature2012 labels; kappa between sources. | M1 |
| Overlapping CV folds make SDs look too small | Nadeau-Bengio corrected t-test; bootstrap CIs on out-of-fold AUC. | M4 |
| Proportional hazards fails in METABRIC | Schoenfeld check; 10-year truncation; separate 5-year model; stratify if needed. | M2 |
| Transfer metrics move with the threshold or prevalence | Freeze before looking; AUC primary; threshold metrics at both the frozen and the prevalence-matched threshold. | M4 |
| METABRIC download fails | Start on Wed 21 Oct with the sha256-checking script; git-lfs fallback; CPTAC 2020 fallback. | M1 |
| PyDESeq2 runtime on 765 samples | Filter genes first and set `n_cpus`. Fallback: balanced 150 vs 150 subsamples, repeated 3 times. The 60-sample test took 70 s. | M3 |
| Package conflicts | Python 3.12, pandas 2.3.x, one `requirements.txt`, install on day 1. | M1 |
| Diwali and end-semester exams eat working days | Calendar check on 21 Oct; results freeze 8 Nov; report v1 13 Nov; submit 21–22 Nov; stretch cancelled unless METABRIC is done by 4 Nov. | All |
| Talk runs over time | Scope freeze 27 Oct; one message per slide; two timed dry runs. | M2, M4 |
| Partly verified literature numbers (Huo & Tseng ARIs, Parker per-subtype HRs) | Check against the paper before citing; otherwise cite the finding without the number. | M5 |
| Hallmark GMT is a third-party mirror | Cite MSigDB v2024.1 and Liberzon 2015; compare once with the official file. | M3 |
| AI-use documentation missing | Log from day 1, including `ai_logs/planning/`; AI-use statement in the report. | All, M1 |

---

## Appendix: helper scripts

These are in `scripts/` of this repository. They use only the Python standard library or pandas.

- `fetch_cbio_lfs.py`: downloads cBioPortal datahub files stored in Git LFS and checks sha256. This is the main METABRIC route.
- `labels.py`, `her2.py`, `groups.py`: rebuild the ER/PR/HER2 labels, the HER2 rule and the clinical-group counts in section 1 from `BRCA_clinicalMatrix`.
- `smoke.py`: a small PyDESeq2 and gseapy test (30 Basal vs 30 LumA) that checks the install works.

---

## Changes made after review

1. **METABRIC labels:** switched to ER_IHC (1,498/439) and HER2_SNP6 GAIN (438/1,537). The contrast is renamed ER-IHC−/HER2-not-gained (286 vs 1,220), with a matching TCGA ER−/HER2− target and DE contrast. The "independent labels" claim is narrowed to "other assays", and the array-derived columns are reported only as agreement.
2. **Leakage:** the freeze ships an unscaled matrix; all preprocessing runs inside CV folds; M3's DE lists are not used as CV features (only for interpretation and external tests); leakage audit added.
3. **Classification scope:** the grid drops from 162 configurations to 16 (2 targets × 4 models × 2 ablations; k tuned in the inner CV). The PAM50 benchmark uses L1-LR only, in 4-class and 5-class versions.
4. **Rebalanced team:** survival → M2, METABRIC DE replication → M3, Q&A bank → M5. The paired tumour-vs-normal DE run is dropped.
5. **DE statistics:** threshold test (`lfc_null=1`, `greaterAbs`), shrunk LFCs for plots and rankings, histology in the design, IDC-only sensitivity run (402 vs 136).
6. **STAR samples:** documented vial-01A dedup (11 patients → 1,095 tumours). Corrected to 113 normals and pairs, and to STAR group sizes (TNBC 160).
7. **Gene IDs:** HGNC harmonisation table for all datasets and the GMT, a coverage table, a PAM50 check in METABRIC, and defined ORA universes.
8. **Amplicon ablation:** removal by coordinates (±1.5 Mb around ERBB2, about 87 HiSeqV2 genes, checked at ±2 Mb). Post-ablation ER features are interpreted as the ER regulon.
9. **Tests:** Nadeau-Bengio corrected t-test and bootstrap AUC CIs; odds ratio, kappa and % overlap instead of bare Fisher p-values; RPPA concordance by Spearman with a label-permutation null and phospho forms collapsed; NES with leading edges plus a phenotype-permutation check.
10. **Survival:** PH checks, 10-year truncation and a 5-year model, cause-specific hazard for breast-cancer death, TCGA stage collapsed with a covariate budget.
11. **Transfer design:** freeze-before-looking CHANGELOG entry, a single refit on shared genes, AUC as the primary metric, metrics at both thresholds, and a quantile-to-TCGA check.
12. **Calendar:** Diwali and exam check added; results freeze moved to 8 Nov, report v1 to 13 Nov, submission to 21–22 Nov; stretch gated on 4 Nov; report sections assigned on 31 Oct.
13. **Scripts and download:** the helper scripts are in the shared project files for a day-1 commit; `fetch_cbio_lfs.py` is the primary METABRIC route; the git-lfs fallback uses a local config.
14. **Low-severity fixes:**
    - PanCan Subtype_mRNA is the main PAM50 comparator (all 1,097 tumours), with bootstrap CIs on ARI.
    - IHC-vs-RNA error-analysis panel added.
    - Contiguous presenter blocks (4 handovers) and a "done / for the final report" line.
    - Slide 8 merged into one panel.
    - Like-for-like 5-class benchmark.
    - `ai_logs/planning/`, an AI-use statement and the ODbL attribution.