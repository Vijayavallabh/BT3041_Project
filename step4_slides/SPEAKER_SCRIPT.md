# Speaker script

Total 7:45 of the 8:00 limit; methods and results (slides 4–9) take 5:10, 67% of the talk (the guidelines ask for at least 40%). About 130 words a minute.

| Slide | Presenter | Time |
|---|---|---|
| 1 | M1 | 0:20 |
| 2 | M1 | 0:45 |
| 3 | M1 | 0:45 |
| 4 | M2 | 0:40 |
| 5 | M2 | 0:45 |
| 6 | M3 | 0:50 |
| 7 | M3 | 0:50 |
| 8 | M4 | 1:05 |
| 9 | M5 | 1:00 |
| 10 | M5 | 0:45 |

## Slide 1 (M1, 0:20)

Good morning. We ask what tumour gene expression reveals about the three clinical receptor groups of breast cancer. We use 1,088 TCGA tumours, and we check every finding in protein from the same tumours and in an independent cohort, METABRIC.

## Slide 2 (M1, 0:45)

Clinically, breast cancer is grouped by three tests: ER and PR by immunohistochemistry, and HER2 by IHC or FISH. These tests decide treatment, and triple-negative tumours have no targeted therapy. Many studies use the PAM50 subtypes instead, but PAM50 is computed from expression, so that would be circular. We use the clinical results as labels and PAM50 only as a comparator. Our three questions: does unsupervised structure recover the groups; which genes and pathways differ; and can expression predict clinical status better than the receptor gene alone. We validate in protein from the same tumours and in METABRIC.

## Slide 3 (M1, 0:45)

Data preparation. We keep one primary tumour per patient and exclude 9 that failed an RNA identity check or had redacted clinical records, leaving 1,088 tumours for unsupervised work. For HER2 we use FISH first, then the IHC score, and set equivocal or contradictory results to missing. Women with all three receptors known give 591 HR+/HER2−, 160 HER2+ and 157 triple-negative tumours. Differential expression uses raw STAR counts, keeping genes with at least 10 counts in 157 tumours. RPPA adds protein for 730 tumours, and METABRIC adds 1,980 tumours with independent labels.

## Slide 4 (M2, 0:40)

First, exploration without labels. On the UMAP of the 2,000 most variable genes, one island separates completely. It holds 82% of all triple-negative tumours. HER2-positive tumours, in orange, mostly sit among the hormone-receptor-positive tumours, at one edge. The violins show why: ESR1 is the strongest single difference between TNBC and HR+, with a rank-biserial effect of −0.95. ERBB2 is the clearest marker of HER2+, while FOXC1, a basal transcription factor, and MKI67, proliferation, are higher in TNBC. All four differences are significant after Benjamini–Hochberg correction.

## Slide 5 (M2, 0:45)

Next, clustering, with the PAM50 genes left out. Consensus clustering is most stable at k = 2, with a PAC of 0.02: a basal-like group of 217 tumours splits off. Its agreement with IHC ER status is an adjusted Rand index of 0.67. At k = 5 a HER2-enriched cluster appears, but it is the least stable, so we call it a tendency, and luminal tumours split differently from PAM50 LumA and LumB, a known continuum. So for question 1: expression recovers the ER-negative, basal axis strongly, HER2 only weakly.

## Slide 6 (M3, 0:50)

Differential expression, TNBC versus HR+/HER2−, with PyDESeq2 on raw counts, adjusting for histology. Instead of filtering by fold change after the test, we test the 2-fold threshold directly. A plain Wald test calls 79% of genes significant, which is not informative; with the threshold test, 2,741 genes differ more than 2-fold at 5% FDR. ESR1 and PGR are positive controls. Beyond them, the ER regulon (GATA3, FOXA1, AGR3) is lower in TNBC, and basal transcription factors (FOXC1, EN1, SOX10) are higher. On the right, the list is robust: it keeps 91% of genes if we drop lobular tumours, 94% adding age, and 87% with strict labels. HER2+ vs HR+ gives only 148 genes, led by genes on the 17q12 amplicon around ERBB2.

## Slide 7 (M3, 0:50)

Pathways. GSEA on the ranked genes finds 28 Hallmark sets at FDR 5%. TNBC ranks higher on proliferation (E2F targets, G2M), MYC and immune programmes, and lower on estrogen response. These gene-permutation q-values are known to be optimistic when genes are co-expressed, so we ran two tests that account for correlation: CAMERA and label-permutation GSEA. They agree on the direction of every set, but only early estrogen response stays significant. So we report effect sizes, not q ≈ 0. Top right: clinical TNBC agrees with PAM50 Basal with kappa 0.78, close to the 2012 TCGA paper. Bottom right: the same pathway shifts replicate in METABRIC, with a correlation of 0.94.

## Slide 8 (M4, 1:05)

Question 3: can expression predict clinical status? We trained four models (L1 logistic regression, linear SVM, Random Forest and PLS-DA) in nested cross-validation, with every filtering, scaling and selection step inside the training folds, so there is no leakage. Shuffled labels give an MCC of about zero. The key result: models match, but do not beat, the receptor gene alone. ESR1 alone has an AUC of 0.959, the best model 0.954; ERBB2 alone 0.845, the best model 0.851. The corrected t-tests, which account for overlapping folds, give p = 0.69 and 0.66. Then we removed the receptor genes, every gene within 1.5 Mb of ERBB2, and the PAM50 genes. ER is unchanged, at 0.95–0.95: ER status is written across a whole regulon, not one gene. HER2 drops to 0.78–0.82: some HER2 signal lies outside the amplicon, largely the lower-ER, proliferative character of HER2+ tumours.

## Slide 9 (M5, 1:00)

Validation. Left: for 85 genes with a protein antibody, the mRNA fold change predicts the protein difference with Spearman rho 0.80; a permutation test gives p = 0.001. Middle: where does the ER model err? Of the 15 tumours that IHC calls ER-negative but the model calls positive, 8 of the 10 with protein data have ER-alpha protein at ER-positive levels. These are probably IHC false negatives, not model errors. Right: METABRIC. We froze the models on TCGA before reading any METABRIC label. The ER model reaches AUC 0.96, HER2 0.91, and ER−/HER2− versus rest 0.93; HER2 transfers even better than it cross-validated, because METABRIC's HER2 label is copy number. And 86% of the TCGA DEGs replicate with the same sign.

## Slide 10 (M5, 0:45)

Three take-home messages. First, triple-negative breast cancer is a distinct transcriptional state: ER regulon off, basal and proliferation programmes on. HER2-positive tumours differ far less, led by the amplicon. Second, the receptor gene already carries almost all the information about clinical status; ER signal is spread over a regulon, HER2 depends more on the amplicon. Third, these results hold in protein and in an independent cohort. Limitations: pathway significance shrinks once gene correlation is counted, we did not adjust for purity, and TCGA survival is exploratory. Thank you; we are happy to take questions.
