# Speaker script

Total 7:45 of the 8:00 limit; methods and results (slides 4–9) take 5:10, 67% of the talk (the guidelines ask for at least 40%). About 130 words a minute.

| Slide | Presenter | Time |
|---|---|---|
| 1 | Vijayavallabh | 0:20 |
| 2 | Vijayavallabh | 0:45 |
| 3 | Vijayavallabh | 0:45 |
| 4 | Vedant | 0:40 |
| 5 | Vedant | 0:45 |
| 6 | Harshit | 0:50 |
| 7 | Harshit | 0:50 |
| 8 | Navya | 1:05 |
| 9 | Anitha | 1:00 |
| 10 | Anitha | 0:45 |

## Slide 1 (Vijayavallabh, 0:20)

Good morning. Our project asks what tumour gene expression tells us about the three clinical receptor groups of breast cancer. We use 1,088 TCGA tumours, and check each finding in protein from the same tumours and in an independent cohort, METABRIC.

## Slide 2 (Vijayavallabh, 0:45)

In the clinic, breast cancer is grouped by three tests: ER and PR by immunohistochemistry, and HER2 by IHC or FISH. These tests decide treatment, and triple-negative tumours have no targeted therapy. Many studies use the PAM50 subtypes instead, but PAM50 is computed from expression, so that would be circular. We use the clinical results as labels and keep PAM50 for comparison. We ask three questions. Does unsupervised structure recover the groups? Which genes and pathways differ? And can expression predict clinical status better than the receptor gene alone? We validate in protein from the same tumours and in METABRIC.

## Slide 3 (Vijayavallabh, 0:45)

Data preparation. We keep one primary tumour per patient and exclude 9 that failed an RNA identity check or had redacted clinical records. That leaves 1,088 tumours for unsupervised work. For HER2 we use FISH first, then the IHC score, and set equivocal or contradictory results to missing. Women with all three receptors known give 591 HR+/HER2−, 160 HER2+ and 157 triple-negative tumours. Differential expression uses raw STAR counts, keeping genes with at least 10 counts in 157 tumours. RPPA adds protein for 730 tumours, and METABRIC adds 1,980 tumours with independent labels.

## Slide 4 (Vedant, 0:40)

First, exploration without labels. On the UMAP of the 2,000 most variable genes, one island separates completely, and it holds 82% of all triple-negative tumours. HER2-positive tumours, in orange, mostly sit at one edge of the hormone-receptor-positive tumours. The violins explain this. ESR1 is the largest single difference between TNBC and HR+, with a rank-biserial effect of −0.95. ERBB2 marks HER2+ most clearly, and FOXC1, a basal transcription factor, and MKI67, for proliferation, are higher in TNBC. All four differences are significant after Benjamini–Hochberg correction.

## Slide 5 (Vedant, 0:45)

Next, clustering, with the PAM50 genes left out. Consensus clustering is most stable at k = 2, with a PAC of 0.02, and a basal-like group of 217 tumours splits off. Its agreement with IHC ER status is an adjusted Rand index of 0.67. At k = 5 a HER2-enriched cluster appears, but it is the least stable, so we call it a tendency. Luminal tumours also split differently from PAM50 LumA and LumB, a known continuum. So for question 1, expression recovers the ER-negative, basal axis strongly, and HER2 only weakly.

## Slide 6 (Harshit, 0:50)

Differential expression, TNBC against HR+/HER2−, with PyDESeq2 on raw counts, adjusting for histology. Rather than filtering by fold change after the test, we test the 2-fold threshold directly. A plain Wald test calls 79% of genes significant, which tells us little. With the threshold test, 2,741 genes differ more than 2-fold at 5% FDR. ESR1 and PGR are positive controls. Beyond them, the ER regulon (GATA3, FOXA1, AGR3) is lower in TNBC, and basal transcription factors (FOXC1, EN1, SOX10) are higher. On the right, the list is stable: it keeps 91% of genes without lobular tumours, 94% with age added, and 87% with strict labels. HER2+ against HR+ gives only 148 genes, led by the 17q12 amplicon around ERBB2.

## Slide 7 (Harshit, 0:50)

Pathways. GSEA on the ranked genes finds 28 Hallmark sets at FDR 5%. TNBC scores higher on proliferation (E2F targets, G2M), MYC and immune programmes, and lower on estrogen response. Gene-permutation q-values are optimistic when genes are co-expressed, so we ran two tests that account for correlation: CAMERA and label-permutation GSEA. Both agree on the direction of every set, but only early estrogen response stays significant. So we report effect sizes, not q ≈ 0. Top right: clinical TNBC agrees with PAM50 Basal with kappa 0.78, close to the 2012 TCGA paper. Bottom right: the same pathway shifts replicate in METABRIC, with a correlation of 0.94.

## Slide 8 (Navya, 1:05)

Question 3: can expression predict clinical status? We trained four models (L1 logistic regression, linear SVM, Random Forest and PLS-DA) in nested cross-validation, with every filtering, scaling and selection step inside the training folds, so nothing leaks. Shuffled labels give an MCC of about zero. The main result: the models match the receptor gene alone but do not beat it. ESR1 alone has an AUC of 0.959, the best model 0.954; ERBB2 alone 0.845, the best model 0.851. Corrected t-tests, which allow for overlapping folds, give p = 0.69 and 0.66. We then removed the receptor genes, every gene within 1.5 Mb of ERBB2, and the PAM50 genes. ER is unchanged, at 0.95–0.95, so the ER signal is spread across a regulon, not one gene. HER2 drops to 0.78–0.82. The HER2 signal left outside the amplicon mostly reflects the lower ER and higher proliferation of HER2+ tumours.

## Slide 9 (Anitha, 1:00)

Validation. Left: for 85 genes with a protein antibody, the mRNA fold change predicts the protein difference with Spearman rho 0.80; a permutation test gives p = 0.001. Middle: where does the ER model go wrong? Of the 15 tumours that IHC calls ER-negative but the model calls positive, 8 of the 10 with protein data have ER-alpha protein at ER-positive levels. These are probably IHC false negatives, not model errors. Right: METABRIC. We froze the models on TCGA before reading any METABRIC label. ER reaches AUC 0.96, HER2 0.91, and ER−/HER2− versus rest 0.93. HER2 does even better than in cross-validation, most likely because METABRIC's HER2 label is copy number. And 86% of the TCGA DEGs replicate with the same sign.

## Slide 10 (Anitha, 0:45)

Three take-home messages. First, triple-negative breast cancer is a distinct transcriptional state: ER regulon down, basal and proliferation programmes up. HER2-positive tumours differ much less, mainly at the amplicon. Second, the receptor gene alone carries almost all the information about clinical status. ER signal is spread over a regulon, while HER2 depends more on the amplicon. Third, these results hold in protein and in an independent cohort. Limitations: pathway significance shrinks once gene correlation is counted, we did not adjust for purity, and TCGA survival is exploratory. Thank you; we are happy to take questions.
