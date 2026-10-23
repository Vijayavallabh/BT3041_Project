# Q&A prep (2 minutes after the talk)

Answer in one or two sentences and then stop. Numbers are for the final labels, signed off on 23 Oct (HER2 rule v1 following
ASCO/CAP 2018, no filter: HR+/HER2− 591, HER2+ 160, TNBC 157); † marks numbers that would move if the labels changed. Suggested owner in brackets (Navya, M4, takes methods questions, per the plan).

## Design and labels

**1. Why not just use the PAM50 subtypes?** [Vijayavallabh]
PAM50 is computed from the same expression data, so asking whether expression separates PAM50 would be circular. We
use the clinical IHC/FISH results as labels and PAM50 only for comparison (TNBC vs Basal κ 0.78, close to the 2012
TCGA paper).

**2. How did you define HER2 status?** [Vijayavallabh]
We followed ASCO/CAP 2018. FISH comes first, then the IHC score (3+ positive, 0/1+ negative). We set three cases to
missing: 2+ without FISH, IHC 3+ with a negative FISH (5 tumours, which ASCO/CAP says should be re-reviewed), and any
contradiction between status and score. The rule agrees with the TCGA 2012 calls at κ 0.88. The noisy calls are the
ones based only on status text; dropping the low-confidence HER2+ calls raises ERBB2 alone from 0.845 to 0.894 AUC †.

**3. Why were 9 tumours excluded?** [Vijayavallabh]
Five failed an identity check because their HiSeqV2 and STAR profiles did not match, and five had redacted clinical
records in TCGA-CDR. One tumour had both problems.

**4. Why are the METABRIC groups not called TNBC?** [Anitha]
METABRIC has no IHC PR, and its HER2 label is DNA copy number, so the accurate name is ER−/HER2-not-gained. We used
only ER by IHC and HER2 by SNP6, because METABRIC's other status columns are derived from the arrays themselves.

**5. Why were male patients left out?** [Vijayavallabh]
The 12 male tumours are in the unsupervised analyses but not in the supervised set. Receptor biology and labels
differ in male breast cancer, and 12 tumours are too few to model.

## Methods

**6. Why test |log2FC| > 1 instead of filtering by fold change afterwards?** [Harshit]
With 748 samples almost everything is significant; the plain Wald test calls 79% of genes. Testing the threshold
directly (lfc_null = 1) controls the FDR for "more than 2-fold", and filtering afterwards does not.

**7. Why Benjamini–Hochberg and not Bonferroni?** [Harshit]
We want a candidate gene list with a controlled false discovery rate, not a guarantee of zero false positives.
Bonferroni over 23,000 genes would throw away real effects.

**8. Why does only one Hallmark set survive CAMERA?** [Harshit]
Gene-permutation GSEA treats genes as independent, so co-expressed sets (E2F, inter-gene ρ 0.25) look more
significant than they are. CAMERA and label-permutation GSEA agree on the direction of all 28 sets, but only early
estrogen response stays at q < 0.05. We therefore report direction and effect size and do not claim q ≈ 0.

**9. How did you prevent data leakage in classification?** [Navya]
We used nested CV. Gene filtering, scaling, feature selection and tuning are refitted inside each training fold, and
the DE lists are never used as features. Shuffled labels give MCC ≈ 0 (−0.01 ER, +0.04 HER2).

**10. If models do not beat ESR1, why build them?** [Navya]
That is itself the finding: the receptor gene already carries almost all of the clinical test result. The models are
useful in the ablation. ER is still predicted at 0.95 without ESR1, PGR, the amplicon and PAM50, so ER status is spread
across a regulon, whereas HER2 falls to 0.78–0.82.

**11. Why remove the whole amplicon and not just ERBB2?** [Navya]
Genes next to ERBB2 (GRB7, STARD3, PGAP3…) are co-amplified, so a model could "predict HER2" from any of them. We
removed every gene within 1.5 Mb and also checked 2 Mb, with the same result.

**12. Why the corrected t-test?** [Navya]
CV folds share training data, so the ordinary paired t-test is overconfident. The Nadeau–Bengio correction inflates
the variance to allow for that overlap, and a paired bootstrap is our second check.

**13. Why does HER2 do better in METABRIC (0.91) than in TCGA CV (0.85)?** [Anitha]
METABRIC's HER2 label is DNA copy number, which follows ERBB2 mRNA more closely than IHC/FISH does, and ERBB2 alone
also rises to 0.90 there. So it is most likely a label effect, not a better model.

**14. How did you handle RNA-seq vs microarray?** [Navya]
The models use within-cohort z-scores on the 18,936 shared genes and were frozen on TCGA before any METABRIC label was
read. Quantile mapping to TCGA instead changes AUC by less than 0.012. AUC is the main metric because it depends on
neither the threshold nor the prevalence.

**15. Did you adjust for tumour purity?** [Harshit]
No. Part of the immune and proliferation difference is probably tissue composition. We list this as a limitation and
as future work.

**16. Why is k = 2 the chosen number of clusters?** [Vedant]
Consensus clustering gives the lowest PAC at k = 2 (0.02 vs 0.14 at k = 5). The k = 5 HER2-enriched cluster is the
least stable one (Jaccard 0.40), so we only call it a tendency.

## Validation and biology

**17. Are the "misclassified" tumours really errors?** [Anitha]
Not all of them. 8 of the 10 IHC-ER− tumours that the model calls ER+ have ERα protein above the 95th percentile of
correctly called ER− tumours, so they look like IHC false negatives.

**18. Does any of this relate to outcome?** [Vedant] (backup figure: `../step3_analysis/figures/slideS1_survival.png`)
In TCGA, PFI differs by group (log-rank p = 0.007), with TNBC at HR 2.69 (1.70–4.26) after adjusting for age and
stage. With only 100 events and 2.4 years of median follow-up, this is exploratory. In METABRIC, ER−/HER2-not-gained
tumours have higher breast-cancer death in years 0–5 (HR 2.46) and lower after that (0.51), which is the known
early-relapse pattern.

**19. How do your numbers compare with published work?** [Anitha]
TNBC called Basal: 83% vs 86% (TCGA 2012). Single-gene accuracy: ESR1 94% / ERBB2 89% vs 88–96% / 89–93% (Gong 2007).
The full 21-row table with DOIs is in `../step3_analysis/results/literature_table.md`.

**20. What is the biological take-home?** [Anitha]
TNBC is a distinct state, with the ER regulon down and basal factors and proliferation up. HER2+ tumours are mostly
luminal tumours with an amplicon, which is why expression separates them much less.

**21. How did you use AI?** [Vijayavallabh]
For planning, code and review. Every prompt and chat log is in the supplementary material submitted with the report,
as the guidelines ask, and all code is public at https://github.com/Vijayavallabh/BT3041_Project.
