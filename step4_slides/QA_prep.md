# Q&A prep (2 minutes after the talk)

Answer in one or two sentences, then stop. Numbers are for the current labels (HER2 rule v1, no filter: HR+/HER2− 591,
HER2+ 160, TNBC 157); if the team changes the labels and reruns step 3, recheck the numbers marked with † against the
step 3 README. Suggested owner in brackets (M4 takes methods questions, per the plan).

## Design and labels

**1. Why not just use the PAM50 subtypes?** [M1]
PAM50 is computed from the same expression data, so asking whether expression separates PAM50 is circular. We use the
clinical IHC/FISH results as labels and PAM50 only as a comparator (TNBC vs Basal κ 0.78, close to the 2012 TCGA paper).

**2. How did you define HER2 status?** [M1]
FISH first; otherwise the IHC score (3+ positive, 0/1+ negative); 2+ without FISH, and any contradiction between
status and score, set to missing. It agrees with the TCGA 2012 calls at κ 0.88. Calls based only on the status text
are the noisy ones: dropping the low-confidence HER2+ calls raises ERBB2 alone from 0.845 to 0.894 AUC †.

**3. Why were 9 tumours excluded?** [M1]
Five failed an identity check (their HiSeqV2 and STAR profiles did not match each other), five had redacted clinical
records in TCGA-CDR, one had both.

**4. Why are the METABRIC groups not called TNBC?** [M5]
METABRIC has no IHC PR, and its HER2 label is DNA copy number, so the honest name is ER−/HER2-not-gained. We used only
ER by IHC and HER2 by SNP6; METABRIC's other status columns are derived from the arrays themselves.

**5. Why were male patients left out?** [M1]
12 male tumours are in the unsupervised analyses but not in the supervised set, because receptor biology and labels
in male breast cancer differ and 12 is too few to model.

## Methods

**6. Why test |log2FC| > 1 instead of filtering by fold change afterwards?** [M3]
With 748 samples almost everything is significant: the plain Wald test calls 79% of genes. Testing the threshold
directly (lfc_null = 1) gives FDR control for "more than 2-fold", which post-hoc filtering does not.

**7. Why Benjamini–Hochberg and not Bonferroni?** [M3]
We want a list of candidate genes with a controlled false discovery rate, not a guarantee of zero false positives;
Bonferroni over 23,000 genes would throw away real effects.

**8. Why does only one Hallmark set survive CAMERA?** [M3]
Gene-permutation GSEA assumes genes are independent; co-expressed sets (E2F, inter-gene ρ 0.25) look more significant
than they are. CAMERA and label-permutation GSEA agree on the direction of all 28 sets, but only early estrogen response
stays at q < 0.05. So we report the direction and effect size and do not claim q ≈ 0.

**9. How did you prevent data leakage in classification?** [M4]
Nested CV: gene filtering, scaling, feature selection and tuning are refitted inside each training fold; the DE lists
are never used as features. Shuffled labels give MCC ≈ 0 (−0.01 ER, +0.04 HER2).

**10. If models do not beat ESR1, why build them?** [M4]
That is the finding: the clinical test result is already almost fully written in the receptor gene. The models become
useful in the ablation: ER is still predicted at 0.95 without ESR1, PGR, the amplicon and PAM50, so ER status is a
regulon-wide state, whereas HER2 falls to 0.78–0.82.

**11. Why remove the whole amplicon and not just ERBB2?** [M4]
Genes next to ERBB2 (GRB7, STARD3, PGAP3…) are co-amplified, so a model could "predict HER2" from any of them. We removed
every gene within 1.5 Mb and checked 2 Mb: the result is the same.

**12. Why the corrected t-test?** [M4]
CV folds share training data, so the ordinary paired t-test is too confident. The Nadeau–Bengio correction inflates
the variance for that overlap; a paired bootstrap is the second check.

**13. Why does HER2 do better in METABRIC (0.91) than in TCGA CV (0.85)?** [M5]
METABRIC's HER2 label is DNA copy number, which tracks ERBB2 mRNA more closely than IHC/FISH; ERBB2 alone also rises
to 0.90 there. It is a label effect, not a better model.

**14. How did you handle RNA-seq vs microarray?** [M4]
Models use within-cohort z-scores on the 18,936 shared genes and were frozen on TCGA before any METABRIC label was
read. Quantile mapping to TCGA instead changes AUC by less than 0.012. AUC is the main metric because it does not
depend on the threshold or the prevalence.

**15. Did you adjust for tumour purity?** [M3]
No. Some of the immune and proliferation difference is probably tissue composition; that is a stated limitation and
future work.

**16. Why is k = 2 the chosen number of clusters?** [M2]
Consensus clustering gives the lowest PAC at k = 2 (0.02 vs 0.14 at k = 5). The k = 5 HER2-enriched cluster is the
least stable one (Jaccard 0.40), so we call it a tendency.

## Validation and biology

**17. Are the "misclassified" tumours really errors?** [M5]
Partly not: 8 of 10 IHC-ER− tumours the model calls ER+ have ERα protein above the 95th percentile of correctly called
ER− tumours, so they look like IHC false negatives.

**18. Does any of this relate to outcome?** [M2] (backup figure: `../step3_analysis/figures/slideS1_survival.png`)
In TCGA, PFI differs by group (log-rank p = 0.007); TNBC HR 2.69 (1.70–4.26) adjusted for age and stage, but with only
100 events and 2.4 years median follow-up it is exploratory. In METABRIC, ER−/HER2-not-gained has higher breast-cancer
death in years 0–5 (HR 2.46) and lower after (0.51), the known early-relapse pattern.

**19. How do your numbers compare with published work?** [M5]
TNBC called Basal 83% vs 86% (TCGA 2012); single-gene accuracy ESR1 94% / ERBB2 89% vs 88–96% / 89–93% (Gong 2007);
the 21-row table with DOIs is in `../step3_analysis/results/literature_table.md`.

**20. What is the biological take-home?** [M5]
TNBC is a distinct state (ER regulon off, basal factors and proliferation on); HER2+ tumours are mostly luminal tumours
with an amplicon, which is why expression separates them far less.

**21. How did you use AI?** [M1]
For planning, code and review; every prompt and chat log is in the supplementary material submitted with the report, as the guidelines
ask. All code is public at https://github.com/Vijayavallabh/BT3041_Project.
