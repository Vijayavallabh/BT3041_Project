# BT3041 term project: what tumour gene expression reveals about breast cancer receptor groups

Term project for BT3041 (Analysis and Interpretation of Biological Data), IIT Madras, 2026.
Team 16: Vijayavallabh, Vedant, Harshit, Navya and Anitha.

We ask what TCGA breast tumour RNA-seq shows about the **clinical receptor groups** (ER, PR and HER2 by IHC/FISH:
HR+/HER2−, HER2+ and triple-negative). We use those clinical labels rather than PAM50 subtypes, because PAM50 is
computed from the same RNA-seq and would make the analysis circular. Results are checked against TCGA protein data
(RPPA), the independent METABRIC cohort and published studies.

## Main results

- **Unsupervised structure.** The most stable clustering is k = 2: a basal-like group that holds most triple-negative
  tumours, against the rest. Assigning METABRIC tumours to the TCGA k = 2 centroids gives a similar basal-like group
  (15.5% of tumours, 91% ER-negative).
- **Differential expression.** 2,741 genes differ between triple-negative and HR+/HER2− tumours (PyDESeq2,
  |log2 fold change| > 1 tested directly, FDR 5%). The result is robust across sensitivity runs. In METABRIC, 86% of the
  DEGs measured on its array (1,599 of 1,857, ER−/HER2− vs ER+/HER2−) replicate with the same sign.
- **Classification.** Whole-transcriptome models match, but do not beat, the single receptor gene. For ER, ESR1
  alone gives AUC 0.959 and the best model (random forest) 0.954. For HER2, ERBB2 alone gives 0.845 and the best
  model 0.851 (nested cross-validation, Nadeau–Bengio corrected tests).
- **Validation.** Models frozen on TCGA transfer to METABRIC with AUC 0.959 (ER) and 0.906 (HER2, which METABRIC
  defines by copy-number gain rather than IHC/FISH). RNA and protein differences between TNBC and HR+/HER2− tumours
  agree across the 85 genes that have a protein antibody (Spearman ρ 0.80, permutation p = 0.001; ρ 0.58 for HER2+ vs
  HR+/HER2−).
- **Survival.** Triple-negative tumours have shorter progression-free interval in TCGA (Cox hazard ratio 2.69, 95% CI
  1.70–4.26, 100 events; follow-up is short).

Pathway caveat: gene-permutation GSEA flags 28 Hallmark sets, but only the early estrogen response survives CAMERA,
which accounts for gene–gene correlation. The report says so.

## Repository layout

| Folder | What it holds |
|---|---|
| `project_plan.md` | Team plan: research questions, five workstreams, timeline. |
| `scripts/` | Small helper scripts used while planning (label rules, METABRIC download helper). |
| `step2_preprocessing/` | Download with checksums, sample QC, receptor labels, normalisation, PCA/UMAP, clustering, METABRIC projection. Its `README.md` is the data dictionary and lists every decision. |
| `step3_analysis/` | Differential expression, pathway enrichment, classifiers, METABRIC transfer, RPPA protein validation, survival. `config.json` holds the label choices. |
| `step4_slides/` | Presentation (`BT3041_talk.pptx` / `.pdf`), its build script, speaker script and Q&A prep. |
| `step5_report/` | Final report, pre-submission draft (`BT3041_final_report.docx` / `.pdf`; roll numbers are still placeholders), figure and report build scripts, `report_numbers.md` (the key numbers used in the report, generated from the result files). |

`step2_preprocessing/` and `step3_analysis/` each have `scripts/`, `results/` (with a `params.json` of every parameter
and package version), `figures/`, `logs/` from the last run, a README and a CHANGELOG.

## Reproduce

Python 3.12 on Linux (tested), with `curl` for the downloads. LibreOffice and `pdfinfo` are only needed to rebuild the
slides and report PDFs.

```bash
python3.12 -m venv .venv && . .venv/bin/activate
pip install -r step3_analysis/requirements.txt python-docx==1.2.0 python-pptx==1.0.2

cd step2_preprocessing && bash run_all.sh   # downloads ~0.9 GB and checks SHA-256; then about 25 min
cd ../step3_analysis   && bash run_all.sh   # about 1.5 h on 4 cores
cd ../step4_slides     && python scripts/build_deck.py
cd ../step5_report     && python scripts/report_figures.py && python scripts/build_report.py
```

The receptor-label rule is set in `step3_analysis/config.json` (or the environment variables `BT3041_HER2_RULE` and
`BT3041_HER2_FILTER`). The report's numbers, all figures and most numbers on the slides are read from the result
files, so a rerun updates them. The exceptions are a few numbers typed into `step4_slides/scripts/build_deck.py` and the
group sizes on the slides and in Table 1, which come from step 2's default rule (v1, no filter). All randomness is
seeded (base seed 42; a few steps use fixed offsets of it).

Build outputs (`step4_slides/previews/`, `step5_report/supplementary/`) and the large data folders
(`step2_preprocessing/data/raw/`, `data/processed/`) are not committed; the scripts recreate them. The report scripts
read `step2_preprocessing/data/processed/`, so run step 2 once before rebuilding the report.

## Data

The large inputs are not stored here; `step2_preprocessing/scripts/00_download.sh` fetches them and checks them
against `step2_preprocessing/data/SHA256SUMS.freeze`. Small inputs that change over time (clinical tables, HGNC gene
names, gene sets) are frozen in `step2_preprocessing/data/raw_small/` so the run is exact.

| Source | Used for | Terms |
|---|---|---|
| TCGA-BRCA via UCSC Xena (TCGA hub, GDC hub, PanCanAtlas) | RNA-seq, clinical labels, RPPA protein, survival (TCGA-CDR) | Open-access TCGA data |
| METABRIC (Curtis et al. 2012; Pereira et al. 2016; Rueda et al. 2019) via the cBioPortal datahub (`brca_metabric`) | External validation cohort | Open Database License (ODbL) v1.0, http://opendatacommons.org/licenses/odbl/1.0/. The METABRIC files here and the results derived from them are made available under the ODbL. |
| MSigDB Hallmark v2024.1 (Liberzon et al. 2015; © Broad Institute, MIT and Regents of the University of California; gsea-msigdb.org), unmodified copy | Pathway enrichment | CC BY 4.0, https://creativecommons.org/licenses/by/4.0/ |
| HGNC complete set (genenames.org) and GENCODE v36 annotation | Gene symbol mapping | Freely available without restrictions |
| PAM50 centroids (Parker et al. 2009), from github.com/perou-lab/PAM50_classifier | Comparator only | No licence stated by the source; published centroids included unmodified for reproducibility |

Full source URLs are in `step2_preprocessing/scripts/00_download.sh`; `step2_preprocessing/README.md` describes each
source. Licences, citations and the METABRIC (ODbL) notice for every included data file are in
`step2_preprocessing/data/DATA_LICENSES.md`.

The results shown here are in whole or part based upon data generated by the TCGA Research Network:
https://www.cancer.gov/tcga.

## Use of AI

This project was done with AI assistance for planning, code, review and drafting. The AI-use
statement is in the report. The AI prompts and chat logs are submitted with the report as supplementary material and
are not part of this public repository.
