# Step 3 changelog

## metabric-freeze (2026-10-21T02:04:04Z)

Transfer models frozen on TCGA before any METABRIC label was read (script 13). Shared genes: 18936. Threshold 0.5. Models and SHA-256:

- `ER`: k = 50, C = 0.01, 6 non-zero weights, sha256 `615e0a5f3097f0cf…`
- `ER_ablated`: k = 50, C = 0.1, 20 non-zero weights, sha256 `b9d1983a843094bf…`
- `HER2`: k = 50, C = 0.1, 19 non-zero weights, sha256 `fec973a0eb99d9f6…`
- `HER2_ablated`: k = 1000, C = 0.1, 143 non-zero weights, sha256 `2e728ecb5e97c30e…`
- `ERneg_HER2neg`: k = 50, C = 0.1, 13 non-zero weights, sha256 `3c9af0c5156897e2…`

## slide-figure legibility (2026-10-21, step 4)

`18_figures_validation.py`: on slide 8 the ROC legends use short model names (RF, SVM) and the panels are spaced wider,
so the HER2 legend no longer runs into the AUC panel's labels; on slide 9 the PAM50 legend moved to the empty upper left,
away from the METABRIC panel's "ER−/HER2− vs rest" label. Figure content and numbers are unchanged; only this script
was rerun.
