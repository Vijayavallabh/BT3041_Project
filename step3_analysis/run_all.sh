#!/usr/bin/env bash
# Step 3: core analysis and validation. Needs step 2's data/processed/ next to this folder (../step2_preprocessing).
# Label inputs: config.json (or BT3041_HER2_RULE / BT3041_HER2_FILTER). About 1.5 h on 4 cores.
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-python}
mkdir -p logs
run() { echo "== $1"; $PY "scripts/$1" > "logs/${1%.py}.log" 2>&1; }
run 10_de.py                 # PyDESeq2, 7 runs
run 11_enrichment.py         # GSEA prerank, ORA, CAMERA, label-permutation GSEA, PAM50 agreement
run 12_classify.py           # nested CV, 4 models x ablations, controls, sensitivity, PAM50 benchmark
run 13_freeze_transfer.py    # freeze METABRIC models on TCGA (writes the metabric-freeze CHANGELOG entry)
run 14_metabric.py           # transfer + DE replication
run 15_rppa.py               # protein validation + ER error analysis
run 15b_classifier_features.py  # stable classifier features + Hallmark ORA
run 16_survival.py           # TCGA and METABRIC survival
run 17_figures_de.py
run 18_figures_validation.py
run 19_literature_table.py
echo "done"
