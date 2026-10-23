#!/usr/bin/env bash
# Reproduce step 2 end to end: about 25 min on a 4-core machine (03_cluster.py is ~15 of them); METABRIC download is 689 MB.
# Usage: bash run_all.sh            (uses the python on PATH; override with PY=/path/to/python bash run_all.sh)
set -euo pipefail
cd "$(dirname "$0")"
export PY=${PY:-python}
mkdir -p logs
bash scripts/00_download.sh 2>&1 | tee logs/00_download.log
for s in 01_preprocess 02_explore 03_cluster 03b_cluster_figures 04_metabric_explore 05_sample_flow 06_slide_figures 07_data_dictionary; do
  echo "== $s"
  "$PY" "scripts/$s.py" 2>&1 | tee "logs/$s.log"
done
