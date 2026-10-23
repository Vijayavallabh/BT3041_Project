#!/usr/bin/env bash
# Step 0: download every raw input, then check it against the data-freeze checksums in data/SHA256SUMS.freeze.
# Usage: bash scripts/00_download.sh [--update]   (run from the project root)
#   Re-runnable: files already present are skipped. --update rewrites the freeze file (only when moving the freeze).
set -euo pipefail
cd "$(dirname "$0")/.."
PY=${PY:-python}
T=data/raw/tcga; R=data/raw/ref; M=data/raw/metabric
mkdir -p "$T" "$R" "$M"
# seed the small freeze copies (clinical tables, references, HGNC set) so the freeze rebuilds exactly; never overwrite
if [ -d data/raw_small ]; then
  find data/raw_small -type f | while read -r f; do
    d="data/raw/${f#data/raw_small/}"
    [ -e "$d" ] || { mkdir -p "$(dirname "$d")"; cp "$f" "$d"; echo "seed $d"; }
  done
fi
get() { [ -s "$2" ] && { echo "have $2"; return; }; echo "get  $2"; curl -fsSL --retry 4 -o "$2.part" "$1" && mv "$2.part" "$2"; }
XENA=https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download
GDC=https://gdc-hub.s3.us-east-1.amazonaws.com/download
get "$XENA/TCGA.BRCA.sampleMap%2FHiSeqV2.gz"            $T/HiSeqV2.gz
get "$XENA/TCGA.BRCA.sampleMap%2FBRCA_clinicalMatrix"   $T/BRCA_clinicalMatrix.tsv
get "$GDC/TCGA-BRCA.star_counts.tsv.gz"                 $T/TCGA-BRCA.star_counts.tsv.gz
get "$GDC/gencode.v36.annotation.gtf.gene.probemap"     $R/gencode.v36.probemap.tsv
get "$XENA/survival%2FBRCA_survival.txt"                $T/BRCA_survival.tsv
get "https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/TCGASubtype.20170308.tsv.gz" $T/TCGASubtype.20170308.tsv.gz
get "$XENA/TCGA.BRCA.sampleMap%2FRPPA_RBN.gz"           $T/RPPA_RBN.gz
# MSigDB Hallmark v2024.1 (public mirror; the official copy at gsea-msigdb.org needs a free login, same content)
get "https://raw.githubusercontent.com/jokergoo/GSEAtopics/master/h.all.v2024.1.Hs.symbols.gmt" $R/h.all.v2024.1.Hs.symbols.gmt
get "https://raw.githubusercontent.com/perou-lab/PAM50_classifier/main/PAM50_R/pam50_centroids.txt" $R/pam50_centroids.txt
# HGNC complete set (official HGNC Google Cloud bucket; re-released monthly, so the freeze copy is kept in the project files)
get "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt" $R/hgnc_complete_set.txt
# METABRIC (cBioPortal datahub, Git LFS; sha256 also checked by the fetch script)
for f in data_clinical_patient.txt data_clinical_sample.txt data_mrna_illumina_microarray.txt; do
  [ -s "$M/$f" ] && { echo "have $M/$f"; continue; }
  "$PY" -I scripts/fetch_cbio_lfs.py brca_metabric "$M" "$f"
done
"$PY" -I scripts/verify_sha256.py data/SHA256SUMS.freeze "$@"
echo done
