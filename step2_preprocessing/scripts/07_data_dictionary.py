"""Step 7: data dictionary and checksums for every processed file (runs in under a minute).

Writes data/processed/DATA_DICTIONARY.tsv (file, row key, shape, orientation, units, consumers) and
data/processed/SHA256SUMS (the data-freeze checksums recorded in CHANGELOG.md).
Run from the project root:  python scripts/07_data_dictionary.py   (after all other steps)
"""
import hashlib, sys
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.common import PROC, ROOT

# file -> (orientation and units, consumers in the team plan); {n} is filled with the file's row count
DESC = {
    "tcga_labels.tsv": ("rows = all {n:,} primary tumours (15-character barcode, e.g. TCGA-3C-AAAU-01), cols = labels and QC flags; "
                        "drop rows with exclude == True (or use src.common.drop_excluded); supervised work also needs supervised_set == True",
                        "everyone"),
    "tcga_hiseq_tumour_log2.parquet": ("rows = analysed tumours (excluded removed), cols = HGNC symbols; log2(RSEM upper-quartile + 1), "
                                       "unscaled; only zero-variance genes removed", "M2 exploration/clustering, M4 classifiers"),
    "tcga_hiseq_normal_log2.parquet": ("rows = {n} adjacent normals (-11 barcodes; normals of TCGA-CDR redacted patients removed), cols = same genes as the tumour matrix; log2(RSEM + 1)",
                                       "M2 (tumour-normal PCA), M3 (optional)"),
    "tcga_star_counts_all.parquet": ("rows = all STAR primary tumours after vial dedup (-01, excluded tumours included; filter with "
                                     "tcga_labels.exclude) and adjacent normals (-11; the normal of the TCGA-CDR redacted case removed); "
                                     "cols = all versioned Ensembl IDs (GENCODE v36), no gene filter; integer counts",
                                     "M3 (custom filters, tumour-normal checks)"),
    "tcga_star_counts_tumour_filtered.parquet": ("rows = analysed tumours with STAR (15-character barcode), cols = versioned Ensembl IDs "
                                                 "passing >=10 counts in >= smallest-group-size tumours; integer counts", "M3 PyDESeq2"),
    "tcga_rppa.parquet": ("rows = primary tumours with RPPA (15-character barcode; filter with tcga_labels.exclude), cols = antibody "
                          "names; replicate-based normalised (RBN) RPPA values",
                          "M5 protein validation"),
    "metabric_log2.parquet": ("rows = {n:,} METABRIC tumours (MB-xxxx), cols = HGNC symbols; Illumina HT-12 log2 intensity", "M2, M3, M4"),
    "metabric_labels.tsv": ("rows = METABRIC tumours, cols = ER_IHC, HER2_SNP6, contrast_group, survival, clinical", "M2, M3, M4"),
    "metabric_clusters.tsv": ("rows = METABRIC tumours, cols = nearest TCGA k-means centroid (Euclidean) per k; Pearson columns are "
                              "a labelled sensitivity only", "M2 survival, M4"),
    "tcga_clusters.tsv": ("rows = analysed tumours, cols = cluster IDs (k-means, Ward, average-linkage, consensus; sensitivity runs); "
                          "C0 = largest cluster", "M2, M3 (optional), M5"),
    "gene_map.tsv": ("one row per input gene name per platform -> HGNC approved symbol and the mapping method", "everyone"),
    "gencode_v36_annotation.tsv": ("rows = versioned Ensembl IDs (GENCODE v36), cols = GENCODE gene name, chrom, chromStart, chromEnd (GRCh38), HGNC "
                                   "approved symbol and mapping method",
                                   "M3 (gene names for DE output)"),
    "erbb2_amplicon_genes.tsv": ("rows = GENCODE v36 genes within 2 Mb of the ERBB2 start (GRCh38), with HGNC symbol, distance and "
                                 "within_1000kb / 1500kb / 2000kb flags",
                                 "M4 ablation (drop the 17q12 amplicon)"),
    "pam50_centroids_hgnc.tsv": ("rows = 50 PAM50 genes (HGNC symbols), cols = subtype centroids", "M2, M4 (PAM50 ablation)"),
    "h.all.v2024.1.hgnc.gmt": ("MSigDB Hallmark gene sets with every gene mapped to HGNC approved symbols", "M3 GSEA/ORA"),
    "rppa_antibody_map.tsv": ("rows = RPPA antibodies, cols = gene(s), isoform_ambiguous, phospho_or_cleaved, total_protein", "M5"),
    "sex_specific_genes.txt": ("one HGNC symbol per line: chrY genes plus XIST and TSIX (left out of unsupervised gene sets)", "M2"),
}


def shape(p):
    """(row key, n_rows, n_cols). n_cols counts every column in the file, the key column included for TSVs."""
    if p.suffix == ".parquet":
        d = pd.read_parquet(p)
        return d.index.name, d.shape[0], d.shape[1]
    if p.suffix == ".tsv":
        d = pd.read_csv(p, sep="\t", low_memory=False)
        key = "dataset + row_id" if p.name == "gene_map.tsv" else d.columns[0]
        return key, d.shape[0], d.shape[1]
    return None, sum(1 for _ in open(p)), None


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        while chunk := fh.read(1 << 20):
            h.update(chunk)
    return h.hexdigest()


rows, sums = [], []
for p in sorted(PROC.iterdir()):
    if p.name in ("DATA_DICTIONARY.tsv", "SHA256SUMS"):
        continue
    if p.name.startswith("tcga_kmeans_k") and p.name.endswith("_centroids_z.tsv"):
        d, c = ("rows = the {n:,} main clustering genes, cols = cluster centroids (C0 = largest) in TCGA z-units; assign other "
                "cohorts by Euclidean nearest centroid after z-scoring within the cohort", "M2, M4 (METABRIC assignment)")
    else:
        d, c = DESC.get(p.name, ("(not described)", ""))
    key, r, k = shape(p)
    rows.append((p.name, key, r, k, d.format(n=r), c))
    sums.append(f"{sha(p)}  data/processed/{p.name}")
dd = pd.DataFrame(rows, columns=["file", "row_key", "n_rows", "n_cols", "orientation_and_units", "consumers"]).astype({"n_cols": "Int64"})
dd.to_csv(PROC / "DATA_DICTIONARY.tsv", sep="\t", index=False)
(PROC / "SHA256SUMS").write_text("\n".join(sums) + "\n")
missing = [f for f in DESC if not (PROC / f).exists()]
assert not missing, f"described but missing: {missing}"
print(dd[["file", "row_key", "n_rows", "n_cols"]].to_string(index=False))
print("done")
