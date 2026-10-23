# Changelog

Every change to `data/processed/` after the data freeze gets an entry here, with the reason and the files affected.
Checksums of the processed files for each version are in `data/processed/SHA256SUMS` (written by `scripts/07_data_dictionary.py`).

## data-freeze-v1 signed off (23 Oct 2026)

The HER2 rule v1, the 9 exclusions and the supervised groups (591 / 160 / 157) were reviewed against ASCO/CAP 2018 HER2
guidance and kept unchanged. `data/processed/` is not modified; checksums below still apply. Reasons:
`../step3_analysis/README.md`, "Decisions and deviations" item 0.

## data-freeze-v1 (proposed 20 Oct 2026)

First complete version of the processed data, built by `bash run_all.sh` from the raw files listed in
`data/SHA256SUMS.freeze`.

Changes from the plan's preliminary numbers (all explained in README, "Decisions and deviations"):

- HER2 label rule v1 uses the IHC score as well as the status text (plan rule kept as `HER2_v0`). IHC 3+ with a
  negative FISH, 2+ without FISH, and status/score contradictions without FISH are missing. Calls that rest on the
  status text alone are flagged `HER2_confidence == "low"` (239 supervised tumours: 30 of the 160 HER2+, 160 HR+/HER2−,
  46 TNBC, 3 incomplete; README decision 1 gives the sensitivity analyses for step 3).
- 9 tumours excluded from all analyses (5 sample-identity failures, 5 TCGA-CDR redacted, 1 in both). The adjacent
  normal of the redacted case is dropped as well, from both normal sets (HiSeqV2 113 normals, STAR 112).
- The 12 male patients are left out of the supervised set (2 tumours with blank sex are female by XIST/RPS4Y1 expression and stay in).
- STAR gene filter uses the size of the smallest clinical group (now TNBC, n = 157) instead of the fixed 160.
- CPTAC 2020 fallback cohort descoped: METABRIC downloaded and passed every check, so CPTAC is not needed.

Processed-file checksums at this version:

```
51d8018f2233a9aba1292d9c7d02bcd5f3117151a81dd0541d33fed9288c58ca  data/processed/erbb2_amplicon_genes.tsv
8017c52d911311f805544e4c8aa420ad578076fa624003069694f0c543797306  data/processed/gencode_v36_annotation.tsv
a81ee878d0adc9a5ebf9d3e74e268b28ecc04248a0f0297ba329f347fd533a79  data/processed/gene_map.tsv
214202d2f98cb6874e2f3a38ceef6618636be4e11d847d87d86d6d95e045f606  data/processed/h.all.v2024.1.hgnc.gmt
cce6fb197f3a0b28ecab0ab6c4cfd8471a711c5dd35f44c1f4770e528a1a451c  data/processed/metabric_clusters.tsv
d014d135b2b009d69b3c0293cd5f6ea36c83e5cc9b90335fe301b550587fab49  data/processed/metabric_labels.tsv
36fb01afb3b1927cba11933e0ab2f70fa7b92d1bd5aae68c97874a2973632645  data/processed/metabric_log2.parquet
60e0e57aa4c4828b4e3702d2e1f099c5fb8bc3e0d6cbb5473306459428d74227  data/processed/pam50_centroids_hgnc.tsv
bfcd08cf792b6ac5c74cb7418de478b7db0b600e31bbe661abb079a80f7a457f  data/processed/rppa_antibody_map.tsv
93dca89a782e7e0df58dee23f43daab1123580502983ffb3f53b14bcee41a55a  data/processed/sex_specific_genes.txt
5658f8be91d09c33c9d9c9aba8cfd3f1d1a0d169b579747ec7f3fa987bd56506  data/processed/tcga_clusters.tsv
8a8286b5a2e85ff3a1786c840c2179b8c324192fc2dd406ba5fa4cef342c565a  data/processed/tcga_hiseq_normal_log2.parquet
71f2a458d5e680acc20c3ad7acc80dfe6009d7312b91c853266a6ad6b0e545a3  data/processed/tcga_hiseq_tumour_log2.parquet
08dd3a15617aa41a3fc4136011f4840fad9e5d9b04dce6ff81fd5ec820cf534f  data/processed/tcga_kmeans_k2_centroids_z.tsv
f4d656f13c1fd38d9322dc1d0bf9faee05110cc906f874641dc20d02f97d5912  data/processed/tcga_kmeans_k4_centroids_z.tsv
fb2d0f371f7fb601a90d1238102ddd328acd2c868a13a8e2242ff149af60f753  data/processed/tcga_kmeans_k5_centroids_z.tsv
bda3f8aaa6ea9d49b24b7568da5cd46ff165014cf81bebb0562e43aa805d700f  data/processed/tcga_labels.tsv
9c973508bfab4b43152b5daa3672a526abecf790179238524c29ee6eaa5d31a8  data/processed/tcga_rppa.parquet
e4106eff4ae6ed135e1bbf45ac4befaee6a03cb4fb02d20d736489723ddf24e8  data/processed/tcga_star_counts_all.parquet
b03218c8a1294e5e1d2b7cda5e920745c7834c78463722bac9c1f5884fb68ec0  data/processed/tcga_star_counts_tumour_filtered.parquet
```
