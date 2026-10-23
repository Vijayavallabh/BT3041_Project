import sys, time, numpy as np, pandas as pd
counts_gz, probemap, subtype_gz, gmt, surv = sys.argv[1:6]
t0 = time.time()
X = pd.read_csv(counts_gz, sep='\t', index_col=0)
print('star_counts shape (genes x samples):', X.shape, f'{time.time()-t0:.0f}s')
C = np.round(np.power(2.0, X.values) - 1).astype(np.int64)
frac_nonint = np.mean(np.abs((np.power(2.0, X.values) - 1) - C) > 1e-3)
print('fraction non-integer after 2^x-1:', round(float(frac_nonint), 5), 'min count', C.min())
C = pd.DataFrame(C, index=X.index, columns=X.columns)
pm = pd.read_csv(probemap, sep='\t', index_col=0)
sub = pd.read_csv(subtype_gz, sep='\t', index_col=0)
sub = sub[sub['Subtype_Selected'].astype(str).str.startswith('BRCA')]
meta = pd.DataFrame({'sample': C.columns})
meta['s15'] = meta['sample'].str[:15]; meta['type'] = meta['sample'].str[13:15]
meta['pam50'] = meta['s15'].map(sub['Subtype_mRNA'])
tum = meta[meta['type'] == '01']
print('primary tumours with PanCan PAM50 label:', tum['pam50'].value_counts(dropna=False).to_dict())
rng = np.random.default_rng(0)
pick = pd.concat([tum[tum.pam50 == g].sample(30, random_state=0) for g in ['Basal', 'LumA']])
cm = C[pick['sample']].T
cm = cm.loc[:, (cm >= 10).sum(0) >= 10]
cm = cm.loc[:, ~cm.columns.duplicated()]
print('DE subset:', cm.shape)
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
md = pd.DataFrame({'condition': pick['pam50'].values}, index=pick['sample'].values)
t1 = time.time()
dds = DeseqDataSet(counts=cm, metadata=md, design='~condition', quiet=True)
dds.deseq2()
st = DeseqStats(dds, contrast=['condition', 'Basal', 'LumA'], quiet=True); st.summary()
res = st.results_df.copy()
res['symbol'] = pm.reindex(res.index)['gene'].values
print(f'PyDESeq2 done in {time.time()-t1:.0f}s; padj<0.05:', int((res.padj < 0.05).sum()))
for g in ['ESR1', 'FOXA1', 'KRT5', 'FOXC1', 'PGR', 'MKI67']:
    r = res[res.symbol == g]
    if len(r): print(f'  {g}: log2FC(Basal/LumA)={r.log2FoldChange.iloc[0]:.2f} padj={r.padj.iloc[0]:.2e}')
import gseapy as gp
rnk = res.dropna(subset=['stat']).groupby('symbol')['stat'].mean().sort_values(ascending=False)
t2 = time.time()
pre = gp.prerank(rnk=rnk, gene_sets=gmt, permutation_num=200, seed=0, threads=4, outdir=None, verbose=False, min_size=10, max_size=500)
r2 = pre.res2d.sort_values('NES')
print(f'gseapy prerank (offline, local GMT) {time.time()-t2:.0f}s, sets tested:', len(r2))
print(r2[['Term', 'NES', 'FDR q-val']].head(3).to_string(index=False))
print(r2[['Term', 'NES', 'FDR q-val']].tail(3).to_string(index=False))
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
s = pd.read_csv(surv, sep='\t', index_col=0)
lab = sub['Subtype_mRNA']
s = s[s.index.str[13:15] == '01'].join(lab, how='inner')
a = s[s.Subtype_mRNA == 'Basal'].dropna(subset=['OS', 'OS.time']); b = s[s.Subtype_mRNA == 'LumA'].dropna(subset=['OS', 'OS.time'])
lr = logrank_test(a['OS.time'], b['OS.time'], a['OS'], b['OS'])
print(f'lifelines logrank OS Basal (n={len(a)}) vs LumA (n={len(b)}): p={lr.p_value:.3g}')
lr = logrank_test(a['PFI.time'], b['PFI.time'], a['PFI'], b['PFI'])
print(f'lifelines logrank PFI Basal vs LumA: p={lr.p_value:.3g}')
