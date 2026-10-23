import sys, csv, gzip, collections
clin, expr = sys.argv[1], sys.argv[2]
with gzip.open(expr,'rt') as f:
    cols = f.readline().rstrip('\n').split('\t')[1:]
rna=set(cols)
allr = list(csv.DictReader(open(clin), delimiter='\t'))
rows = [r for r in allr if r['sampleID'] in rna and r['sample_type']=='Primary Tumor']
ihc='lab_proc_her2_neu_immunohistochemistry_receptor_status'; fish='lab_procedure_her2_neu_in_situ_hybrid_outcome_type'
er='breast_carcinoma_estrogen_receptor_status'; pr='breast_carcinoma_progesterone_receptor_status'
def h2(r):
    i,f=r[ihc],r[fish]
    if f=='Positive' or (i=='Positive' and f!='Negative'): return 'Pos'
    if f=='Negative' or (i=='Negative'): return 'Neg'
    return 'NA'
g=collections.Counter()
for r in rows:
    H=h2(r); E,P=r[er],r[pr]
    if H=='NA' or E not in('Positive','Negative') or P not in('Positive','Negative'): g['incomplete']+=1; continue
    if H=='Pos': g['HER2+']+=1
    elif E=='Positive' or P=='Positive': g['HR+/HER2-']+=1
    else: g['TNBC']+=1
print(g)
pat=collections.defaultdict(set)
for c in cols: pat[c[:12]].add(c[13:15])
print("pairs 01+11:", sum(1 for p,s in pat.items() if '01' in s and '11' in s), "patients", len(pat), "dup 01?", len(cols)-len(set(cols)))
print("normals in RNA", sum(1 for c in cols if c[13:15]=='11'))
k='SigClust_Unsupervised_mRNA_nature2012'
print(k, sum(1 for r in rows if r[k].strip()), collections.Counter(r[k] for r in rows).most_common(6))
print('PAM50_mRNA_nature2012', collections.Counter(r['PAM50_mRNA_nature2012'] for r in rows).most_common(7))
