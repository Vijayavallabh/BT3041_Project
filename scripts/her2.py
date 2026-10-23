import sys, csv, gzip, collections
clin, expr = sys.argv[1], sys.argv[2]
with gzip.open(expr,'rt') as f:
    rna = set(f.readline().rstrip('\n').split('\t')[1:])
rows = [r for r in csv.DictReader(open(clin), delimiter='\t') if r['sampleID'] in rna and r['sample_type']=='Primary Tumor']
ihc='lab_proc_her2_neu_immunohistochemistry_receptor_status'; fish='lab_procedure_her2_neu_in_situ_hybrid_outcome_type'
print("FISH", collections.Counter(r[fish] for r in rows).most_common())
print("IHC x FISH", collections.Counter((r[ihc], r[fish]) for r in rows).most_common(15))
def h2(r):
    i,f=r[ihc],r[fish]
    if f=='Positive' or (i=='Positive' and f!='Negative'): return 'Pos'
    if f=='Negative' or (i=='Negative'): return 'Neg'
    return 'NA'
er='breast_carcinoma_estrogen_receptor_status'; pr='breast_carcinoma_progesterone_receptor_status'
c=collections.Counter()
for r in rows:
    H=h2(r); c['HER2_'+H]+=1
    if r[er]=='Negative' and r[pr]=='Negative' and H=='Neg': c['TNBC']+=1
    if r[er] in('Positive','Negative') and r[pr] in('Positive','Negative') and H!='NA': c['complete_ER_PR_HER2']+=1
print(c)
print("TNBC x PAM50", collections.Counter(r['PAM50Call_RNAseq'] for r in rows if r[er]=='Negative' and r[pr]=='Negative' and h2(r)=='Neg'))
print("HER2pos x PAM50", collections.Counter(r['PAM50Call_RNAseq'] for r in rows if h2(r)=='Pos'))
print("histology", collections.Counter(r.get('histological_type','?') for r in rows).most_common(5))
