import sys, csv, gzip, collections
clin, expr = sys.argv[1], sys.argv[2]
with gzip.open(expr,'rt') as f:
    hdr = f.readline().rstrip('\n').split('\t')[1:]
rna = set(hdr)
rows = list(csv.DictReader(open(clin), delimiter='\t'))
print("clin rows", len(rows), "rna cols", len(rna))
cols = ['ER_Status_nature2012','PR_Status_nature2012','HER2_Final_Status_nature2012','PAM50Call_RNAseq','breast_carcinoma_estrogen_receptor_status','breast_carcinoma_progesterone_receptor_status','lab_proc_her2_neu_immunohistochemistry_receptor_status','HER2_Final_Status_nature2012']
print([c for c in rows[0].keys() if 'receptor' in c.lower() or 'her2' in c.lower() or 'nature2012' in c.lower()][:40])
tum = [r for r in rows if r['sampleID'] in rna and r['sample_type']=='Primary Tumor']
print("primary tumour with RNA", len(tum))
for c in cols:
    if c in rows[0]:
        print(c, collections.Counter(r[c] for r in tum).most_common(8))
print("PAM50 x sample_type (RNA)", collections.Counter((r['PAM50Call_RNAseq'], r['sample_type']) for r in rows if r['sampleID'] in rna).most_common(20))
# TNBC by EMR fields
er='breast_carcinoma_estrogen_receptor_status'; pr='breast_carcinoma_progesterone_receptor_status'; h='lab_proc_her2_neu_immunohistochemistry_receptor_status'
if er in rows[0]:
    print("ER EMR x PAM50", collections.Counter((r[er], r['PAM50Call_RNAseq']) for r in tum).most_common(20))
