"""Download cBioPortal datahub study files that are stored in Git LFS.
Usage: python3 -I fetch_cbio_lfs.py <study_id> <out_dir> <file1> [file2 ...]
e.g.   python3 -I fetch_cbio_lfs.py brca_metabric ./metabric data_clinical_patient.txt data_clinical_sample.txt data_mrna_illumina_microarray.txt
Stdlib only. Reads the LFS pointer from raw.githubusercontent.com, asks the datahub LFS
server (from the repo's .lfsconfig) for a pre-signed S3 URL, streams the file to <file>.part, checks
size and sha256, and only then renames it into place. Exits with code 1 on any mismatch."""
import sys, os, json, hashlib, urllib.request
RAW = "https://raw.githubusercontent.com/cBioPortal/datahub/master/public/{study}/{f}"
LFS = "https://nsssw8k94d.execute-api.us-east-1.amazonaws.com/objects/batch"
study, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
failed = []
for f in sys.argv[3:]:
    dest = os.path.join(out, f)
    ptr = urllib.request.urlopen(RAW.format(study=study, f=f), timeout=60).read().decode()
    if not ptr.startswith("version https://git-lfs"):
        open(dest, "w").write(ptr); print(f, "plain file", len(ptr)); continue
    kv = dict(l.split(" ", 1) for l in ptr.strip().splitlines())
    oid, size = kv["oid"].split(":", 1)[1], int(kv["size"])
    body = json.dumps({"operation": "download", "transfers": ["basic"], "objects": [{"oid": oid, "size": size}]}).encode()
    req = urllib.request.Request(LFS, data=body, method="POST", headers={
        "Accept": "application/vnd.git-lfs+json", "Content-Type": "application/vnd.git-lfs+json"})
    href = json.load(urllib.request.urlopen(req, timeout=60))["objects"][0]["actions"]["download"]["href"]
    h, n = hashlib.sha256(), 0
    with urllib.request.urlopen(href, timeout=600) as r, open(dest + ".part", "wb") as fh:
        while chunk := r.read(1 << 20):
            fh.write(chunk); h.update(chunk); n += len(chunk)
    if h.hexdigest() == oid and n == size:
        os.replace(dest + ".part", dest)
        print(f, n, "bytes, sha256 OK")
    else:
        print(f, n, "bytes, CHECKSUM MISMATCH (left as .part)")
        failed.append(f)
sys.exit(1 if failed else 0)
