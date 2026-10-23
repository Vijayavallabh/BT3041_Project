"""Check downloaded raw files against the data-freeze checksums.
Usage: python3 -I scripts/verify_sha256.py data/SHA256SUMS.freeze [--update | --allow-volatile]
Exit code 1 if any file is missing or differs. A file that is re-released upstream (the HGNC set) also fails, because the
freeze copy ships in data/raw_small; pass --allow-volatile to accept a newer upstream copy with a warning.
--update rewrites the freeze file from the files on disk (use only when deliberately moving the data freeze)."""
import hashlib, sys
from pathlib import Path

# HGNC publishes a new complete set every month at the same URL; a mismatch there changes a few gene mappings only.
VOLATILE = {"data/raw/ref/hgnc_complete_set.txt"}
freeze = Path(sys.argv[1])
rows = [l.split(None, 1) for l in freeze.read_text().splitlines() if l.strip()]


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        while chunk := fh.read(1 << 20):
            h.update(chunk)
    return h.hexdigest()


if "--update" in sys.argv:
    freeze.write_text("".join(f"{sha(p.strip())}  {p.strip()}\n" for _, p in rows))
    print("freeze file rewritten:", freeze)
    sys.exit(0)
bad = 0
for want, p in rows:
    p = p.strip()
    if not Path(p).exists():
        print("MISSING ", p); bad += 1; continue
    got = sha(p)
    if got == want:
        print("ok      ", p)
    elif p in VOLATILE and "--allow-volatile" in sys.argv:
        print("WARNING ", p, "differs from the freeze (monthly upstream release); gene mapping may differ slightly")
    elif p in VOLATILE:
        print("MISMATCH", p, "(monthly upstream release; copy data/raw_small/ref/ to data/raw/ref/ or pass --allow-volatile)"); bad += 1
    else:
        print("MISMATCH", p); bad += 1
sys.exit(1 if bad else 0)
