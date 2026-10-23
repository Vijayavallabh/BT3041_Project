"""Assemble the BT3041 supplementary package (code, results, logs, AI prompts and chat histories).

Usage: python scripts/make_supplementary.py
Writes step5_report/supplementary/bt3041_supplementary/ and step5_report/supplementary/bt3041_supplementary.zip.
Large data (step2_preprocessing/data/processed and raw, ~600 MB) is left out: run_all.sh re-downloads it and
checks it against data/SHA256SUMS.freeze. The small frozen copies in data/raw_small (23 MB; clinical tables, HGNC,
gene sets) are kept, because 00_download.sh seeds from them and some sources are re-released over time. Re-run this script after any rerun or rebuild.
"""
import shutil
import zipfile
from pathlib import Path

STEP5 = Path(__file__).resolve().parents[1]
PLAN = STEP5.parent
OUT = STEP5 / "supplementary" / "bt3041_supplementary"
ZIP = OUT.with_suffix(".zip")

SKIP_DIRS = {"__pycache__", ".ipynb_checkpoints", "previews", "supplementary", ".venv"}
SKIP_SUFFIX = {".pyc", ".lock"}

# (source relative to bt3041-plan, destination relative to the package)
ITEMS = [
    ("project_plan.md", "project_plan.md"),
    ("roll_numbers.json", "roll_numbers.json"),  # private, read by the report and slide builds; absent in the public repo
    ("scripts", "scripts"),
    ("step2_preprocessing", "step2_preprocessing"),
    ("step3_analysis", "step3_analysis"),
    ("step4_slides", "step4_slides"),
    ("step5_report", "step5_report"),
]
# Inside step 2 only the checksum list and the small frozen raw copies of the data folder are kept.
STEP2_DATA_KEEP = {"SHA256SUMS.freeze", "raw_small", "DATA_LICENSES.md"}


def keep(src: Path) -> bool:
    if any(p in SKIP_DIRS for p in src.parts):
        return False
    if src.suffix in SKIP_SUFFIX:
        return False
    rel = src.relative_to(PLAN).parts
    if rel[:2] == ("step2_preprocessing", "data") and len(rel) > 2 and rel[2] not in STEP2_DATA_KEEP:
        return False
    return True


def copy_item(src: Path, dst: Path) -> int:
    n = 0
    files = [src] if src.is_file() else sorted(p for p in src.rglob("*") if p.is_file())
    for f in files:
        if not keep(f):
            continue
        target = dst if src.is_file() else dst / f.relative_to(src)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(f, target)
        n += 1
    return n


README = """# BT3041 term project: supplementary material

Team project for BT3041 (Analysis and Interpretation of Biological Data, IIT Madras), 2026.
Question: what TCGA breast tumour gene expression reveals about the clinical receptor groups (ER/PR/HER2),
validated on TCGA RPPA protein data and the independent METABRIC cohort.
Code repository (public, without the AI chat logs): https://github.com/Vijayavallabh/BT3041_Project

## Where to look

| Folder | What it holds |
|---|---|
| `AI_PROMPTS_AND_CHATS/` | **Start here for AI use.** `AI_CHAT_HISTORY.md` (every prompt the team typed, verbatim, and every AI reply) and an index of the detailed per-step AI logs. |
| `project_plan.md`, `scripts/` | Team plan (workstreams, timeline) and small helper scripts used while planning. |
| `step2_preprocessing/` | Download, QC, labels, normalisation, PCA/UMAP, clustering, METABRIC projection. `README.md` is the data dictionary and lists every decision. `ai_logs/` has the AI log, the planning record and three independent review rounds. |
| `step3_analysis/` | Differential expression, enrichment, classifiers, METABRIC transfer, RPPA, survival. `config.json` holds the label choices. |
| `step4_slides/` | Deck build script, slides (.pptx/.pdf), speaker script, Q&A prep. |
| `step5_report/` | Report build scripts, the report (.docx/.pdf), `report_numbers.md` (every number in the report and its source). |

Each step folder has `scripts/`, `results/` or outputs, `logs/` (console output of the last run), a README and,
where relevant, a CHANGELOG.

## Reproduce everything

Python 3.12. Install `step2_preprocessing/requirements.txt` and
`step3_analysis/requirements.txt` (plus `python-docx`, `python-pptx`; LibreOffice for PDFs).

```
cd step2_preprocessing && bash run_all.sh      # downloads ~0.7 GB, checks SHA-256, ~25 min
cd ../step3_analysis   && bash run_all.sh      # ~1 h on 4 cores (nested CV is the slow part)
cd ../step4_slides          && python scripts/build_deck.py
cd ../step5_report                && python scripts/report_figures.py && python scripts/build_report.py
```

The large raw and processed data are not included (about 600 MB; TCGA and METABRIC are public). `run_all.sh`
downloads them and verifies them against `step2_preprocessing/data/SHA256SUMS.freeze`; the small frozen inputs it
seeds from (clinical tables, HGNC, gene sets, 23 MB) are included in `step2_preprocessing/data/raw_small/`.

## Data licences

TCGA data via UCSC Xena (open access). METABRIC via cBioPortal; its terms apply. See the step 2 README for sources.
"""


def main():
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True)
    counts = {}
    for src, dst in ITEMS:
        counts[dst] = copy_item(PLAN / src, OUT / dst)
    # AI prompts and chat histories, gathered in one place
    ai = OUT / "AI_PROMPTS_AND_CHATS"
    ai.mkdir()
    chat = STEP5 / "ai_logs" / "AI_CHAT_HISTORY.md"  # not in the public repo
    if chat.exists():
        shutil.copy2(chat, ai / "AI_CHAT_HISTORY.md")
    (ai / "INDEX.md").write_text(
        "# AI prompts and chat histories\n\n"
        "- `AI_CHAT_HISTORY.md`: all prompts typed by the team (verbatim) and all AI replies, in order.\n"
        "- Planning: `../step2_preprocessing/ai_logs/planning/`\n"
        "- Step 2: `../step2_preprocessing/ai_logs/step2_preprocessing_exploration.md`\n"
        "- Step 2 review rounds (reviewer prompts + verified findings): `../step2_preprocessing/ai_logs/step2_review/`\n"
        "- Step 3: `../step3_analysis/ai_logs/step3_analysis_validation.md`\n"
        "- Step 4: `../step4_slides/ai_logs/step4_slides.md`\n"
        "- Step 5: `../step5_report/ai_logs/step5_report.md`\n\n"
        "TEAM TO DO before upload: add exports of any other AI chats the team used (one file per member is fine).\n"
    )
    (OUT / "README.md").write_text(README)
    if ZIP.exists():
        ZIP.unlink()
    with zipfile.ZipFile(ZIP, "w", zipfile.ZIP_DEFLATED) as z:
        for f in sorted(OUT.rglob("*")):
            if f.is_file():
                z.write(f, Path(OUT.name) / f.relative_to(OUT))
    for k, v in counts.items():
        print(f"{v:5d} files  {k}")
    print(f"zip: {ZIP} ({ZIP.stat().st_size / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
