# Step 4a: in-class presentation (3–6 Nov 2026)

| File | What it is |
|---|---|
| `BT3041_talk.pptx` | The 10-slide 16:9 deck, editable (imports into Google Slides via File → Import slides). Speaker notes on every slide start with `[presenter · time]`. |
| `BT3041_talk.pdf` | PDF export of the same deck (the safe copy to present from). |
| `SPEAKER_SCRIPT.md` | All speaker notes with the timing table, for rehearsal. 7:45 in total, about 130 words a minute. |
| `QA_prep.md` | 21 likely questions with short answers and a suggested owner. |
| `previews/` | Low-resolution PNG of each slide (not committed; `scripts/build_deck.py` creates it). |
| `scripts/build_deck.py` | Rebuilds everything above from step 2 and step 3 outputs. |

## Storyline (guideline check: ≤ 10 slides, intro + data preparation, ≥ 40% methods/results)

| # | Slide | Presenter | Time |
|---|---|---|---|
| 1 | Title and question, group sizes | Vijayavallabh (M1) | 0:20 |
| 2 | Background, 3 research questions, validation design (why not PAM50) | Vijayavallabh (M1) | 0:45 |
| 3 | Data and preparation (sample flow) | Vijayavallabh (M1) | 0:45 |
| 4 | Exploration: UMAP + marker violins | Vedant (M2) | 0:40 |
| 5 | Clustering vs labels (k = 2 stable, ARI) | Vedant (M2) | 0:45 |
| 6 | Differential expression (2,741 DEGs, robustness) | Harshit (M3) | 0:50 |
| 7 | Pathways (GSEA with correlation-aware checks), TNBC vs Basal | Harshit (M3) | 0:50 |
| 8 | Classification: models match, not beat, ESR1/ERBB2; ablation | Navya (M4) | 1:05 |
| 9 | Validation: RPPA, ER error analysis, METABRIC transfer | Anitha (M5) | 1:00 |
| 10 | Conclusions, limitations, done / for the final report | Anitha (M5) | 0:45 |

Methods and results (slides 4–9) are 67% of the talk. Presenters follow the M1–M5 workstreams in the team plan
(contiguous blocks, 4 handovers). Survival is not in the 10 slides; if the instructor allows a backup slide, use
`../step3_analysis/figures/slideS1_survival.png` (Q&A answer 18 covers it otherwise).

## Before the talk

- Team 16's names are on slide 1 and in the notes. Roll numbers on slide 1 come from `../roll_numbers.json`, which
  is kept out of the public repository; without it the slide lists names only. Presenters follow the plan's workstreams in the order the
  names were given (M1 Vijayavallabh, M2 Vedant, M3 Harshit, M4 Navya, M5 Anitha); edit `TEAM` in
  `scripts/build_deck.py` to reassign.
- Do two timed dry runs (8:00 hard limit); cut words, not slides.
- Edits made by hand in the .pptx are lost if the deck is rebuilt; make wording changes in `scripts/build_deck.py`
  until the last rebuild (the labels were signed off on 23 Oct, so the numbers are final).

## If the labels ever change (HER2 rule, exclusions, groups)

1. Rerun step 2 if its labels change, then step 3: set `../step3_analysis/config.json` and `bash run_all.sh`
   (figures 6–9 and every report number regenerate).
2. `python scripts/build_deck.py` (needs `pip install python-pptx` and pandas; LibreOffice for the PDF). Slides 3–9 are
   the regenerated figures; every number on slides 1, 2 and 10 and in the notes is read from the result JSON files.
3. Recheck `QA_prep.md`: it is hand-written, and the numbers that move with the labels are marked †.
