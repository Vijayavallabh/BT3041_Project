# Step 5: final report

Course limits: at most 8 pages, Arial 11, single spacing. This is the pre-submission draft. The submitted copy also
carries the team's roll numbers and the link to the AI-log supplement, which the build reads from files kept out of
this repository (`../roll_numbers.json` and `supplement_link.txt`); without them it lists names only and says the
AI logs are shared privately.

## Files

| File | What it is |
|---|---|
| `BT3041_final_report.docx` | Editable report (A4, 2 cm margins, Arial 11, single spacing). |
| `BT3041_final_report.pdf` | PDF of the same, 8 pages including references (the build checks the page count). |
| `report_numbers.md` | The key numbers used in the text, generated from the step 2 and step 3 result files. |
| `figures/report_fig1…5` | The five report figures (.png at 300 dpi and .pdf). |
| `scripts/report_figures.py` | Makes the figures from the step 2/3 results, reusing the step 3 plotting code. |
| `scripts/build_report.py` | Writes the .docx from the result files, converts it to PDF and fails if it is over 8 pages. |
| `scripts/make_supplementary.py` | Packages code, results and logs for the course submission (output is not committed). |

## Report structure

Abstract · 1 Introduction (three questions) · 2 Data and preprocessing (Table 1) · 3 Methods · 4 Results:
4.1 exploration and clustering, 4.2 differential expression and pathways, 4.3 classification (whole-transcriptome
models match, but do not beat, ESR1/ERBB2), 4.4 validation on RPPA, METABRIC and published studies (Table 2),
4.5 survival · 5 Conclusion (interpretation, limitations, future work) · Data and code availability · Use of AI ·
26 references.

## Rebuilding

```
# Python 3.12 with step3_analysis/requirements.txt + python-docx + python-pptx; LibreOffice (soffice) and pdfinfo
python scripts/report_figures.py
python scripts/build_report.py          # prints the page count; exits with an error above 8 pages
```

Both scripts also read `../step2_preprocessing/data/processed/` (labels, clusters, METABRIC labels, ERBB2 amplicon
genes), which is not in the repository: run `../step2_preprocessing/run_all.sh` once first. Every number in the text
comes from the result files, so after a change to the label rule, rerun step 2 and step 3 and then these two scripts.
Edit the text in `build_report.py`, not in the .docx, because a rebuild overwrites the .docx.

The PDF was made with Liberation Sans, which has the same metrics as Arial, so the layout is identical.
