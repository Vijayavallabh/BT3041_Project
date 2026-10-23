"""Shared paths, label loading and helpers for step 3 (DE, enrichment, classification, validation, survival).

Step 3 reads step 2's processed data and never writes to it. Plotting style, colours and the exclusion rule come from
step 2's `src/common.py`, so the figures of both steps look the same.

Label inputs are parameterised (the team signed off HER2 rule v1, the exclusions and the group counts on 23 Oct 2026).
Change them in one place, `config.json` next to this folder, or per run with environment variables:
    BT3041_HER2_RULE   = v1 (default; step 2 rule) | v0 (plan rule, `clinical_group_v0`)
    BT3041_HER2_FILTER = none (default) | drop_low_pos (drop HER2+ called on status text only) | drop_low_all
If step 2 changes `tcga_labels.tsv` (e.g. different exclusions), rerunning `run_all.sh` picks it up; no number in
step 3 is typed by hand.
"""
import json, os, sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT3 = Path(__file__).resolve().parents[1]
STEP2 = ROOT3.parent / "step2_preprocessing"
# step 2's package is also called `src`, so load its common.py by path rather than by import name
import importlib.util  # noqa: E402
_spec = importlib.util.spec_from_file_location("step2_common", STEP2 / "src" / "common.py")
_c2 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_c2)
PROC, SEED, GROUP_ORDER, GROUP_COL, PAM50_ORDER, PAM50_COL = (_c2.PROC, _c2.SEED, _c2.GROUP_ORDER, _c2.GROUP_COL,
                                                              _c2.PAM50_ORDER, _c2.PAM50_COL)
CLUSTER_COL, INK, INK2, MUTED, GRID, DIVERGING = _c2.CLUSTER_COL, _c2.INK, _c2.INK2, _c2.MUTED, _c2.GRID, _c2.DIVERGING
style, drop_excluded, RES2 = _c2.style, _c2.drop_excluded, _c2.RES
import matplotlib as mpl  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

RES3 = ROOT3 / "results"
FIG3 = ROOT3 / "figures"
CFG = json.load(open(ROOT3 / "config.json"))
NAME = {"HR+/HER2-": "HR+/HER2−", "HER2+": "HER2+", "TNBC": "TNBC"}
SLIDE = (40 / 3, 7.5)
UP_COL, DOWN_COL = "#1baf7a", "#2a78d6"     # higher in TNBC (TNBC hue) / higher in HR+ (HR+ hue)


def label_config(her2_rule=None, her2_filter=None):
    """Resolve the label configuration: argument > environment variable > config.json."""
    rule = her2_rule or os.environ.get("BT3041_HER2_RULE") or CFG["her2_rule"]
    filt = her2_filter or os.environ.get("BT3041_HER2_FILTER") or CFG["her2_filter"]
    assert rule in ("v1", "v0") and filt in ("none", "drop_low_pos", "drop_low_all"), (rule, filt)
    return rule, filt


def load_labels(her2_rule=None, her2_filter=None, supervised=True):
    """Analysed tumours (step 2 exclusions applied) with a `group` column under the chosen label rule.

    group: HR+/HER2- | HER2+ | TNBC | Incomplete.  her2: Positive | Negative | NaN under the same rule.
    The HER2 confidence filter sets the dropped tumours' group (and HER2) to NaN, so they leave every supervised
    analysis but stay visible as rows.
    """
    rule, filt = label_config(her2_rule, her2_filter)
    lab = drop_excluded(pd.read_csv(PROC / "tcga_labels.tsv", sep="\t", index_col=0))
    if supervised:
        lab = lab[lab.supervised_set]
    lab = lab.copy()
    lab["group"] = lab["clinical_group"] if rule == "v1" else lab["clinical_group_v0"]
    lab["her2"] = lab["HER2"] if rule == "v1" else lab["HER2_v0"]
    low = lab.HER2_confidence.eq("low")
    if filt == "drop_low_pos":
        drop = low & lab.her2.eq("Positive")
    elif filt == "drop_low_all":
        drop = low
    else:
        drop = pd.Series(False, index=lab.index)
    lab.loc[drop, ["group", "her2"]] = np.nan
    lab["label_rule"], lab["her2_filter"] = rule, filt
    # ER-/HER2- vs ER+/HER2- (PR ignored), the contrast matched to METABRIC
    er, h = lab.ER, lab.her2
    lab["erher2"] = np.where(h.eq("Negative") & er.eq("Negative"), "ERneg_HER2neg",
                             np.where(h.eq("Negative") & er.eq("Positive"), "ERpos_HER2neg", None))
    return lab


def write_params3(out_dir, params, name="params.json"):
    """Parameters + run metadata (python, package versions, UTC time, label configuration)."""
    import datetime, importlib.metadata as md
    pk = {}
    for p in ["numpy", "pandas", "scipy", "scikit-learn", "statsmodels", "matplotlib", "seaborn", "pyarrow", "pydeseq2",
              "gseapy", "lifelines", "anndata"]:
        try:
            pk[p] = md.version(p)
        except md.PackageNotFoundError:
            pk[p] = None
    rule, filt = label_config()
    meta = {"python": sys.version.split()[0], "packages": pk, "seed": SEED, "label_rule": rule, "her2_filter": filt,
            "run_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")}
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    json.dump({"params": params, "run": meta}, open(Path(out_dir) / name, "w"), indent=2, default=str)


def dump(obj, path):
    """JSON with numpy types converted."""
    def conv(o):
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            return None if np.isnan(o) else float(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, (pd.Series,)):
            return o.to_dict()
        return str(o)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    json.dump(obj, open(path, "w"), indent=2, default=conv)


def save3(fig, name, exact=False):
    """PNG (300 dpi) + PDF into step 3's figures/. exact=True keeps the canvas (16:9 slides)."""
    FIG3.mkdir(parents=True, exist_ok=True)
    with mpl.rc_context({"savefig.bbox": "standard"} if exact else {}):
        fig.savefig(FIG3 / f"{name}.png")
        fig.savefig(FIG3 / f"{name}.pdf")
    plt.close(fig)
    print("saved", name)


def slide_style():
    style()
    plt.rcParams.update({"font.size": 14, "axes.titlesize": 15, "axes.labelsize": 13, "xtick.labelsize": 12,
                         "ytick.labelsize": 12})


def gene_symbols_star():
    """Versioned Ensembl ID -> HGNC symbol (GENCODE v36 name if unmapped)."""
    ann = pd.read_csv(PROC / "gencode_v36_annotation.tsv", sep="\t", index_col=0)
    return ann.hgnc.fillna(ann.gene), ann


def read_gmt(path=PROC / "h.all.v2024.1.hgnc.gmt"):
    sets = {}
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        sets[p[0]] = [g for g in p[2:] if g]
    return sets


def bh(p):
    from statsmodels.stats.multitest import multipletests
    p = np.asarray(p, float)
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    out[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return out
