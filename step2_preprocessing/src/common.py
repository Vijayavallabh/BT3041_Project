"""Shared paths, constants and plotting style for the BT3041 TCGA-BRCA project."""
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROC = ROOT / "data" / "processed"
FIG = ROOT / "figures"
RES = ROOT / "results"
SEED = 42

# One fixed colour per clinical group, used in every figure (validated all-pairs for CVD).
GROUP_ORDER = ["HR+/HER2-", "HER2+", "TNBC"]
GROUP_COL = {"HR+/HER2-": "#2a78d6", "HER2+": "#eb6834", "TNBC": "#1baf7a",
             "Incomplete": "#c3c2b7", "Normal": "#52514e"}
PAM50_ORDER = ["LumA", "LumB", "Her2", "Basal", "Normal"]
# Luminal subtypes share the blue family, Her2-enriched/Basal reuse the HER2+/TNBC hues they mostly overlap with.
PAM50_COL = {"LumA": "#184f95", "LumB": "#6da7ec", "Her2": "#eb6834",
             "Basal": "#1baf7a", "Normal": "#52514e"}
# Clusters and methods use hexes that differ from the clinical-group / PAM50 colours and never appear in the same
# figure as them, so a cluster is not mistaken for a receptor group.
CLUSTER_COL = ["#6250d6", "#eda100", "#e87ba4", "#008300", "#e34948", "#7a5c3e", "#00a3b4", "#898781"]
# one fixed colour per clustering run, shared by fig9 and slide 5 (k = 5 matches its colour in fig7)
METHOD_COL = {"kmeans_k2": CLUSTER_COL[0], "ward_k2": CLUSTER_COL[1], "avgcorr_k2": CLUSTER_COL[2], "kmeans_k5": CLUSTER_COL[3],
              "consensus_k2": CLUSTER_COL[4]}
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
DIVERGING = "RdBu_r"

MARKERS = ["ESR1", "PGR", "ERBB2", "GATA3", "FOXA1", "FOXC1", "KRT5", "MKI67"]


def write_params(out_dir, params, report=None, name="params.json"):
    """Write the parameters of a step together with run metadata (python, package versions, UTC time)."""
    import json, sys, datetime, importlib.metadata as md
    pk = {}
    for p in ["numpy", "pandas", "scipy", "scikit-learn", "statsmodels", "umap-learn", "matplotlib", "seaborn", "pyarrow"]:
        try:
            pk[p] = md.version(p)
        except md.PackageNotFoundError:
            pk[p] = None
    meta = {"python": sys.version.split()[0], "packages": pk, "seed": SEED,
            "run_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")}
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    json.dump({"params": params, "run": meta}, open(Path(out_dir) / name, "w"), indent=2)


def drop_excluded(lab):
    """Tumours used in every analysis: drop identity-QC failures and TCGA-CDR redacted cases."""
    return lab[~lab.exclude]


def style():
    mpl.rcParams.update({
        "figure.dpi": 110, "savefig.dpi": 300, "savefig.bbox": "tight",
        "font.family": "DejaVu Sans", "font.size": 10, "axes.titlesize": 11,
        "axes.titleweight": "bold", "axes.labelsize": 10,
        "axes.edgecolor": "#c3c2b7", "axes.labelcolor": INK2, "axes.titlecolor": INK,
        "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False,
        "axes.spines.right": False, "axes.grid": False, "grid.color": GRID,
        "grid.linewidth": 0.6, "legend.frameon": False, "figure.facecolor": "white",
        "axes.facecolor": "white", "lines.linewidth": 2,
    })


def save(fig, name, exact=False):
    """Save PNG (300 dpi) and PDF. exact=True keeps the figure's own canvas size (no tight crop), for 16:9 slide figures."""
    FIG.mkdir(parents=True, exist_ok=True)
    with mpl.rc_context({"savefig.bbox": "standard"} if exact else {}):
        fig.savefig(FIG / f"{name}.png")
        fig.savefig(FIG / f"{name}.pdf")
    plt.close(fig)
    print("saved", name)
