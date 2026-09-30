"""
spatial_compass — composite module scoring for CosMx SMI spatial transcriptomics.

Computes a signed, multi-module composite activity score per cell and registers it
onto the matched CellComposite tissue image.

Method summary
--------------
1. Per-cell expression and metadata (cell centroids) are joined per field of view.
2. For each gene in a module, a per-gene threshold is computed by StepMiner
   step-fitting and used to binarize expression across cells.
3. Binarized genes are accumulated into a per-cell score, with each module
   contributing at a user-specified signed weight (e.g. programme A +1, programme B -1).
4. Scores are overlaid on the tissue image after coordinate registration
   (flip / rotation), so the score is read in situ rather than in embedding space.

StepMiner thresholding is an established published method and is used here as a
dependency; see docs/METHOD.md for the citation. The contribution in this
repository is the signed multi-module composite score and its spatial registration.
"""
from __future__ import annotations

import glob
import os

import numpy as np
import pandas as pd

try:
    import StepMiner
except ImportError:  # pragma: no cover
    StepMiner = None


# ----------------------------------------------------------------------------- IO
def list_samples(raw_dir: str) -> list[str]:
    """Sample names discovered from *_exprMat_file.csv.gz in raw_dir."""
    expr_files = sorted(glob.glob(os.path.join(raw_dir, "*_exprMat_file.csv.gz")))
    return [os.path.basename(x).replace("_exprMat_file.csv.gz", "") for x in expr_files]


def image_folder_name(sample_name: str) -> str:
    """GSM7473684_HC_c -> GSM7473684_HCc (GEO image folders drop the underscore)."""
    gsm, disease, rep = sample_name.split("_")
    return f"{gsm}_{disease}{rep}"


def load_expr_meta(raw_dir: str, sample_name: str) -> pd.DataFrame:
    """Join per-cell expression with metadata on (fov, cell_ID)."""
    expr = pd.read_csv(
        os.path.join(raw_dir, f"{sample_name}_exprMat_file.csv.gz"), compression="gzip"
    )
    meta = pd.read_csv(
        os.path.join(raw_dir, f"{sample_name}_metadata_file.csv.gz"), compression="gzip"
    )
    required = ["fov", "cell_ID", "CenterX_local_px", "CenterY_local_px"]
    missing = [c for c in required if c not in meta.columns]
    if missing:
        raise ValueError(f"Missing required metadata columns: {missing}")
    return expr.merge(meta, on=["fov", "cell_ID"], how="inner")


# ------------------------------------------------------------------- registration
def transform_coords(x, y, img_shape, mode: str = "flip_y"):
    """Register cell centroids onto the tissue image coordinate frame."""
    h, w = img_shape[:2]
    x = np.asarray(x, dtype=float).copy()
    y = np.asarray(y, dtype=float).copy()
    if mode == "none":
        return x, y
    if mode == "flip_x":
        return (w - x), y
    if mode == "flip_y":
        return x, (h - y)
    if mode == "flip_xy":
        return (w - x), (h - y)
    if mode == "rot90":
        return (h - y), x
    if mode == "rot180":
        return (w - x), (h - y)
    if mode == "rot270":
        return y, (w - x)
    raise ValueError(f"Unknown transform mode: {mode}")


# ------------------------------------------------------------------------ scoring
def stepminer_threshold(values) -> float:
    """Per-gene binarization threshold via StepMiner step-fitting."""
    arr = np.array([float(v) for v in values if pd.notnull(v)], dtype=float)
    arr.sort()
    if arr.size == 0:
        return float("nan")
    if arr.size < 3:
        return float(np.median(arr))
    if StepMiner is None:
        raise ImportError(
            "StepMiner is required for thresholding. See docs/METHOD.md."
        )
    return float(StepMiner.fitstep(arr)["threshold"])


def score_module(expr_cells_by_genes: pd.DataFrame, genes) -> np.ndarray:
    """Binarize each gene at its StepMiner threshold and accumulate per cell.

    Genes absent from the panel, or with zero variance across cells, are skipped
    rather than imputed — the score reflects only genes the assay actually measured.
    """
    acc = np.zeros(expr_cells_by_genes.shape[0], dtype=float)
    cols_upper = {c.upper(): c for c in expr_cells_by_genes.columns}
    for g_raw in genes:
        g = str(g_raw).upper()
        if g not in cols_upper:
            continue
        v = expr_cells_by_genes[cols_upper[g]].astype(float).values
        sd = float(np.std(v))
        if not np.isfinite(sd) or sd <= 0:
            continue
        thr = stepminer_threshold(v)
        if not np.isfinite(thr):
            continue
        acc += (v > thr).astype(float)
    return acc


def composite_score(
    expr_cells_by_genes: pd.DataFrame, modules: list[list[str]], weights: list[float]
) -> np.ndarray:
    """Signed weighted sum of module scores — one composite value per cell.

    Parameters
    ----------
    modules : list of gene lists, e.g. [["GENE_A", "GENE_B"], ["GENE_C"]]
    weights : signed weight per module, e.g. [1.0, -1.0]

    A positive weight marks a program the score should read as present; a negative
    weight marks one it should read as opposing. The result is directional: high
    scores indicate the positively weighted program dominating in that cell.
    """
    if len(modules) != len(weights):
        raise ValueError("modules and weights must be the same length")
    total = np.zeros(expr_cells_by_genes.shape[0], dtype=float)
    for genes, w in zip(modules, weights):
        total += float(w) * score_module(expr_cells_by_genes, genes)
    return total


def to_grayscale(img: np.ndarray) -> np.ndarray:
    """Luminance conversion, so score colour is the only chromatic signal."""
    if img.ndim == 2:
        return img
    return (0.2989 * img[..., 0] + 0.5870 * img[..., 1] + 0.1140 * img[..., 2])
