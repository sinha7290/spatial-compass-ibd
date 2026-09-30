# spatial-compass-ibd

**Signed, multi-module composite activity scoring for CosMx SMI spatial transcriptomics.**

A small toolkit that computes a **signed composite score per cell** from user-defined gene modules and registers that score onto the matched tissue image, so a program can be read in situ rather than only in embedding space.

> **Status:** method and tooling only. Results from applying this method are part of a manuscript currently under review and are not included here.

---

## Method

For each gene in a module, a per-gene threshold is computed by **StepMiner step-fitting** and used to binarize expression across cells in a field of view. Binarized genes accumulate into a per-cell module score. Modules are then combined as a **signed weighted sum**, producing one directional value per cell.

A positive module weight marks a program the score should read as present; a negative weight marks an opposing one. This makes the output directional: rather than "how much of this program is here", the score expresses which of several weighted programs dominates a given cell.

Scores are overlaid on the CellComposite image after coordinate registration (flip / rotation), with the tissue rendered in grayscale so that colour carries only the score.

Genes absent from the panel, or with zero variance across cells, are **skipped rather than imputed** — the score reflects only what the assay actually measured. Binarizing per gene and per field of view, rather than applying a global cutoff, keeps the score robust to the intensity variation between fields and samples that otherwise dominates spatial panels.

Full detail, parameters and citations: [`docs/METHOD.md`](docs/METHOD.md).

## Data

Built against the public dataset **[GSE234713](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234713)** (CosMx SMI, human colon). **No data is committed to this repository** — the raw download is tens of GB.

```bash
bash scripts/download_GSE234713.sh ./RAW
```

## Install

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

`StepMiner` must be importable on your `PYTHONPATH` — see [`docs/METHOD.md`](docs/METHOD.md) for the citation and source.

## Use

```python
from src.spatial_compass import load_expr_meta, composite_score

df  = load_expr_meta("./RAW", "GSM7473684_HC_c")
fov = df[df["fov"] == 1]

gene_cols = [c for c in fov.columns if c not in
             ("fov", "cell_ID", "CenterX_local_px", "CenterY_local_px")]

# two modules, opposing weights -> one directional score per cell
score = composite_score(
    fov[gene_cols],
    modules=[["GENE_A", "GENE_B"], ["GENE_C"]],
    weights=[1.0, -1.0],
)
```

An interactive notebook — sample / field-of-view / module / weight selection with live overlay — is in [`notebooks/`](notebooks/).

## Layout

```
src/spatial_compass.py   scoring, registration and IO
notebooks/               interactive composite-score explorer
scripts/                 GSE234713 download helper
docs/METHOD.md           method detail, parameters and citations
```

## Citation

Please cite the dataset (GSE234713) and the StepMiner method as listed in [`docs/METHOD.md`](docs/METHOD.md).

## License

MIT — see [LICENSE](LICENSE).
