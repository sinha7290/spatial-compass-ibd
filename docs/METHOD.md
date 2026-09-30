# Method

## 1. Input

CosMx SMI per-cell output, joined on `(fov, cell_ID)`:

- `*_exprMat_file.csv.gz` — cells × genes counts
- `*_metadata_file.csv.gz` — must contain `CenterX_local_px`, `CenterY_local_px`
- `CellComposite/CellComposite_F<fov>.jpg` — matched tissue image per field of view

## 2. Per-gene binarization

For gene *g* with expression vector **v** across cells in the field of view:

1. Genes not present in the panel are skipped.
2. Genes with `std(v) <= 0` are skipped (no information to threshold).
3. Otherwise a threshold *t_g* is obtained by **StepMiner step-fitting** on the sorted values, and the gene is binarized as `v > t_g`.

Fewer than three non-null values falls back to the median.

Binarizing per gene, per field of view, rather than applying a global cutoff, keeps the score robust to the intensity variation between fields and samples that otherwise dominates spatial panels.

## 3. Module score

A module *M* is a gene list. Its score for a cell is the count of binarized-positive genes in *M* for that cell.

## 4. Signed composite

Given modules `M1..Mk` and signed weights `w1..wk`:

```
composite(cell) = Σ  w_i · score_{M_i}(cell)
```

A positive weight marks a program the score should read as present; a negative weight marks an opposing one. The output is directional — high values indicate the positively weighted program dominating that cell.

Example shape: `modules=[["GENE_A","GENE_B"], ["GENE_C"]]`, `weights=[1.0, -1.0]` contrasts two opposing programs in a single value.

## 5. Spatial registration

Cell centroids are given in local pixel coordinates but do not necessarily share the image's origin convention. `transform_coords` applies one of `none`, `flip_x`, `flip_y`, `flip_xy`, `rot90`, `rot180`, `rot270`. For this dataset `flip_y` aligns centroids to the CellComposite frame.

The tissue image is converted to luminance so that colour carries only the score. Scores are clipped at an upper quantile (default 0.995) before mapping to colour, so a small number of extreme cells do not flatten the visible dynamic range.

## 6. Parameters

| Parameter | Default | Purpose |
|---|---|---|
| `transform_mode` | `flip_y` | centroid-to-image registration |
| `clip_quantile` | `0.995` | upper clip before colour mapping |
| `crop_quantile` | `0.02` | trims empty margin around occupied tissue |
| `background_mode` | `gray` | luminance tissue background |
| `point_size` | `8` | centroid marker size |

## Citations

- **Dataset** — GSE234713, CosMx SMI of human colon. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234713
- **StepMiner** — Sahoo D, Dill DL, Tibshirani R, Plevritis SK. *Extracting binary signals from microarray time-course data.* Nucleic Acids Research 35(11):3705–3712, 2007.

StepMiner is used here as an established published dependency. The contribution in this repository is the signed multi-module composite score and its registration onto tissue images.
