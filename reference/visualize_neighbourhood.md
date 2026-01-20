# visualize_neighbourhood

Visualizes the neighbours of cells belonging to a specific group on a
dimensionality reduction plot (e.g., UMAP, t-SNE, PCA). Can display
neighbours either as simple highlighted points or as density contours
showing the spatial distribution patterns. This helps understand how
cell types interact spatially and identify neighbourhood structures.

## Usage

``` r
visualize_neighbourhood(
  seu,
  meta_data_column,
  meta_data_highlight,
  reduction = NULL,
  density = F,
  graph = "RNA_nn",
  percent = 95
)
```

## Arguments

- seu:

  A Seurat object containing single-cell data with nearest-neighbour
  graphs and dimensionality reductions.

- meta_data_column:

  Name of the column in seu@meta.data to pull values from for
  identifying the cells of interest.

- meta_data_highlight:

  The specific value within meta_data_column to analyze neighbours for
  (e.g., a specific cluster or cell type).

- reduction:

  Name of the dimensionality reduction to plot (e.g., "umap", "tsne",
  "pca").

- density:

  Logical indicating whether to plot density contours (TRUE) or simple
  point highlights (FALSE). Default is FALSE.

- graph:

  Name of the nearest-neighbour graph to use from seu@graphs (default:
  "RNA_nn"). Use "RNA_snn" for shared nearest neighbours or "SCT_nn" if
  using SCTransform.

- percent:

  Percentile level for density contours when density=TRUE (default: 95).
  Only used when density=TRUE.

## Value

A ggplot2 object showing the dimensionality reduction with neighbours
highlighted either as red points (density=FALSE) or as density contours
(density=TRUE).

## Details

Visualize Cell Neighbourhoods on Dimensionality Reduction
