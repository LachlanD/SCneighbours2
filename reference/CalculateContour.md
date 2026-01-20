# CalculateContour

Calculates kernel density estimate contour lines for a specific cell
population in a dimensionality reduction space. Uses 2D kernel density
estimation to compute contours at a specified percentile level, useful
for visualizing the spatial distribution of cell populations.

## Usage

``` r
CalculateContour(
  seu,
  meta_data_column,
  meta_data_highlight,
  reduction = "umap",
  percent = 95
)
```

## Arguments

- seu:

  A Seurat object containing single-cell data with dimensionality
  reductions stored in the reductions slot.

- meta_data_column:

  Name of the column in seu@meta.data to pull values from for subsetting
  cells.

- meta_data_highlight:

  The specific value within meta_data_column to calculate contours for.

- reduction:

  Name of the dimensionality reduction to use (default: "umap").

- percent:

  Percentile level for the contour (default: 95). Higher values create
  tighter contours around the densest regions.

## Value

A data frame with x and y coordinates defining the contour line.

## Details

Calculate Density Contour Lines for Cell Populations
