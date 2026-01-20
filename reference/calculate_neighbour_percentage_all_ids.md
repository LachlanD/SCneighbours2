# calculate_neighbour_percentage_all_ids

Calculates the percentage composition of neighbours for all cell
types/groups in the dataset. This creates a matrix showing what
percentage of each cell type's neighbourhood is composed of each other
cell type, useful for understanding global neighbourhood composition
patterns.

## Usage

``` r
calculate_neighbour_percentage_all_ids(
  seu,
  meta_data_column,
  graph = "RNA_nn",
  reduction = NULL
)
```

## Arguments

- seu:

  A Seurat object containing single-cell data with nearest-neighbour
  graphs stored in the graphs slot.

- meta_data_column:

  Name of the column in seu@meta.data to pull values from for grouping
  and analysis.

- graph:

  Name of the nearest-neighbour graph to use from seu@graphs (e.g.,
  "RNA_nn", "RNA_snn", or "SCT_nn").

## Value

A data frame where the first column 'ids' contains all unique cell type
labels, and subsequent columns (named by cell type) contain the
percentage of neighbours belonging to each cell type.

## Details

Calculate Neighbour Percentages for All Cell Types
