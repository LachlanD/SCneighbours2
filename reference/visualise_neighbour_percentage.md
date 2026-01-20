# visualise_neighbour_percentage

Creates a clustered heatmap showing the percentage of shared neighbours
between different cell types or groups. The heatmap is ordered by
hierarchical clustering to group similar neighbourhood patterns
together. This visualization helps identify which cell types have
similar spatial distributions or microenvironments.

## Usage

``` r
visualise_neighbour_percentage(seu, meta_data_column, graph)
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

A ggplot2 object showing a heatmap of shared neighbour percentages
between cell types, with hierarchical clustering on both axes.

## Details

Visualize Neighbour Percentage Heatmap
