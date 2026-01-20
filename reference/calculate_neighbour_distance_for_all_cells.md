# calculate_neighbour_distance_for_all_cells.

The function to calculate the average variance in coordinates of
neighbour cells of all cells in the Seurat object based on the provided
reduction map.

## Usage

``` r
calculate_neighbour_distance_for_all_cells(seu, reduction, colname, graph)
```

## Arguments

- seu:

  The Seurat object.

- reduction:

  The reduction map used to calculate the coordinate variance.

- colname:

  The column name to store the neighbourhood distance value in metadata.

- graph:

  Name of the nearest-neighbour graph to use from seu@graphs.

## Value

A Seurat or SingleCellExperiment object with a new metadata column
storing the variance in coordinates of neighbour cells for each cell
