# neighbour_distance.scaled.

The function to calculate the average variance in coordinates of
neighbour cells of a certain cell i based on the provided reduction map.

## Usage

``` r
neighbour_distance.scaled(i, reduction, seu, graph = "RNA_nn")
```

## Arguments

- i:

  The cell index.

- reduction:

  The reduction map used to calculate the coordinate variance.

- seu:

  The Seurat object.

- graph:

  Name of the nearest-neighbour graph to use from seu@graphs (default:
  "RNA_nn").

## Value

A number, the average variance in coordinates of neighbour cells of a
certain cell i based on the provided reduction map.
