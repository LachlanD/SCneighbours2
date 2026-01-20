# check_single_cell_object

Validates that the input object is either a Seurat object or a
SingleCellExperiment object. Returns a named list containing
standardized accessors for metadata, graphs, and other components needed
by package functions.

## Usage

``` r
check_single_cell_object(obj, graph, reduction = NULL)
```

## Arguments

- obj:

  An object to validate (should be either Seurat or
  SingleCellExperiment).

## Value

A named list with the following elements:

- `type` - Character string: "Seurat" or "SingleCellExperiment"

- `embeddings` - Tibble of the selected reduction embeddings

- `key` - Character string of embedding dimension name key

- `metadata` - Data frame of cell metadata

- `graphs` - Seurat style nearest neighbour graph

- `n_cells` - Integer, number of cells in the dataset

- `cell_names` - Character vector of cell identifiers

## Details

Check Single-Cell Object Type and Extract Components
