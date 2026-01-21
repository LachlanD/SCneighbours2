# Tutorial

``` r
library(SCneighbours)
library(Seurat)
#> Loading required package: SeuratObject
#> Loading required package: sp
#> 'SeuratObject' was built under R 4.5.0 but the current version is
#> 4.5.2; it is recomended that you reinstall 'SeuratObject' as the ABI
#> for R may have changed
#> 
#> Attaching package: 'SeuratObject'
#> The following objects are masked from 'package:base':
#> 
#>     intersect, t
```

## Example

Here we use 2 different PBMC datasets from SeuratData to demonstrate
usage.

``` r
SeuratData::InstallData("ifnb")
#> Installing package into '/home/runner/work/_temp/Library'
#> (as 'lib' is unspecified)
SeuratData::InstallData("pbmc3k")
#> Installing package into '/home/runner/work/_temp/Library'
#> (as 'lib' is unspecified)

ifnb <- SeuratData::LoadData("ifnb")
#> Validating object structure
#> Updating object slots
#> Ensuring keys are in the proper structure
#> Warning: Assay RNA changing from Assay to Assay
#> Ensuring keys are in the proper structure
#> Ensuring feature names don't have underscores or pipes
#> Updating slots in RNA
#> Validating object structure for Assay 'RNA'
#> Object representation is consistent with the most current Seurat version
#> Warning: Assay RNA changing from Assay to Assay5
pbmc <- SeuratData::LoadData("pbmc3k")
#> Validating object structure
#> Updating object slots
#> Ensuring keys are in the proper structure
#> Warning: Assay RNA changing from Assay to Assay
#> Ensuring keys are in the proper structure
#> Ensuring feature names don't have underscores or pipes
#> Updating slots in RNA
#> Validating object structure for Assay 'RNA'
#> Object representation is consistent with the most current Seurat version
#> Warning: Assay RNA changing from Assay to Assay5

ifnb <- subset(ifnb, stim == 'CTRL')
```

The PBMCs are merged without integration and Run through a standard
pipeline to generate a nearest neighbour graph, clusters, and UMAP
reduction.

``` r
seu <- merge(ifnb, pbmc)
seu <- NormalizeData(seu)
#> Normalizing layer: counts.ifnb
#> Normalizing layer: counts.pbmc3k
seu <- FindVariableFeatures(seu)
#> Finding variable features for layer counts.ifnb
#> Finding variable features for layer counts.pbmc3k
seu <- ScaleData(seu)
#> Centering and scaling data matrix
seu <- RunPCA(seu)
#> PC_ 1 
#> Positive:  FCER1G, TYROBP, FTL, TIMP1, TYMP, C15orf48, LGALS3, ANXA5, SOD2, LGALS1 
#>     NPC2, PSAP, SPI1, PLAUR, CST3, S100A11, ANXA2, CD63, IL8, FTH1 
#>     CTSB, KYNU, CD14, S100A8, SAT1, S100A10, GPX1, HLA-DRA, CFP, S100A6 
#> Negative:  LTB, GIMAP7, CCR7, IL32, CD2, SELL, TSC22D3, TRAT1, GIMAP5, STK17A 
#>     CD247, ITM2A, CD69, CD27, CCL5, GIMAP4, C12orf57, S1PR4, CTSW, ISG20 
#>     PASK, SNHG8, ACAP1, PPA1, CD8B, BIRC3, CREM, NKG7, AQP3, CD8A 
#> PC_ 2 
#> Positive:  FTH1, C15orf48, IL8, CREM, LMNA, GPR183, CCR7, PLA2G7, CCL2, MARCKSL1 
#>     CTSL, CXCL3, PHLDA1, THBS1, IER3, KYNU, CYP1B1, MGST1, ANXA1, TXN 
#>     CXCL2, BIRC3, SOD2, HSPA5, PPIF, INSIG1, ANXA5, CD63, SERPINB2, EIF5 
#> Negative:  FOS, LY6E, IFITM2, ITGB2, DUSP1, LGALS2, HLA-DRB5, PYCARD, AIF1, LST1 
#>     CDA, MT2A, FAM26F, IFITM3, CEBPD, IFI6, CTSS, COTL1, S100A4, LINC00936 
#>     PLAC8, TNFSF10, TNFSF13B, LYZ, CST3, SERPINA1, CFD, MS4A6A, JUN, ISG15 
#> PC_ 3 
#> Positive:  HLA-DQA1, HLA-DQB1, CD79A, HLA-DPA1, CD74, HLA-DRA, HLA-DRB1, HLA-DPB1, MS4A1, CD83 
#>     HLA-DMA, HERPUD1, MIR155HG, SYNGR2, IRF8, CD79B, TCL1A, HLA-DMB, REL, HLA-DQA2 
#>     LINC00926, BLNK, SPIB, HVCN1, HSP90AB1, HLA-DRB5, NME1, FABP5, ID3, EBI3 
#> Negative:  NKG7, CCL5, GZMA, CTSW, GNLY, PRF1, IL32, FGFBP2, ANXA1, GIMAP7 
#>     GZMH, TMSB4X, GZMB, CD247, S100A4, CD2, S100A6, KLRD1, CST7, SPON2 
#>     CLIC3, GIMAP4, GIMAP5, HOPX, S100A8, XCL2, KLRG1, PRSS23, CD8A, S100A9 
#> PC_ 4 
#> Positive:  LTB, LGALS2, S100A9, S100A8, FOS, SELL, RBP7, CDA, AIF1, FCN1 
#>     FOLR3, MS4A6A, LINC00936, LYZ, TSC22D3, TCL1A, CD79A, LST1, LINC00926, CD27 
#>     PASK, MAL, DUSP1, HLA-DRB5, C1orf162, NCF1, FPR1, MS4A1, MNDA, CD79B 
#> Negative:  GZMB, NKG7, CST7, GNLY, FGFBP2, PRF1, GZMA, GZMH, CCL5, CTSW 
#>     KLRD1, CLIC3, SPON2, APOBEC3G, HOPX, CCL4, XCL2, PRSS23, RAMP1, SRGN 
#>     C12orf75, CD247, TTC38, KLRC1, XCL1, IGFBP7, TNFRSF18, AKR1C3, ID2, CTSC 
#> PC_ 5 
#> Positive:  VMO1, MS4A4A, MS4A7, TNFSF10, CXCL10, TMSB4X, FCGR3A, PPM1N, HN1, WARS 
#>     IFITM3, GBP5, GBP1, TNFSF13B, CXCL16, SMPDL3A, OAS1, ATP1B3, FAM26F, FGL2 
#>     CREM, HSPB1, C3AR1, IFIT3, ADA, GBP4, CD86, CDKN1C, SLAMF7, FTH1 
#> Negative:  CD79A, TCL1A, MS4A1, CD79B, CD74, LINC00926, S100A8, NKG7, HLA-DQA2, HLA-DQA1 
#>     HLA-DRB5, FOS, HLA-DQB1, CTSW, HLA-DRA, S100A9, RNASE6, HLA-DPB1, GZMA, CAPG 
#>     PRF1, FGFBP2, GZMB, HLA-DMA, FCRLA, TSPAN13, GNLY, HLA-DMB, CST7, HLA-DRB1

seu <- FindNeighbors(seu, dims = 1:30, reduction = "pca")
#> Computing nearest neighbor graph
#> Computing SNN
seu <- FindClusters(seu, resolution = 0.5)
#> Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
#> 
#> Number of nodes: 9248
#> Number of edges: 391141
#> 
#> Running Louvain algorithm...
#> Maximum modularity in 10 random starts: 0.9226
#> Number of communities: 15
#> Elapsed time: 1 seconds


seu <- RunUMAP(seu, dims = 1:30, reduction = "pca")
#> Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
#> To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
#> This message will be shown once per session
#> 03:07:23 UMAP embedding parameters a = 0.9922 b = 1.112
#> 03:07:23 Read 9248 rows and found 30 numeric columns
#> 03:07:23 Using Annoy for neighbor search, n_neighbors = 30
#> 03:07:23 Building Annoy index with metric = cosine, n_trees = 50
#> 0%   10   20   30   40   50   60   70   80   90   100%
#> [----|----|----|----|----|----|----|----|----|----|
#> **************************************************|
#> 03:07:24 Writing NN index file to temp file /tmp/RtmpySuLRy/file1d057fd9977a
#> 03:07:24 Searching Annoy index using 1 thread, search_k = 3000
#> 03:07:26 Annoy recall = 100%
#> 03:07:27 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
#> 03:07:28 Initializing from normalized Laplacian + noise (using RSpectra)
#> 03:07:28 Commencing optimization for 500 epochs, with 400930 positive edges
#> 03:07:28 Using rng type: pcg
#> 03:07:37 Optimization finished
```

``` r
DimPlot(seu, group.by = c("orig.ident", "seurat_clusters"))
```

![](Tutorial_files/figure-html/unnamed-chunk-4-1.png)

``` r
FeaturePlot(seu, "CD3E")
```

![](Tutorial_files/figure-html/unnamed-chunk-4-2.png)

Here the T cells are split into 2 clusters for the seperate datasets.

Using SC neighbours we can see that there is neighbourhood sharing
between the 2 seperate T cell clusters.

``` r
visualize_neighbourhood(seu, meta_data_column = 'seurat_clusters', meta_data_highlight = 1, 'umap', density = T) +
visualize_neighbourhood(seu, meta_data_column = 'seurat_clusters', meta_data_highlight = 1, 'umap', density = F)
```

![](Tutorial_files/figure-html/unnamed-chunk-5-1.png)

``` r

visualize_neighbourhood(seu, 'seurat_clusters', 2, density = T, percent = 90) +
visualize_neighbourhood(seu, 'seurat_clusters', 2, density = F)
```

![](Tutorial_files/figure-html/unnamed-chunk-5-2.png)

A heatmap showing all neighbourhood sharing of clusters can be generated
with visualise_neighbour_percentage.

``` r
visualise_neighbour_percentage(seu, graph = 'RNA_nn', meta_data_column = 'seurat_clusters')
```

![](Tutorial_files/figure-html/unnamed-chunk-6-1.png)

For each cell we can calculate the percentage of cells from its
neighbours share the cells cluster (or other metadata column).

``` r
seu <- calculate_outside_neighbours_cell(seu, 'seurat_clusters',colname = 'outside_neighbourhood')
FeaturePlot(seu, 'outside_neighbourhood')
```

![](Tutorial_files/figure-html/unnamed-chunk-7-1.png)

Here we see can see the high percentage cells form the boundaries
between clusters, even when those boundaries don’t appear close together
on the UMAP.

For each cell we can calculate the variance in the distance to it it
neighbours based on a reductions embeddings. This can be useful to see
if there are any problems with your UMAP or there are some problamtic
cells in your dataset.

``` r
seu <- calculate_neighbour_distance_for_all_cells(seu, reduction = 'umap', colname = 'neighbour_distance')
FeaturePlot(seu, 'neighbour_distance')
```

![](Tutorial_files/figure-html/unnamed-chunk-8-1.png)
