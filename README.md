# DA-seq (Detecting regions of differential abundance  between scRNA-seq  datasets)

## Introduction
DA-seq is a method to detect cell subpopulations with differential abundance between single cell RNA-seq (scRNA-seq) datasets from different samples, described in the preprint, "Detecting regions of differential abundance between scRNA-Seq datasets" available [here](http://biorxiv.org/cgi/content/short/711929v1). Given a low dimensional transformation, for example principal component analysis (PCA), of the merged gene expression matrices from different samples (cell states, condition, etc.), DA-seq first computes a score vector for each cell to represent the DA behavior in the neighborhood to select cells in the most DA areas; then groups these cells into distinct DA regions.

This repository contains codes for running DA-seq in R.


## Dependencies
Required packages: RANN, tclust, ggplot2, cowplot, RColorBrewer, scales

Suggested package: diffusionMap


## Usage
DA-seq can be used as follows:

Let X be a N-by-p matrix of the PCA embeddings of merged scRNA-seq datasets A (A1 and A2) and B (B1 and B2); X.label be a vector of N specifying the original of each cell ('A1', 'A2', 'B1', or 'B2'); X.2d be the 2D embedding of the cells.

For marker detection of each DA region, let SeuratObj be a Seurat V2.3.0 ([tutorial](https://satijalab.org/seurat/v2.4/pbmc3k_tutorial.html)) object containing the merged datasets (N cells), after proper preprocessing steps (normalization, scaling, etc.).

~~~~
X.da.cells <- getDAcells(
  X = X, 
  cell.labels = X.label, 
  labels.1 = c("A1","A2"), 
  labels.2 = c("B1","B2"), 
  k.vector = seq(50,500,50), 
  plot.embedding = X.2d
)

X.da.regions <- getDAregion(
  X = X, 
  cell.idx = X.da.cells$da.cell.idx, 
  k = 4, alpha = 0.05, 
  cell.labels = X.label, 
  labels.1 = c("A1","A2"), 
  labels.2 = c("B1","B2"), 
  plot.embedding = X.2d
)

X.da.markers <- findMarkersForDAregion(
  cell.idx = X.da.cells$da.cell.idx,
  da.region.label = X.da.regions$cluster.res,
  obj = SeuratObj,
  only.pos = T, min.pct = 0.1, min.diff.pct = 0.09
)
~~~~


A comprehensive example of using DA-seq on scRNA-seq data from [Sade-Feldman, Moshe, et al. (Cell. 2018)](https://www.sciencedirect.com/science/article/pii/S0092867418313941) can be found in `DA-seq_melanoma.R`.

