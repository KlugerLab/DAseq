# DA-seq (Detecting regions of differential abundance  between scRNA-seq  datasets)

## Introduction
DA-seq is a method to detect cell subpopulations with differential abundance between single cell RNA-seq (scRNA-seq) datasets from different samples. Given a low dimensional transformation, for example principal component analysis (PCA), of the merged gene expression matrices, DA-seq first computes a score vector for each cell to represent the DA behavior in the neighborhood to select cells in the most DA areas; then groups these cells into distinct DA regions.

This repository contains codes for running DA-seq in R.


## Dependencies
Required packages: RANN, tclust, ggplot2, cowplot, RColorBrewer, scales

Suggested package: diffusionMap


## Usage
DA-seq can be used as follows:

Let X be a N-by-p matrix of the PCA embeddings of merged scRNA-seq datasets A and B; X.label be a vector of N specifying the original of each cell ('A' or 'B'); X.2d be the 2D embedding of the cells.

~~~~
X.da.cells <- getDAcells(
  X = X, 
  cell.labels = X.label, 
  labels.1 = "A", 
  labels.2 = "B", 
  k.vector = seq(50,500,10), 
  plot.embedding = X.2d
)

X.da.regions <- getDAregion(
  X = X, 
  cell.idx = X.da.cells$da.cell.idx, 
  k = 4, alpha = 0.05, 
  cell.labels = X.label, 
  labels.1 = "A", 
  labels.2 = "B", 
  plot.embedding = X.2d
)
~~~~
