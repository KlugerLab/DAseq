# DA-seq (Detecting regions of differential abundance between scRNA-seq datasets)

## Introduction
DA-seq is a method to detect cell subpopulations with differential abundance between single cell RNA-seq (scRNA-seq) datasets from different samples, described in the preprint, "Detecting regions of differential abundance between scRNA-Seq datasets" available [here](https://www.biorxiv.org/content/10.1101/711929v2). Given a low dimensional transformation, for example principal component analysis (PCA), of the merged gene expression matrices from different samples (cell states, condition, etc.), DA-seq first computes a score vector for each cell to represent the DA behavior in the neighborhood to select cells in the most DA areas; then groups these cells into distinct DA regions.

This repository contains the DA-seq package.

For more information, please visit our website for DAseq [here](https://klugerlab.github.io/DAseq).