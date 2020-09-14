# DA-seq (Detecting regions of differential abundance between scRNA-seq datasets)

## Introduction
DA-seq is a method to detect cell subpopulations with differential abundance between single cell RNA-seq (scRNA-seq) datasets from different samples, described in the preprint, "Detecting regions of differential abundance between scRNA-Seq datasets" available [here](https://www.biorxiv.org/content/10.1101/711929v2). Given a low dimensional transformation, for example principal component analysis (PCA), of the merged gene expression matrices from different samples (cell states, condition, etc.), DA-seq first computes a score vector for each cell to represent the DA behavior in the neighborhood to select cells in the most DA areas; then groups these cells into distinct DA regions.

[This](https://github.com/KlugerLab/DAseq) repository contains the DA-seq package.

## R Dependencies
Required packages: RANN, glmnet, caret, Seurat, e1071, reticulate, ggplot2, cowplot, scales, ggrepel

## Python Dependencies
Python 3 or above ([Miniconda](https://docs.conda.io/en/latest/miniconda.html) is recommended.)

Required modules: numpy, pandas, sklearn, tensorflow, keras


## Installation
The DAseq package can be installed from this GitHub repository.

```r
devtools::install_github("KlugerLab/DAseq")
```


## References
References of DAseq functions can be found [here](https://klugerlab.github.io/DAseq/reference/index.html).


## Usage
Please check DA-seq [tutorial](https://klugerlab.github.io/DAseq/articles/tutorial.html).

Data used in the tutorial is from [Sade-Feldman, Moshe, et al. (Cell. 2018)](https://www.sciencedirect.com/science/article/pii/S0092867418313941).
