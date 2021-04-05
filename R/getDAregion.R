#' DA-seq Step 3: get DA regions
#'
#' Cluster the DA cells retained from Step 1 and Step 2 of DA-seq to obtain spatially
#' coherent DA regions.
#'
#' @param X size N-by-p matrix, input merged dataset of interest after dimension reduction
#' @param da.cells output from getDAcells() or updateDAcells()
#' @param cell.labels size N vector, labels for each input cell
#' @param labels.1 vector, label name(s) that represent condition 1
#' @param labels.2 vector, label name(s) that represent condition 2
#' @param prune.SNN parameter for Seurat function FindNeighbors(), default 1/15
#' @param resolution parameter for Seurat function FindClusters(), default 0.05
#' @param group.singletons parameter for Seurat function FindClusters(), default True
#' @param min.cell integer, number of cells below which a DA region will be removed as outliers, default NULL, use minimum k value in k-vector
#' @param do.plot a logical value to indicate whether to return ggplot objects showing the results, default True
#' @param plot.embedding size N-by-2 matrix, 2D embedding for the cells
#' @param size cell size to use in the plot, default 0.5
#' @param do.label a logical value to indicate whether to label each DA region with text, default False

#' @param ... other parameters to pass to Seurat FindClusters()
#'
#' @import Seurat
#' @import scales
#'
#' @return a list of results
#' \describe{
#'   \item{da.region.label}{DA region label for each cell from the whole dataset,
#'   '0' represents non-DA cells.}
#'   \item{DA.stat}{a table showing DA score and p-value for each DA region}
#'   \item{da.region.plot}{ggplot object showing DA regions on plot.embedding}
#' }
#'
#' @export
#'
getDAregion <- function(
  X, da.cells,
  cell.labels, labels.1, labels.2,
  prune.SNN = 1/15, resolution = 0.05, group.singletons = F, min.cell = NULL,
  do.plot = T, plot.embedding = NULL, size = 0.5, do.label = F,
  ...
){
  if(!inherits(x = X, what = "matrix")){
    cat("Turning X to a matrix.\n")
    X <- as.matrix(X)
  }
  n.cells <- nrow(X)
  n.dims <- ncol(X)
  if(is.null(rownames(X))){
    rownames(X) <- paste("C",c(1:n.cells), sep = "")
  }
  # check label input
  if(!inherits(cell.labels, "character") |
     !inherits(labels.1, "character") | !inherits(labels.2, "character")){
    stop("Input parameters cell.labels, labels.1 and labels.2 must be character")
  }
  if(length(setdiff(cell.labels, c(labels.1, labels.2))) > 0){
    stop("Input parameter cell.labels contain labels not from labels.1 or labels.2")
  }
  if(is.null(min.cell)){
    min.cell <- as.integer(colnames(da.cells$da.ratio)[1])
    cat("Using min.cell = ", min.cell, "\n", sep = "")
  }

  seurat.version <- substr(packageVersion("Seurat"),1,1)
  if(seurat.version == "3" | seurat.version == "4"){
    X.S <- CreateSeuratObject(counts = t(X))
    X.S@reductions$pca <- new(
      "DimReduc", cell.embeddings = X,
      assay.used = DefaultAssay(X.S), key = "PC_"
    )
    X.S <- FindNeighbors(X.S, reduction = "pca", dims = 1:n.dims, prune.SNN = prune.SNN, verbose = F)
    if(length(da.cells$da.up) > 1){
      up.S <- CreateSeuratObject(
        counts = t(X[da.cells$da.up,])
      )
      up.S@reductions$pca <- new(
        "DimReduc", cell.embeddings = X[da.cells$da.up,],
        assay.used = DefaultAssay(up.S), key = "PC_"
      )
      up.S <- FindNeighbors(up.S, reduction = "pca", dims = 1:n.dims, verbose = F)
      up.snn <- X.S@graphs$RNA_snn[da.cells$da.up,da.cells$da.up]
      up.S@graphs$RNA_snn@i <- up.snn@i
      up.S@graphs$RNA_snn@p <- up.snn@p
      up.S@graphs$RNA_snn@x <- up.snn@x
      up.S <- FindClusters(up.S, resolution = resolution, group.singletons = group.singletons, verbose = F, ...)
      up.clusters <- as.numeric(up.S@active.ident)
      up.clusters[up.S@active.ident == "singleton"] <- 0
    } else {
      up.clusters <- NULL
    }
    n.up.clusters <- length(unique(up.clusters)) - as.numeric(0 %in% up.clusters)

    if(length(da.cells$da.down) > 1){
      down.S <- CreateSeuratObject(
        counts = t(X[da.cells$da.down,])
      )
      down.S@reductions$pca <- new(
        "DimReduc", cell.embeddings = X[da.cells$da.down,],
        assay.used = DefaultAssay(down.S), key = "PC_"
      )
      down.S <- FindNeighbors(down.S, reduction = "pca", dims = 1:n.dims, verbose = F)
      down.snn <- X.S@graphs$RNA_snn[da.cells$da.down,da.cells$da.down]
      down.S@graphs$RNA_snn@i <- down.snn@i
      down.S@graphs$RNA_snn@p <- down.snn@p
      down.S@graphs$RNA_snn@x <- down.snn@x
      down.S <- FindClusters(
        down.S, resolution = resolution, group.singletons = group.singletons, verbose = F, ...
      )
      down.clusters <- as.numeric(down.S@active.ident) + n.up.clusters
      down.clusters[down.S@active.ident == "singleton"] <- 0
    } else {
      down.clusters <- NULL
    }
  }
  da.region.label <- rep(0, n.cells)
  da.region.label[da.cells$da.up] <- up.clusters
  da.region.label[da.cells$da.down] <- down.clusters

  # remove small clusters with cells < min.cell
  da.region.label.tab <- table(da.region.label)
  if(min(da.region.label.tab) < min.cell){
    da.region.to.remove <- as.numeric(names(da.region.label.tab)[which(da.region.label.tab < min.cell)])
    cat("Removing ", length(da.region.to.remove), " DA regions with cells < ", min.cell, ".\n", sep = "")
    da.region.label.old <- da.region.label
    for(ii in da.region.to.remove){
      da.region.label[da.region.label.old == ii] <- 0
      da.region.label[da.region.label.old > ii] <- da.region.label[da.region.label.old > ii] - 1
    }
  }

  X.n.da <- length(unique(da.region.label)) - 1
  X.da.stat <- matrix(0, nrow = X.n.da, ncol = 3)
  colnames(X.da.stat) <- c("DA.score","pval.wilcoxon","pval.ttest")
  if(X.n.da > 0){
    for(ii in 1:X.n.da){
      X.da.stat[ii,] <- getDAscore(
        cell.labels = cell.labels, cell.idx = which(da.region.label == ii),
        labels.1 = labels.1, labels.2 = labels.2
      )
    }
  } else {
    warning("No DA regions found.")
  }

  if(do.plot & is.null(plot.embedding)){
    warning("plot.embedding must be provided by user if do.plot = T")
    X.region.plot <- NULL
  } else if(do.plot & !is.null(plot.embedding)){
    X.da.label <- da.region.label
    X.da.order <- order(X.da.label, decreasing = F)
    X.region.plot <- plotCellLabel(
      X = plot.embedding[X.da.order,], label = as.factor(X.da.label[X.da.order]),
      size = size, do.label = do.label, label.plot = as.character(c(1:X.n.da))
    ) + scale_color_manual(values = c("gray",hue_pal()(X.n.da)), breaks = c(1:X.n.da))
  } else {
    X.region.plot <- NULL
  }

  return(list(
    da.region.label = da.region.label,
    DA.stat = X.da.stat,
    da.region.plot = X.region.plot
  ))
}



# Calculate DA statistics
getDAscore <- function(cell.labels, cell.idx, labels.1, labels.2){
  labels.1 <- labels.1[labels.1 %in% cell.labels]
  labels.2 <- labels.2[labels.2 %in% cell.labels]

  idx.label <- cell.labels[cell.idx]
  ratio.1 <- sum(idx.label %in% labels.1) / sum(cell.labels %in% labels.1)
  ratio.2 <- sum(idx.label %in% labels.2) / sum(cell.labels %in% labels.2)
  ratio.diff <- (ratio.2 - ratio.1) / (ratio.2 + ratio.1)

  cell.label.name <- sort(unique(cell.labels))
  cell.label.tab <- table(factor(cell.labels, levels = cell.label.name))

  idx.label.ratio <- table(factor(idx.label, levels = cell.label.name)) / cell.label.tab
  # print(idx.label.ratio)
  # score <- (mean(idx.label.ratio[labels.2]) - mean(idx.label.ratio[labels.1]))
  # score.n <- (mean(idx.label.ratio[labels.2]) - mean(idx.label.ratio[labels.1])) / sum(idx.label.ratio)

  if(length(labels.1) > 1 & length(labels.2) > 1){
    pval.wilcox <- wilcox.test(x = idx.label.ratio[labels.2], idx.label.ratio[labels.1])$p.value
    pval.ttest <- t.test(x = idx.label.ratio[labels.2], idx.label.ratio[labels.1])$p.value
  } else {
    pval.wilcox <- NA
    pval.ttest <- NA
  }

  return(c("DA.score" = ratio.diff, "pval.wilcoxon" = pval.wilcox, "pval.ttest" = pval.ttest))
}


