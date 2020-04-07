#' DA-seq Step 3: get DA regions
#'
#' Cluster the DA cells retained from Step 1 and Step 2 of DA-seq with outlier removal to obtain spatially
#' coherent DA regions.
#'
#' @param X size N-by-p matrix, input merged dataset of interest after dimension reduction
#' @param cell.idx result "da.cell.idx" from the output of function getDAcells
#' @param k number of DA regions to find for cells from function getDAcells
#' @param alpha estimated ratio of outliers of cells from function getDAcells
#' @param restr.fact parameter inherited from function "tclust"
#' @param cell.labels size N vector, labels for each input cell
#' @param labels.1 vector, label name(s) that represent condition 1
#' @param labels.2 vector, label name(s) that represent condition 2
#' @param do.plot a logical value to indicate whether to return ggplot objects showing the results, default True
#' @param plot.embedding size N-by-2 matrix, 2D embedding for the cells
#' @param size cell size to use in the plot, default 0.5
#' @param ... other parameters passed to function tclust
#'
#' @import scales
#' @import tclust
#'
#' @return a list of results
#' \describe{
#'   \item{cluster.res}{DA region number for each cell from cell.idx, '0' represents outlier cells}
#'   \item{DA.stat}{a table showing DA score and p-value for each DA region}
#'   \item{da.region.plot}{ggplot object showing DA regions on plot.embedding}
#' }
#'
#' @export

getDAregion <- function(
  X, cell.idx, k, alpha, restr.fact = 50,
  cell.labels, labels.1, labels.2,
  do.plot = T, plot.embedding = NULL, size = 0.5,
  seed = 0, ...
){
  set.seed(seed)
  X.tclust <- runtclust(X, cell.idx, k, alpha, restr.fact, ...)
  X.n.da <- length(unique(X.tclust)) - 1
  X.da.stat <- matrix(0, nrow = X.n.da, ncol = 3)
  colnames(X.da.stat) <- c("DA.score","pval.wilcoxon","pval.ttest")
  for(ii in 1:X.n.da){
    X.da.stat[ii,] <- getDAscore(
      cell.labels = cell.labels, cell.idx = cell.idx[X.tclust == ii],
      labels.1 = labels.1, labels.2 = labels.2
    )
  }

  if(do.plot & is.null(plot.embedding)){
    warning("plot.embedding must be provided by user if do.plot = T")
    X.region.plot <- NULL
  } else if(do.plot & !is.null(plot.embedding)){
    X.da.label <- rep(0,nrow(X))
    X.da.label[cell.idx] <- X.tclust
    X.da.order <- order(X.da.label, decreasing = F)
    X.region.plot <- plotCellLabel(
      X = plot.embedding[X.da.order,], label = as.factor(X.da.label[X.da.order]),
      size = size, do.label = F, return.plot = T,
    ) + scale_color_manual(values = c("gray",hue_pal()(X.n.da)), breaks = c(1:X.n.da))
    # X.region.plot <- plotCellLabel(
    #   X = plot.embedding[cell.idx,], label = as.factor(X.tclust),
    #   size = size, do.label = F, return.plot = T
    # ) + scale_color_manual(values = c(rgb(255,255,255,max = 255,alpha = 0),hue_pal()(X.n.da)),
    #                        breaks = c(1:X.n.da))
  } else {
    X.region.plot <- NULL
  }

  return(list(
    cluster.res = X.tclust,
    DA.stat = X.da.stat,
    da.region.plot = X.region.plot
  ))
}



# Run tclust
runtclust <- function(X, cell.idx, k, alpha, restr.fact = 50, ...){
  X.tclust.res <- tclust(
    x = X[cell.idx,], k = k, alpha = alpha, restr.fact = restr.fact, ...
  )

  return(X.tclust.res$cluster)
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


