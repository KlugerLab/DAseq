library(RANN)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(scales)
library(tclust)



## Step1: select cells with DA neighborhood

#' @param X size N-by-p matrix, input merged dataset of interest after dimension reduction
#' @param cell.labels size N vector, labels for each input cell
#' @param labels.1 vector, label name(s) that represent condition 1
#' @param labels.2 vector, label name(s) that represent condition 2
#' @param k.vector vector, k values to create the score vector
#' @param ratio maximum ratio of cells to keep, default 0.2
#' @param i.max maximum iteration to run in the iterative clustering of the score vector, default 10
#' @param do.diffuse a logical value to indicate whether to calculate diffusion coordinates for X, default False
#' @param neigen number of diffusion coordinates, default 20
#' @param do.plot a logical value to indicate whether to return ggplot objects showing the results, default True
#' @param plot.embedding size N-by-2 matrix, 2D embedding for the cells
#' @param size cell size to use in the plot, default 0.5
#' 
#' @return a list of results
#'         da.ratio: score vector for each cell
#'         da.cell.idx: cell index with the most DA neighborhood
#'         da.plot: ggplot object showing the steps of iterative clustering result on plot.embedding
#'         da.cells.plot: ggplot object highlighting cells of da.cell.idx on plot.embedding

getDAcells <- function(
  X, cell.labels, labels.1, labels.2, k.vector,
  ratio = 0.2, i.max = 10, 
  do.diffuse = F, neigen = 20, do.plot = T, plot.embedding = NULL, size = 0.5
){
  # get diffusion coordinates
  if(do.diffuse){
    library(diffusionMap)
    cat("Calculating diffusion coordinates.\n")
    X.input <- X
    X <- diffuse(D = dist(X.input), neigen = neigen)$X
  }
  
  # get DA score vector for each cell
  cat("Calculating DA score vector.\n")
  X.knn.ratio <- daPerCell(
    X = X, 
    cell.labels = cell.labels, 
    labels.1 = labels.1, 
    labels.2 = labels.2, 
    k.vector = k.vector
  )
  
  # select DA cells
  cat("Selecting top DA cells.\n")
  X.da.res <- removeNeutralCells(
    X = X.knn.ratio, ratio = ratio, i.max = i.max, keep.info = T
  )
  X.da.idx <- X.da.res[[1]]
  
  # get top down clustering process plot
  if(do.plot & is.null(plot.embedding)){
    warning("plot.embedding must be provided by user if do.plot = T")
    X.da.plot <- NULL
    X.da.cells.plot <- NULL
  } else if(do.plot & !is.null(plot.embedding)){
    X.da.plot <- plotCellLabel(
      X = plot.embedding, label = as.factor(X.da.res[[2]]), 
      cell.col = c("black", brewer.pal(length(unique(X.da.res[[2]])),name = "Blues")[-1]),
      size = size, do.label = F, return.plot = T
    )
    X.da.cells.plot <- plotDAsite(
      X = plot.embedding, 
      site.list = list(X.da.idx), 
      size = size
    )
  } else {
    X.da.plot <- NULL
    X.da.cells.plot <- NULL
  }
  
  # return result
  return(list(
    da.ratio = X.knn.ratio,
    da.cell.idx = X.da.idx,
    da.plot = X.da.plot,
    da.cells.plot = X.da.cells.plot
  ))
}



## Step2: get DA regions from DA cells in Step1

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
#' 
#' @return a list of results
#'         cluster.res: DA region number for each cell from cell.idx, '0' represents outlier cells
#'         DA.stat: a table showing DA score and p-value for each DA region
#'         da.region.plot: ggplot object showing DA regions (cells in gray are outliers) on plot.embedding

getDAregion <- function(
  X, cell.idx, k, alpha, restr.fact = 50,
  cell.labels, labels.1, labels.2, 
  do.plot = T, plot.embedding = NULL, size = 0.5
){
  X.tclust <- runtclust(X, cell.idx, k, alpha, restr.fact)
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
    X.region.plot <- plotCellLabel(
      X = plot.embedding[cell.idx,], label = as.factor(X.tclust), 
      size = size, do.label = F, return.plot = T
    ) + scale_color_manual(values = c(rgb(255,255,255,max = 255,alpha = 0),hue_pal()(X.n.da)), breaks = c(1:X.n.da))
  } else {
    X.region.plot <- NULL
  }
  
  return(list(
    cluster.res = X.tclust,
    DA.stat = X.da.stat,
    da.region.plot = X.region.plot
  ))
}



## Step 3: detect genes that characterize DA regions from Step 2

#' @param cell.idx result "da.cell.idx" from the output of function getDAcells
#' @param da.region.label result "cluster.res" from the output of function getDAregion
#' @param obj Seurat object that contain ALL cells in the analysis
#' @param ... parameters for Seurat function FindMarkers()
#' 
#' @return a list of matrices with markers for each DA region

findMarkersForDAregion <- function(
  cell.idx, da.region.label, obj, ...
){
  n.da <- length(unique(da.region.label)) - 1
  obj@meta.data$da <- 0
  obj@meta.data$da[cell.idx] <- da.region.label
  obj <- SetAllIdent(obj, id = "da")
  
  da.markers <- list()
  for(ii in 1:n.da){
    da.markers[[ii]] <- FindMarkers(obj, ident.1 = ii, ...)
    da.markers[[ii]]$pct.diff <- da.markers[[ii]]$pct.1 - da.markers[[ii]]$pct.2
  }
  
  return(da.markers)
}





##=======================================================##
## Other functions

# calculate knn.diff.ratio for each cell
daPerCell <- function(
  X, cell.labels, labels.1, labels.2, k.vector
){
  knn.out <- nn2(data = X, query = X, k = max(k.vector))
  knn.graph <- knn.out$nn.idx
  
  n.cells <- length(cell.labels)
  knn.diff.ratio <- matrix(0, nrow = n.cells, ncol = length(k.vector))
  colnames(knn.diff.ratio) <- as.character(k.vector)
  
  labels.1 <- labels.1[labels.1 %in% cell.labels]
  labels.2 <- labels.2[labels.2 %in% cell.labels]
  
  cell.label.name <- sort(unique(cell.labels))
  cell.label.tab <- table(factor(cell.labels, levels = cell.label.name))
  
  # count labels
  for(ii in 1:n.cells){
    for(kk in k.vector){
      i.kk.label <- cell.labels[knn.graph[ii,1:kk]]
      i.kk.label.ratio <- table(factor(i.kk.label, levels = cell.label.name)) / cell.label.tab
      knn.diff.ratio[ii,as.character(kk)] <- (mean(i.kk.label.ratio[labels.2]) - mean(i.kk.label.ratio[labels.1])) / 
        sum(i.kk.label.ratio)
      #  (mean(i.kk.label.ratio[labels.2]) + mean(i.kk.label.ratio[labels.1]))
    }
  }
  
  return(knn.diff.ratio)
}



# remove neutral cells by cluster knn.diff.ratio
removeNeutralCells <- function(X, ratio = 0.2, i.max = 10, keep.info = F){
  n <- nrow(X)
  idx.out <- c(1:n)
  remove.info <- rep(0,n)
  
  for(ii in 1:i.max){
    if((length(idx.out)/n) <= ratio & ii > 1){
      break
    }
    
    # cluster into 3 groups
    X.clust <- kmeans(X[idx.out,], centers = 3)
    X.mean.by.clust <- by(
      rowMeans(X[idx.out,]), INDICES = X.clust$cluster, FUN = mean
    )
    to.remove <- names(X.mean.by.clust)[order(X.mean.by.clust)[2]]
    # to.remove <- names(X.mean.by.clust)[which.min(abs(X.mean.by.clust))]
    
    remove.info[idx.out][which(X.clust$cluster %in% to.remove)] <- rep(ii, length(to.remove))
    
    idx.hat <- which(X.clust$cluster %in% setdiff(c(1:3),to.remove))
    idx.out <- idx.out[idx.hat]
  }
  
  if(keep.info){
    return(list(idx.out, remove.info))
  } else {
    return(idx.out)
  }
}



# get DA regions with tclust
runtclust <- function(X, cell.idx, k, alpha, restr.fact = 50){
  X.tclust.res <- tclust(
    x = X[cell.idx,], k = k, alpha = alpha, restr.fact = restr.fact
  )
  
  return(X.tclust.res$cluster)
}



# calculate score
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



# plot da site
plotDAsite <- function(X, site.list, size = 0.5, cols = NULL){
  colnames(X) <- c("Dim1","Dim2")
  
  site.label <- rep(0, nrow(X))
  for(ii in 1:length(site.list)){
    site.label[site.list[[ii]]] <- ii
  }
  
  myggplot <- ggplot() + theme_classic() + 
    geom_point(data = data.frame(X), aes(Dim1, Dim2), col = "gray", size = size) + 
    geom_point(
      data = data.frame(X[unlist(site.list),]), 
      aes(Dim1, Dim2, col = as.factor(site.label[unlist(site.list)])), 
      size = size
    ) + theme(legend.position = "none")
  
  if(!is.null(cols) & length(cols) == length(site.list)){
    myggplot <- myggplot + 
      scale_color_manual(values = cols)
  }
  
  return(myggplot)
}



plotCellLabel <- function(X, label, cell.col = NULL, size = 0.5, do.label = T, return.plot = F){
  # Add colnames for X
  colnames(X) <- c("Dim1","Dim2")
  
  # Plot cells with labels
  myggplot <- ggplot() + theme_cowplot() +
    geom_point(data = data.frame(Dim1 = X[,1], Dim2 = X[,2], Group = label, stringsAsFactors = F), 
               aes(x = Dim1, y = Dim2, col = Group), size = size) + 
    guides(colour = guide_legend(override.aes = list(size=3), title = NULL))
  
  # Change cell color
  if(!is.null(cell.col)){
    myggplot <- myggplot + scale_color_manual(values = cell.col)
  }
  
  # Label text
  if(do.label){
    mylabels <- unique(label)
    labeldim1 <- by(X[,1], INDICES = label, FUN = median)[mylabels]
    labeldim2 <- by(X[,2], INDICES = label, FUN = median)[mylabels]
    
    myggplot <- myggplot + 
      annotate("text", x = labeldim1, y = labeldim2, label = as.character(mylabels))
  }
  
  if(return.plot) {return(myggplot)} else {print(myggplot)}
}

