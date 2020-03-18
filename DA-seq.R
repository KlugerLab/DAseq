library(RANN)
library(reticulate)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(scales)
library(tclust)
# library(Seurat)



## Step1: select cells with DA neighborhood

#' @param X size N-by-p matrix, input merged dataset of interest after dimension reduction
#' @param cell.labels size N vector, labels for each input cell
#' @param labels.1 vector, label name(s) that represent condition 1
#' @param labels.2 vector, label name(s) that represent condition 2
#' @param k.vector vector, k values to create the score vector
#' @param k.folds integer, number of data splits used in the neural network, default 10
#' @param n.runs integer, number of times to run the neural network to get the predictions, default 10
#' @param pred.thres length-2 vector, top and bottom threshold on the predictions from the neural network, default c(0.05,0.95)
#' @param do.plot a logical value to indicate whether to return ggplot objects showing the results, default True
#' @param plot.embedding size N-by-2 matrix, 2D embedding for the cells
#' @param size cell size to use in the plot, default 0.5
#' @param python.use character string, the Python to use, default "/usr/bin/python"
#' @param source.code character string, the neural network source code, default "./DA_nn.py"
#' @param GPU which GPU to use, default '', using CPU
#' 
#' @return a list of results
#'         da.ratio: score vector for each cell
#'         da.pred: (mean) prediction from the neural network
#'         da.cell.idx: cell index with the most DA neighborhood
#'         pred.plot: ggplot object showing the predictions of the neural network on plot.embedding
#'         da.cells.plot: ggplot object highlighting cells of da.cell.idx on plot.embedding

getDAcells <- function(
  X, cell.labels, labels.1, labels.2, k.vector,
  k.folds = 10, n.runs = 10, pred.thres = c(0.05,0.95),
  do.plot = T, plot.embedding = NULL, size = 0.5, 
  python.use = "/usr/bin/python", source.code = "./DA_logit.py", GPU = ""
){
  
  # get DA score vector for each cell
  cat("Calculating DA score vector.\n")
  X.knn.ratio <- daPerCell(
    X = X, 
    cell.labels = cell.labels, 
    labels.1 = labels.1, 
    labels.2 = labels.2, 
    k.vector = k.vector
  )
  
  
  # prepare data for python script
  use_python(python = python.use, required = T)
  
  binary.labels <- cell.labels
  binary.labels[cell.labels %in% labels.1] <- 0.0
  binary.labels[cell.labels %in% labels.2] <- 1.0
  
  binary.labels_py <- r_to_py(as.matrix(binary.labels))
  X.knn.ratio_py <- r_to_py(as.matrix(X.knn.ratio))
  
  
  # set GPU device
  py_run_string(paste("os.environ['CUDA_VISIBLE_DEVICES'] = '", GPU, "'", sep = ""))
  
  # get neural network predictions for each cell
  cat("Running neural network classification.\n")
  source_python(file = source.code)
  # py_run_string(paste("epochs = ", epochs, sep = ""))
  py_run_string(paste("k_folds = ", k.folds, sep = ""))
  
  
  # check n.runs
  if(n.runs == 1){
    X.pred <- k_fold_predict_linear(X.knn.ratio_py, binary.labels_py, py$k_folds)
    # X.pred <- logi_predict(X.knn.ratio_py, binary.labels_py, py$epochs)
    # if(linear){
    #   X.pred <- k_fold_predict_linear(X.knn.ratio_py, binary.labels_py, py$k_folds)
    # } else {
    #   X.pred <- k_fold_predict(X.knn.ratio_py, binary.labels_py, py$k_folds)
    # }
    X.pred <- as.numeric(X.pred)
    X.std <- NULL
  } else {
    cat("Running ", n.runs, " runs.\n", sep = "")
    X.pred.all <- list()
    for(ii in 1:n.runs){
      X.pred.all[[ii]] <- as.numeric(k_fold_predict_linear(X.knn.ratio_py, binary.labels_py, py$k_folds))
      # X.pred.all[[ii]] <- as.numeric(logi_predict(X.knn.ratio_py, binary.labels_py, py$epochs))
      # if(linear){
      #   X.pred.all[[ii]] <- as.numeric(k_fold_predict_linear(X.knn.ratio_py, binary.labels_py, py$k_folds))
      # } else {
      #   X.pred.all[[ii]] <- as.numeric(k_fold_predict(X.knn.ratio_py, binary.labels_py, py$k_folds))
      # }
    }
    X.pred.all <- do.call("rbind", X.pred.all)
    X.pred <- colMeans(X.pred.all)
    X.std <- apply(X.pred.all, 2, sd)
  }
  
  
  # select DA cells
  pred.thres <- sort(pred.thres, decreasing = F)
  X.da.idx <- which(
    X.pred < quantile(X.pred, pred.thres[1]) | X.pred > quantile(X.pred, pred.thres[2])
  )
  
  
  # plot results
  if(do.plot & is.null(plot.embedding)){
    warning("plot.embedding must be provided by user if do.plot = T")
    X.pred.plot <- NULL
    X.da.cells.plot <- NULL
  } else if(do.plot & !is.null(plot.embedding)){
    X.pred.plot <- plotCellScore(
      X = plot.embedding, score = X.pred, size = size
    ) + theme(legend.title = element_blank())
    X.da.cells.plot <- plotDAsite(
      X = plot.embedding, 
      site.list = list(
        which(X.pred < quantile(X.pred, pred.thres[1])),
        which(X.pred > quantile(X.pred, pred.thres[2]))
      ), 
      size = size, cols = c("blue","red")
    )
  } else {
    X.pred.plot <- NULL
    X.da.cells.plot <- NULL
  }
  
  # return result
  return(list(
    da.ratio = X.knn.ratio,
    da.pred = X.pred, 
#    da.std = X.std, 
    da.cell.idx = X.da.idx,
    pred.plot = X.pred.plot,
    da.cells.plot = X.da.cells.plot
  ))
}



## Step1.1: play with threshold

#' @param X output from getDAcells()
#' @param pred.thres length-2 vector, top and bottom threshold on the predictions from the neural network, default c(0.05,0.95)
#' @param do.plot a logical value to indicate whether to return ggplot objects showing the results, default True
#' @param plot.embedding size N-by-2 matrix, 2D embedding for the cells
#' @param size cell size to use in the plot, default 0.5
#' 
#' @return a list of results with updated DA cells

updateDAcells <- function(
  X, pred.thres = c(0.05,0.95), do.plot = T, plot.embedding = NULL, size = 0.5
){
  # select DA cells
  X.pred <- X$da.pred
  pred.thres <- sort(pred.thres, decreasing = F)
  X.da.idx <- which(
    X.pred < quantile(X.pred, pred.thres[1]) | X.pred > quantile(X.pred, pred.thres[2])
  )
  
  # plot results
  if(do.plot & is.null(plot.embedding)){
    warning("plot.embedding must be provided by user if do.plot = T")
    X.pred.plot <- NULL
    X.da.cells.plot <- NULL
  } else if(do.plot & !is.null(plot.embedding)){
    X.pred.plot <- plotCellScore(
      X = plot.embedding, score = X.pred, size = size
    ) + theme(legend.title = element_blank())
    X.da.cells.plot <- plotDAsite(
      X = plot.embedding, 
      site.list = list(
        which(X.pred < quantile(X.pred, pred.thres[1])),
        which(X.pred > quantile(X.pred, pred.thres[2]))
      ), 
      size = size, cols = c("blue","red")
    )
  } else {
    X.pred.plot <- NULL
    X.da.cells.plot <- NULL
  }
  
  # return result
  return(list(
    da.ratio = X$da.ratio,
    da.pred = X.pred, 
    da.cell.idx = X.da.idx,
    pred.plot = X.pred.plot,
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



## Step 3: detect genes that characterize DA regions from Step 2

#' @param X matrix, normalized expression matrix of all cells in the dataset, genes are in rows, rownames must be gene names
#' @param cell.idx result "da.cell.idx" from the output of function getDAcells
#' @param da.region.label result "cluster.res" from the output of function getDAregion
#' @param da.regions.to.run numeric (vector), which DA regions to run the marker finder, default is to run all regions
#' @param lambda numeric, regularization parameter that weights the number of selected genes, a larger lambda leads to fewer genes, default 1
#' @param n.runs integer, number of runs to run the model, default 5
#' @param python.use character string, the Python to use, default "/usr/bin/python"
#' @param source.code character string, the neural network source code, default "./STG_model.py"
#' @param return.model a logical value to indicate whether to return the actual model of STG
#' @param GPU which GPU to use, default '', using CPU
#' 
#' @return a list of results:
#'         da.markers: a list of data frame with markers for each DA region
#'         accuracy: a numeric vector showing mean accuracy for each DA region
#'         model: a list of model for each DA region, each model contains:
#'                model: the model of STG of the final run
#'                features: features used to train the model
#'                selected.features: the selected features of the final run
#'                pred: the linear prediction value for each cell from the model

STGmarkerFinder <- function(
  X, cell.idx, da.region.label,
  da.regions.to.run = NULL, 
  lambda = 1.5, n.runs = 5, return.model = F, 
  python.use = "/usr/bin/python", source.code = "./DA_STG.py", GPU = ""
){
  # set Python
  use_python(python = python.use, required = T)
  
  # set GPU device
  py_run_string(paste("os.environ['CUDA_VISIBLE_DEVICES'] = '", GPU, "'", sep = ""))
  
  source_python(file = source.code)
  
  # turn X into Python format
  X.py <- r_to_py(as.matrix(X))
  
  # get DA regions to run
  n.da <- length(unique(da.region.label)) - 1
  if(is.null(da.regions.to.run)){
    da.regions.to.run <- c(1:n.da)
  }
  
  # create DA label vector
  n.cells <- ncol(X)
  da.label <- rep(0, n.cells)
  da.label[cell.idx] <- da.region.label
  
  
  # run model for each da region
  da.markers <- list()
  da.accr <- vector("numeric")
  da.model <- list()
  for(ii in da.regions.to.run){
    # prepare labels
    da.label.bin <- (da.label == ii)
    da.label.bin.py <- r_to_py(as.matrix(da.label.bin))
    da.label.bin.py <- da.label.bin.py$flatten()
    
    py_run_string(sprintf("num_run = %s; lam = %s", n.runs, lambda))
    stg.out <- STG_FS(X.py, da.label.bin.py, py$num_run, py$lam)
    da.markers[[as.character(ii)]] <- rownames(X)[unique(unlist(stg.out[[1]])) + 1]
    da.accr[[as.character(ii)]] <- mean(stg.out[[2]])
    da.model[[as.character(ii)]] <- list()
    da.model[[as.character(ii)]][["model"]] <- stg.out[[3]]
    da.model[[as.character(ii)]][["features"]] <- rownames(X)[stg.out[[4]] + 1]
    da.model[[as.character(ii)]][["selected.features"]] <- rownames(X)[stg.out[[1]][[n.runs]] + 1]
    da.model[[as.character(ii)]][["pred"]] <- stg.out[[5]][[1]][,2] - stg.out[[5]][[1]][,1]
    da.model[[as.character(ii)]][["alpha"]] <- as.numeric(stg.out[[6]])
    names(da.model[[as.character(ii)]][["alpha"]]) <- da.model[[as.character(ii)]][["features"]]
  }
  
  
  ## Get statistics for each gene
  da.markers.result <- list()
  for(ii in da.regions.to.run){
    da.markers.logfc <- sapply(da.markers[[as.character(ii)]], function(x, x.data1, x.data2){
      log2(mean(x.data1[x,] + 1/100000) / mean(x.data2[x,] + 1/100000))
    }, x.data1 = X[,da.label == ii], x.data2 = X[,da.label != ii])
    da.markers.pval <- sapply(da.markers[[as.character(ii)]], function(x, x.data1, x.data2){
      wilcox.test(x = x.data1[x,], y = x.data2[x,])$p.value
    }, x.data1 = X[,da.label == ii], x.data2 = X[,da.label != ii])
    da.markers.result[[as.character(ii)]] <- data.frame(
      gene = da.markers[[as.character(ii)]],
      avg_logFC = da.markers.logfc,
      p_value = da.markers.pval, 
      stringsAsFactors = F
    )
    
    da.markers.result[[as.character(ii)]] <- da.markers.result[[as.character(ii)]][
      order(da.markers.result[[as.character(ii)]][,"p_value"]),
    ]
  }
  
  # return output
  if(return.model){
    return(list(
      da.markers = da.markers.result,
      accuracy = da.accr,
      model = da.model
    ))
  } else {
    return(list(
      da.markers = da.markers.result,
      accuracy = da.accr
    ))
  }
}


#' @param cell.idx result "da.cell.idx" from the output of function getDAcells
#' @param da.region.label result "cluster.res" from the output of function getDAregion
#' @param obj Seurat object that contain ALL cells in the analysis
#' @param version Seurat version of 'obj', default is '2', support '2' or '3'
#' @param ... parameters for Seurat function FindMarkers()
#' 
#' @return a list of matrices with markers for each DA region

findMarkersForDAregion <- function(
  cell.idx, da.region.label, obj, version = "2", ...
){
  if(!version %in% c("2","3")){
    stop("Parameter version can only be '2' or '3'!")
  }
  
  n.da <- length(unique(da.region.label)) - 1
  obj@meta.data$da <- 0
  obj@meta.data$da[cell.idx] <- da.region.label
  
  if(version == "2"){
    obj <- SetAllIdent(obj, id = "da")
  } else if (version == "3"){
    Idents(obj) <- obj@meta.data$da
  }
  
  
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
      i.kk.ratio1 <- sum(i.kk.label %in% labels.1) / sum(cell.labels %in% labels.1)
      i.kk.ratio2 <- sum(i.kk.label %in% labels.2) / sum(cell.labels %in% labels.2)
      knn.diff.ratio[ii,as.character(kk)] <- (i.kk.ratio2 - i.kk.ratio1) / (i.kk.ratio2 + i.kk.ratio1)
      # i.kk.label.ratio <- table(factor(i.kk.label, levels = cell.label.name)) / cell.label.tab
      # knn.diff.ratio[ii,as.character(kk)] <- 
      #   (mean(i.kk.label.ratio[labels.2]) - mean(i.kk.label.ratio[labels.1])) / 
      #   sum(i.kk.label.ratio)
      #  (mean(i.kk.label.ratio[labels.2]) + mean(i.kk.label.ratio[labels.1]))
    }
  }
  
  return(knn.diff.ratio)
}



# get DA regions with tclust
runtclust <- function(X, cell.idx, k, alpha, restr.fact = 50, ...){
  X.tclust.res <- tclust(
    x = X[cell.idx,], k = k, alpha = alpha, restr.fact = restr.fact, ...
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



# plot a score for each cell
plotCellScore <- function(X, score, cell.col = c("blue","white","red"), size = 0.5){
  # Add colnames for X
  colnames(X) <- c("Dim1","Dim2")
  
  # Plot cells with labels
  myggplot <- ggplot() + theme_cowplot() +
    geom_point(data = data.frame(Dim1 = X[,1], Dim2 = X[,2], Score = score), 
               aes(x = Dim1, y = Dim2, col = Score), size = size) + 
    scale_color_gradientn(colours = cell.col)
  
  return(myggplot)
}



# plot da site
plotDAsite <- function(X, site.list, size = 0.5, cols = NULL){
  colnames(X) <- c("Dim1","Dim2")
  
  site.label <- rep(0, nrow(X))
  for(ii in 1:length(site.list)){
    site.label[site.list[[ii]]] <- ii
  }
  
  myggplot <- ggplot() + theme_cowplot() + 
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



plotCellLabel <- function(X, label, cell.col = NULL, size = 0.5, do.label = T, return.plot = T){
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



# STG in general
runSTG <- function(
  X, X.labels, lambda = 1, n.runs = 5, return.model = F, 
  python.use = "/usr/bin/python", source.code = "./STG_model.py"
){
  # set Python
  use_python(python = python.use, required = T)
  source_python(file = source.code)
  
  X.py <- r_to_py(as.matrix(X))
  
  X.label.bin <- (X.labels == 1)
  X.label.bin.py <- r_to_py(as.matrix(X.label.bin))
  X.label.bin.py <- X.label.bin.py$flatten()
  
  py_run_string(sprintf("num_run = %s; lam = %s", n.runs, lambda))
  stg.out <- STG_FS(X.py, X.label.bin.py, py$num_run, py$lam)
  da.markers <- rownames(X)[unique(unlist(stg.out[[1]])) + 1]
  da.accr <- mean(stg.out[[2]])
  da.model <- list()
  da.model[["model"]] <- stg.out[[3]]
  da.model[["features"]] <- rownames(X)[stg.out[[4]] + 1]
  da.model[["selected.features"]] <- rownames(X)[stg.out[[1]][[n.runs]] + 1]
  da.model[["pred"]] <- stg.out[[5]][[1]][,2] - stg.out[[5]][[1]][,1]
  
  da.markers.logfc <- sapply(da.markers, function(x, x.data1, x.data2){
    log2(mean(x.data1[x,] + 1/100000) / mean(x.data2[x,] + 1/100000))
  }, x.data1 = X[,X.labels == 1], x.data2 = X[,X.labels == 0])
  da.markers.pval <- sapply(da.markers, function(x, x.data1, x.data2){
    wilcox.test(x = x.data1[x,], y = x.data2[x,])$p.value
  }, x.data1 = X[,X.labels == 1], x.data2 = X[,X.labels == 0])
  da.markers.result <- data.frame(
    gene = da.markers,
    avg_logFC = da.markers.logfc,
    p_value = da.markers.pval, 
    stringsAsFactors = F
  )
  
  da.markers.result <- da.markers.result[
    order(da.markers.result[,"p_value"]),
    ]
  
  # return results
  if(return.model){
    return(list(
      da.markers = da.markers.result,
      accuracy = da.accr,
      model = da.model
    ))
  } else {
    return(list(
      da.markers = da.markers.result,
      accuracy = da.accr
    ))
  }
}


