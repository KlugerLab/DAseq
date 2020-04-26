#' DA-seq Step 1 & Step 2: select DA cells
#'
#' Step 1: compute a multiscale score measure for each cell of its k-nearest-neighborhood for
#' multiple values of k.
#' Step 2: train a logistic regression classifier based on the multiscale score measure and retain cells
#' that may reside in DA regions.
#'
#' @param X size N-by-p matrix, input merged dataset of interest after dimension reduction
#' @param cell.labels size N vector, labels for each input cell
#' @param labels.1 vector, label name(s) that represent condition 1
#' @param labels.2 vector, label name(s) that represent condition 2
#' @param k.vector vector, k values to create the score vector
#' @param k.folds integer, number of data splits used in the neural network, default 10
#' @param n.runs integer, number of times to run the neural network to get the predictions, default 10
#' @param pred.thres length-2 vector, top and bottom threshold on the predictions from the
#' logistic classification, default c(0.05,0.95)
#' @param do.plot a logical value to indicate whether to return ggplot objects showing the results,
#' default True
#' @param plot.embedding size N-by-2 matrix, 2D embedding for the cells
#' @param size cell size to use in the plot, default 0.5
#' @param python.use character string, the Python to use, default "/usr/bin/python"
# @param source.code character string, the neural network source code, default "./DA_nn.py"
#' @param GPU which GPU to use, default '', using CPU
#'
#' @import reticulate
#' @import RANN
#'
#' @return a list of results
#' \describe{
#'   \item{da.ratio}{score vector for each cell}
#'   \item{da.pred}{(mean) prediction from the neural network}
#'   \item{da.cell.id}{index for DA cells}
#'   \item{pred.plot}{ggplot object showing the predictions of logistic regression on plot.embedding}
#'   \item{da.cells.plot}{ggplot object highlighting cells of da.cell.idx on plot.embedding}
#' }
#'
#' @export

getDAcells <- function(
  X, cell.labels, labels.1, labels.2, k.vector,
  k.folds = 10, n.runs = 10, pred.thres = c(0.05,0.95),
  do.plot = T, plot.embedding = NULL, size = 0.5,
  python.use = "/usr/bin/python", GPU = ""
){
#  cat("Using GPU ", GPU, ".\n", sep = "")

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

  # get neural network predictions for each cell
  cat("Running logistic regression.\n")
  source_python(file = paste(system.file(package="DAseq"), "DA_logit.py", sep = "/"))
  # py_run_string(paste("epochs = ", epochs, sep = ""))
  py_run_string(paste("k_folds = ", k.folds, sep = ""))

  # set GPU device
  py_run_string(paste("os.environ['CUDA_VISIBLE_DEVICES'] = '", GPU, "'", sep = ""))


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



#' Update DA cells
#'
#' Use different threshold to select DA cells based on an output from getDAcells()
#'
#' @param X output from getDAcells()
#' @param pred.thres length-2 vector, top and bottom threshold on the predictions from the
#' logistic classification, default c(0.05,0.95)
#' @param do.plot a logical value to indicate whether to return ggplot objects showing the results,
#' default True
#' @param plot.embedding size N-by-2 matrix, 2D embedding for the cells
#' @param size cell size to use in the plot, default 0.5
#'
#' @return a list of results with updated DA cells
#' @export

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



# Calculate multiscale score vector for each cell
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
