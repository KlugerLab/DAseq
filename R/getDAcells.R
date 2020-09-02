#' DA-seq Step 1 & Step 2: select DA cells
#'
#' Step 1: compute a multiscale score measure for each cell of its k-nearest-neighborhood for
#' multiple values of k.
#' Step 2: train a logistic regression classifier based on the multiscale score measure and retain cells
#' that may reside in DA regions.
#'
#' @param X size N-by-p matrix, input merged dataset of interest after dimension reduction.
#' @param cell.labels size N character vector, labels for each input cell
#' @param labels.1 character vector, label name(s) that represent condition 1
#' @param labels.2 character vector, label name(s) that represent condition 2
#' @param k.vector vector, k values to create the score vector
#' @param save.knn a logical value to indicate whether to save computed kNN result, default False
#' @param alpha numeric, elasticnet mixing parameter passed to glmnet(), default 0 (Ridge)
#' @param k.folds integer, number of data splits used in the neural network, default 10
#' @param n.runs integer, number of times to run the neural network to get the predictions, default 10
#' @param pred.thres length-2 vector, top and bottom threshold on the predictions from the
#' logistic classification, default c(0.05,0.95)
#' @param do.plot a logical value to indicate whether to return ggplot objects showing the results,
#' default True
#' @param plot.embedding size N-by-2 matrix, 2D embedding for the cells
#' @param size cell size to use in the plot, default 0.5
#'
#' @import RANN
#' @importFrom caret createFolds
#' @import glmnet
#'
#' @return a list of results
#' \describe{
#'   \item{da.ratio}{score vector for each cell}
#'   \item{da.pred}{(mean) prediction from the logistic regression}
#'   \item{da.up}{index for DA cells more abundant in condition of labels.2}
#'   \item{da.down}{index for DA cells more abundant in condition of labels.1}
#'   \item{pred.plot}{ggplot object showing the predictions of logistic regression on plot.embedding}
#'   \item{da.cells.plot}{ggplot object highlighting cells of da.cell.idx on plot.embedding}
#' }
#'
#' @export

getDAcells <- function(
  X, cell.labels, labels.1, labels.2, k.vector = NULL, save.knn = F,
  alpha = 0, k.folds = 10, n.runs = 10, pred.thres = c(0.05,0.95),
  do.plot = T, plot.embedding = NULL, size = 0.5
){
  if(!inherits(x = X, what = "matrix")){
    cat("Turning X to a matrix.\n")
    X <- as.matrix(X)
  }

  # check label input
  if(!inherits(cell.labels, "character") |
     !inherits(labels.1, "character") | !inherits(labels.2, "character")){
    stop("Input parameters cell.labels, labels.1 and labels.2 must be character")
  }
  if(length(setdiff(cell.labels, c(labels.1, labels.2))) > 0){
    stop("Input parameter cell.labels contain labels not from labels.1 or labels.2")
  }
  n.cells <- length(cell.labels)

  # get DA score vector for each cell
  cat("Calculating DA score vector.\n")
  X.knn.result <- daPerCell(
    X = X,
    cell.labels = cell.labels,
    labels.1 = labels.1,
    labels.2 = labels.2,
    k.vector = k.vector
  )
  X.knn.ratio <- X.knn.result[["ratio"]]
  X.knn.graph <- X.knn.result[["graph"]]

  # GLM
  cat("Running GLM.\n")
  binary.labels <- cell.labels
  binary.labels[cell.labels %in% labels.1] <- 0.0
  binary.labels[cell.labels %in% labels.2] <- 1.0

  X.pred <- runDAlasso(
    X = X.knn.ratio, y = factor(binary.labels),
    k.folds = k.folds, n.runs = n.runs, alpha = alpha
  )

  # select DA cells
  pred.thres <- sort(pred.thres, decreasing = F)
  X.da.up <- which(X.pred > quantile(X.pred, pred.thres[2]))
  X.da.down <- which(X.pred < quantile(X.pred, pred.thres[1]))

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
        X.da.down, X.da.up
      ),
      size = size, cols = c("blue","red")
    )
  } else {
    X.pred.plot <- NULL
    X.da.cells.plot <- NULL
  }

  # return result
  if(save.knn){
    return(list(
      knn.graph = X.knn.graph,
      da.ratio = X.knn.ratio,
      da.pred = X.pred,
      da.up = X.da.up,
      da.down = X.da.down,
      #    da.std = X.std,
      #    da.cell.idx = X.da.idx,
      pred.plot = X.pred.plot,
      da.cells.plot = X.da.cells.plot
    ))
  } else {
    return(list(
      da.ratio = X.knn.ratio,
      da.pred = X.pred,
      da.up = X.da.up,
      da.down = X.da.down,
      #    da.std = X.std,
      #    da.cell.idx = X.da.idx,
      pred.plot = X.pred.plot,
      da.cells.plot = X.da.cells.plot
    ))
  }
}



#' Update DA cells
#'
#' Use different threshold to select DA cells based on an output from getDAcells().
#'
#' @param X output from getDAcells()
#' @param pred.thres length-2 vector, top and bottom threshold on the predictions from the
#' logistic classification, default c(0.05,0.95)
#' @param alpha set this parameter to not NULL to rerun Logistic regression
#' @param do.plot a logical value to indicate whether to return ggplot objects showing the results,
#' default True
#' @param plot.embedding size N-by-2 matrix, 2D embedding for the cells
#' @param size cell size to use in the plot, default 0.5

#'
#' @return a list of results with updated DA cells
#' @export

updateDAcells <- function(
  X, pred.thres = c(0.05,0.95),
  alpha = NULL, k.folds = 10, n.runs = 10,
  cell.labels = NULL, labels.1 = NULL, labels.2 = NULL,
  do.plot = T, plot.embedding = NULL, size = 0.5
){
  n.cells <- length(X$da.pred)

  if(!is.null(alpha) & (is.null(cell.labels) | is.null(labels.1) | is.null(labels.2))){
    stop("To update DA cells with new alpha, please also specify cell.labels, labels.1, labels.2.")
  }
  if(!is.null(alpha) & (!is.null(cell.labels) & !is.null(labels.1) & !is.null(labels.2))){
    # check label input
    if(!inherits(cell.labels, "character") |
       !inherits(labels.1, "character") | !inherits(labels.2, "character")){
      stop("Input parameters cell.labels, labels.1 and labels.2 must be character")
    }
    if(length(setdiff(cell.labels, c(labels.1, labels.2))) > 0){
      stop("Input parameter cell.labels contain labels not from labels.1 or labels.2")
    }

    binary.labels <- cell.labels
    binary.labels[cell.labels %in% labels.1] <- 0.0
    binary.labels[cell.labels %in% labels.2] <- 1.0

    X.pred <- runDAlasso(
      X = X$da.ratio, y = factor(binary.labels),
      k.folds = k.folds, n.runs = n.runs, alpha = alpha
    )
  }
  if(is.null(alpha)){
    X.pred <- X$da.pred
  }

  # select DA cells based on new pred.thres
  pred.thres <- sort(pred.thres, decreasing = F)
  X.da.up <- which(X.pred > quantile(X.pred, pred.thres[2]))
  X.da.down <- which(X.pred < quantile(X.pred, pred.thres[1]))

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
        X.da.down, X.da.up
      ),
      size = size, cols = c("blue","red")
    )
  } else {
    X.pred.plot <- NULL
    X.da.cells.plot <- NULL
  }

  # return result
  X$da.pred <- X.pred
  X$da.up <- X.da.up
  X$da.down <- X.da.down
  X$pred.plot <- X.pred.plot
  X$da.cells.plot <- X.da.cells.plot
  return(X)
}



# Calculate multiscale score vector for each cell
daPerCell <- function(
  X, cell.labels, labels.1, labels.2, k.vector
){
  knn.out <- nn2(data = X, query = X, k = max(k.vector))
  knn.graph <- knn.out$nn.idx

  n.cells <- length(cell.labels)
  n.k <- length(k.vector)
  # knn.diff.ratio <- matrix(0, nrow = n.cells, ncol = length(k.vector))

  labels.1 <- labels.1[labels.1 %in% cell.labels]
  labels.2 <- labels.2[labels.2 %in% cell.labels]
  cell.labels.bin <- cell.labels
  cell.labels.bin[cell.labels %in% labels.1] <- 0
  cell.labels.bin[cell.labels %in% labels.2] <- 1
  cell.labels.bin <- as.numeric(cell.labels.bin)
  n.label.2 <- sum(cell.labels.bin)
  n.label.1 <- n.cells - n.label.2

  # count labels
  knn.diff.ratio <- matrix(unlist(lapply(seq_len(n.cells), function(ii) lapply(k.vector, function(kk){
    i.kk.label <- cell.labels.bin[knn.graph[ii,1:kk]]
    i.kk.ratio1 <- (kk - sum(i.kk.label)) / n.label.1
    i.kk.ratio2 <- sum(i.kk.label) / n.label.2
    return((i.kk.ratio2 - i.kk.ratio1) / (i.kk.ratio2 + i.kk.ratio1))
  }))), nrow = n.cells, byrow = T)
  colnames(knn.diff.ratio) <- as.character(k.vector)

  return(list(
    "graph" = knn.graph,
    "ratio" = knn.diff.ratio
  ))
}


# LASSO regression with CV and multiple runs
runDAlasso <- function(X, y, k.folds = 10, n.runs = 10, alpha = 0){
  #X.data <- data.frame(X, response = as.factor(y))
  n.obs <- length(y)
  X.pred.all <- list()
  for(ii in 1:n.runs){
    set.seed(ii)
    X.pred.all[[ii]] <- rep(0, n.obs)
    X.folds <- createFolds(y = y, k = k.folds)
    for(jj in 1:k.folds){
      X.glm <- cv.glmnet(x = X[-X.folds[[jj]],], y = y[-X.folds[[jj]]], family = "binomial", alpha = alpha)
      X.glm.pred <- predict(
        object = X.glm, newx = X[X.folds[[jj]],], s = "lambda.1se", type = "response"
      )
      X.pred.all[[ii]][X.folds[[jj]]] <- X.glm.pred
    }
  }
  X.pred.all <- do.call("rbind", X.pred.all)
  X.pred <- colMeans(X.pred.all)
  return(X.pred)
}

