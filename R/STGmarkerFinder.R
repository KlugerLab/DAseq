#' DA-seq Step 4: STG feature selection
#'
#' Use STG (stochastic gates) to select genes that separate each DA region from the rest of the cells.
#' For a full description of the algorithm, see Y. Yamada, O. Lindenbaum, S. Negahban, and Y. Kluger.
#' Feature selection using stochastic gates. arXiv preprint arXiv:1810.04247, 2018.
#'
#' @param X matrix, normalized expression matrix of all cells in the dataset, genes are in rows,
#' rownames must be gene names
#' @param da.regions output from the function getDAregion()
#' @param da.regions.to.run numeric (vector), which DA regions to run the marker finder,
#' default is to run all regions
#' @param lambda numeric, regularization parameter that weights the number of selected genes,
#' a larger lambda leads to fewer genes, default 1.5
#' @param n.runs integer, number of runs to run the model, default 5
#' @param python.use character string, the Python to use, default "/usr/bin/python"
# @param source.code character string, the neural network source code, default "./STG_model.py"
#' @param return.model a logical value to indicate whether to return the actual model of STG
#' @param GPU which GPU to use, default '', using CPU
#'
#' @import reticulate
#' @importFrom e1071 sigmoid
#'
#' @return a list of results:
#' \describe{
#'   \item{da.markers}{a list of data.frame with markers for each DA region}
#'   \item{accuracy}{a numeric vector showing mean accuracy for each DA region}
#'   \item{model}{a list of model for each DA region, each model contains:
#'     \describe{\item{model}{the model of STG of the final run}
#'     \item{features}{features used to train the model}
#'     \item{selected.features}{the selected features of the final run}
#'     \item{pred}{the linear prediction value for each cell from the model}}
#'   }
#' }
#'
#' @export

STGmarkerFinder <- function(
  X, da.regions,
  da.regions.to.run = NULL,
  lambda = 1.5, n.runs = 5, return.model = T,
  python.use = "/usr/bin/python", GPU = ""
){
  if(!inherits(X, what = "matrix") & !inherits(X, what = "Matrix")){
    X <- as.matrix(X)
  }
  # set Python
  use_python(python = python.use, required = T)

  source_python(file = paste(system.file(package="DAseq"), "DA_STG.py", sep = "/"))

  # set GPU device
  py_run_string(paste("os.environ['CUDA_VISIBLE_DEVICES'] = '", GPU, "'", sep = ""))

  # turn X into Python format
  X.py <- r_to_py(as.matrix(X))

  # get DA regions to run
  n.da <- length(unique(da.regions$da.region.label)) - 1
  if(is.null(da.regions.to.run)){
    da.regions.to.run <- c(1:n.da)
  }

  # create DA label vector
  n.cells <- ncol(X)
  da.label <- da.regions$da.region.label


  # run model for each da region
  da.markers <- list()
  da.accr <- vector("numeric")
  da.model <- list()
  for(ii in da.regions.to.run){
    # prepare labels
    da.label.bin <- as.numeric(da.label == ii)
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
    da.model[[as.character(ii)]][["pred"]] <- e1071::sigmoid(stg.out[[5]][[1]][,2] - stg.out[[5]][[1]][,1])
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



#' STG local markers
#' Run STG to find a set of genes that separate a given DA region from a local subset of cells.
#'
#' @param X matrix, normalized expression matrix of all cells in the dataset, genes are in rows,
#' rownames must be gene names
#' @param da.regions output from the function getDAregion()
#' @param da.region.to.run numeric, which (single) DA region to find local markers for
#' @param cell.label.info vector, cell labeling information to select the local subset of cells to compare
#' with input DA region
#' @param cell.label.to.run cell label(s) to select from cell.label.info that represent
#' the local neiborhood for the input DA region
#' @param lambda numeric, regularization parameter that weights the number of selected genes,
#' a larger lambda leads to fewer genes, default 1.5
#' @param n.runs integer, number of runs to run the model, default 5
#' @param python.use character string, the Python to use, default "/usr/bin/python"
#' @param return.model a logical value to indicate whether to return the actual model of STG
#' @param GPU which GPU to use, default '', using CPU
#'
#' @return a list of results:
#' \describe{
#'   \item{markers}{a list of data.frame with markers for each DA region}
#'   \item{accuracy}{a numeric vector showing mean accuracy for each DA region}
#'   \item{model}{a list of model for each DA region, each model contains:
#'     \describe{\item{model}{the model of STG of the final run}
#'     \item{features}{features used to train the model}
#'     \item{selected.features}{the selected features of the final run}
#'     \item{pred}{the linear prediction value for each cell from the model}}
#'   }
#' }
#'
#' @export
#'
STGlocalMarkers <- function(
  X, da.regions, da.region.to.run, cell.label.info, cell.label.to.run, ...
){
  n.cells <- ncol(X)
  if(length(cell.label.info) != n.cells){
    stop("cell.label.info must have the same length as ncol(X).")
  }

  label.da <- rep("notused", n.cells)
  label.da[cell.label.info %in% cell.label.to.run] <- "local"
  label.da[da.regions$da.region.label == da.region.to.run] <- "da"

  output.markers <- runSTG(
    X = X, X.labels = label.da, label.1 = "da", label.2 = "local", ...
  )
  return(output.markers)
}



#' Run STG
#'
#' Run STG to select a set of genes that separate cells with label.1 from label.2 (other labels)
#'
#' @param X matrix, normalized expression matrix of all cells in the dataset, genes are in rows,
#' rownames must be gene names
#' @param X.labels numeric vector, specify labels for each cell, must be 0 or 1
#' @param label.1 cell label to define markers for
#' @param label.2 second cell label to for comparison, if NULL, use all other labels
#' @param lambda numeric, regularization parameter that weights the number of selected genes,
#' a larger lambda leads to fewer genes, default 1.5
#' @param n.runs integer, number of runs to run the model, default 5
#' @param python.use character string, the Python to use, default "/usr/bin/python"
#' @param return.model a logical value to indicate whether to return the actual model of STG
#' @param GPU which GPU to use, default '', using CPU
#'
#' @import reticulate
#'
#' @return a list of results:
#' \describe{
#'   \item{markers}{a list of data.frame with markers for each DA region}
#'   \item{accuracy}{a numeric vector showing mean accuracy for each DA region}
#'   \item{model}{a list of model for each DA region, each model contains:
#'     \describe{\item{model}{the model of STG of the final run}
#'     \item{cells}{cell names/indices used to train the model}
#'     \item{features}{features used to train the model}
#'     \item{selected.features}{the selected features of the final run}
#'     \item{pred}{the linear prediction value for each cell from the model}}
#'   }
#' }
#'
#' @export
#'
runSTG <- function(
  X, X.labels, label.1, label.2 = NULL,
  lambda = 1.5, n.runs = 5, return.model = T,
  python.use = "/usr/bin/python", GPU = ""
){
  if(!inherits(X, what = "matrix") & !inherits(X, what = "Matrix")){
    X <- as.matrix(X)
  }
  if(is.null(rownames(X))){
    stop("rownames of X must be gene names.")
  }

  # set Python
  use_python(python = python.use, required = T)
  source_python(file = paste(system.file(package="DAseq"), "DA_STG.py", sep = "/"))
  py_run_string(paste("os.environ['CUDA_VISIBLE_DEVICES'] = '", GPU, "'", sep = ""))

  if(!is.null(label.2)){
    X.use <- which(X.labels %in% c(label.1,label.2))
    X <- X[,X.use]
    X.labels <- X.labels[X.use]
  } else {
    X.use <- c(1:ncol(X))
  }
  if(!is.null(colnames(X))){
    X.use <- colnames(X)[X.use]
  }

  X.py <- r_to_py(as.matrix(X))

  X.label.bin <- as.numeric(X.labels == label.1)
  X.label.bin.py <- r_to_py(as.matrix(X.label.bin))
  X.label.bin.py <- X.label.bin.py$flatten()

  py_run_string(sprintf("num_run = %s; lam = %s", n.runs, lambda))
  stg.out <- STG_FS(X.py, X.label.bin.py, py$num_run, py$lam)
  da.markers <- rownames(X)[unique(unlist(stg.out[[1]])) + 1]
  da.accr <- mean(stg.out[[2]])
  da.model <- list()
  da.model[["model"]] <- stg.out[[3]]
  da.model[["cells"]] <- X.use
  da.model[["features"]] <- rownames(X)[stg.out[[4]] + 1]
  da.model[["selected.features"]] <- rownames(X)[stg.out[[1]][[n.runs]] + 1]
  da.model[["pred"]] <- e1071::sigmoid(stg.out[[5]][[1]][,2] - stg.out[[5]][[1]][,1])

  da.markers.logfc <- sapply(da.markers, function(x, x.data1, x.data2){
    log2(mean(x.data1[x,] + 1/100000) / mean(x.data2[x,] + 1/100000))
  }, x.data1 = X[,X.label.bin == 1], x.data2 = X[,X.label.bin == 0])
  da.markers.pval <- sapply(da.markers, function(x, x.data1, x.data2){
    wilcox.test(x = x.data1[x,], y = x.data2[x,])$p.value
  }, x.data1 = X[,X.label.bin == 1], x.data2 = X[,X.label.bin == 0])
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
      markers = da.markers.result,
      accuracy = da.accr,
      model = da.model
    ))
  } else {
    return(list(
      markers = da.markers.result,
      accuracy = da.accr
    ))
  }
}




