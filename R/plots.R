#' @import ggplot2
#' @import cowplot
NULL


#' Plot a score for each cell
#'
#' Produce a ggplot object with cells on 2D embedding, colored by a given score for each cell.
#'
#' @param X matrix, 2D embedding of each cell for the plot
#' @param score numeric vector, a single value to color each cell, continuous
#' @param cell.col string vector, color bar to use for "score", defaul c("blue","white","red")
#' @param size numeric, dot size for each cell, default 0.5
#' @param alpha numeric between 0 to 1, dot opacity for each cell, default 1
#'
#' @return a ggplot object
#' @export
#'
plotCellScore <- function(X, score, cell.col = c("blue","white","red"), size = 0.5, alpha = 1){
  # Add colnames for X
  colnames(X) <- c("Dim1","Dim2")

  # Plot cells with labels
  myggplot <- ggplot() + theme_cowplot() +
    geom_point(data = data.frame(Dim1 = X[,1], Dim2 = X[,2], Score = score),
               aes(x = Dim1, y = Dim2, col = Score), size = size, alpha = alpha) +
    scale_color_gradientn(colours = cell.col)

  return(myggplot)
}



#' Plot da site
#'
#' @param X matrix, 2D embedding of each cell for the plot
#' @param site.list list, a list of cell indices for each site to plot
#' @param size numeric, dot size for each cell, default 0.5
#' @param cols string vector, color bar to use for each site, default ggplot default
#'
#' @export
#'
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



#' Plot cell labels
#'
#' Produce a ggplot object with cells on 2D embedding, colored by given labels of each cell.
#'
#' @param X matrix, 2D embedding of each cell for the plot
#' @param label vector, label for each cell
#' @param cell.col string vector, color bar to use for cell labels, default ggplot default
#' @param size numeric, dot size for each cell, default 0.5
#' @param alpha numeric between 0 to 1, dot opacity for each cell, default 1
#' @param do.label a logical value indicating whether to add text to mark each cell label
#' @param label.size numeric, size of text labels, default 4
#'
#' @return a ggplot object
#' @export
#'
plotCellLabel <- function(X, label, cell.col = NULL, size = 0.5, alpha = 1, do.label = T, label.size = 4){
  # Add colnames for X
  colnames(X) <- c("Dim1","Dim2")

  # Plot cells with labels
  myggplot <- ggplot() + theme_cowplot() +
    geom_point(data = data.frame(Dim1 = X[,1], Dim2 = X[,2], Group = label, stringsAsFactors = F),
               aes(x = Dim1, y = Dim2, col = Group), size = size, alpha = alpha) +
    guides(colour = guide_legend(override.aes = list(size=3,alpha=1), title = NULL))

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
      annotate("text", x = labeldim1, y = labeldim2, label = as.character(mylabels), size = label.size)
  }

  return(myggplot)
}

