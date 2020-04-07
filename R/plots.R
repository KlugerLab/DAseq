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
#'
#' @return a ggplot object
#' @export
#'
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



# Plot da site
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



# Plot cell labels
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

