#' DA-seq Step 4: Seurat marker finder to characterize DA regions
#'
#' Use Seurat FindMarkers() function to identify characteristic genes for DA regions
#'
#' @param object input Seurat object
#' @param da.slot character, variable name that represents DA regions in Seurat meta.data, default "da"
#' @param da.regions.to.run numeric (vector), which DA regions to find markers for,
#' default is to run all regions
#' @param ... parameters passed to Seurat FindMarkers() function
#'
#' @import Seurat
#'
#' @return a list of markers for DA regions with statistics
#'
#' @export
#'
SeuratMarkerFinder <- function(
  object, da.slot = "da", da.regions.to.run = NULL, ...
){
  # get DA regions to run
  n.da <- length(unique(object@meta.data[,da.slot])) - 1
  if(is.null(da.regions.to.run)){
    da.regions.to.run <- c(1:n.da)
  }

  seurat.version <- substr(packageVersion("Seurat"),1,1)
  if(seurat.version == "3" | seurat.version == "4"){
    Idents(object) <- object@meta.data[,da.slot]

    output.markers <- list()
    for(ii in da.regions.to.run){
      output.markers[[as.character(ii)]] <- FindMarkers(
        object, ident.1 = ii, ...
      )
      output.markers[[as.character(ii)]]$pct.diff <- output.markers[[as.character(ii)]]$pct.1 -
        output.markers[[as.character(ii)]]$pct.2
    }
  }

  return(output.markers)
}



#' Add DA slot
#'
#' Add DA region information to the meta.data of a Seurat object
#'
#' @param object input Seurat object
#' @param da.regions output from function getDAregion()
#' @param da.slot character, variable name to put in Seurat meta.data, default "da"
#' @param set.ident a logical value to indicate whether to set Idents of the Seurat object to DA information,
#' default False
#'
#' @import Seurat
#'
#' @return updated Seurat object
#'
#' @export
#'
addDAslot <- function(object, da.regions, da.slot = "da", set.ident = F){
  seurat.version <- substr(packageVersion("Seurat"),1,1)
  if(seurat.version == "3" | seurat.version == "4"){
    object@meta.data[,da.slot] <- NA
    object@meta.data[,da.slot][da.regions$cell.idx] <- da.regions$da.region.label
    if(set.ident){
      Idents(object) <- object@meta.data[,da.slot]
    }
  }

  return(object)
}



#' Find local markers
#'
#' Use Seurat FindMarkers() function to identify genes that distinguish a DA region from its local neighborhood
#'
#' @param object input Seurat object
#' @param da.slot character, variable name that represents DA regions in Seurat meta.data, default "da"
#' @param da.region.to.run numeric, which (single) DA region to find local markers for
#' @param cell.label.slot character, variable name that represents cell labeling information in Seurat meta.data
#'  to combine with DA information
#' @param cell.label.to.run cell label(s) that represent the local neiborhood for the input DA region
#' @param ... parameters passed to Seurat FindMarkers() function
#'
#' @import Seurat
#'
#' @return a data.frame of markers and statistics
#'
#' @export
#'
SeuratLocalMarkers <- function(
  object, da.slot = "da", da.region.to.run,
  cell.label.slot, cell.label.to.run, ...
){
  seurat.version <- substr(packageVersion("Seurat"),1,1)
  if(seurat.version == "3" | seurat.version == "4"){
    object@meta.data$"label_da" <- "notused"
    object@meta.data$"label_da"[object@meta.data[,cell.label.slot] %in% cell.label.to.run] <- "local"
    object@meta.data$"label_da"[object@meta.data[,da.slot] == da.region.to.run] <- "da"
    Idents(object) <- object@meta.data$"label_da"
    output.markers <- FindMarkers(
      object, ident.1 = "da", ident.2 = "local", ...
    )
    output.markers$pct.diff <- output.markers$pct.1 - output.markers$pct.2
  }
  return(output.markers)
}


