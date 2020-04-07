#' Top 10 PCs of the melanoma dataset
#'
#' A dataset containing the top 10 PCs (principal components) of the melanoma dataset
#'
#' @source \href{https://www.sciencedirect.com/science/article/pii/S0092867418313941}
#' {Sade-Feldman, Moshe, et al. (Cell. 2018)}
"X.melanoma"

#' Cell sample labels of the melanoma dataset
#'
#' A string vector with the length equal to number of cells, indicating sample labels of each cell: which
#' sample each cell comes from
#'
#' @source \href{https://www.sciencedirect.com/science/article/pii/S0092867418313941}
#' {Sade-Feldman, Moshe, et al. (Cell. 2018)}
"X.label.melanoma"

#' Sample label information
#'
#' A dataset containing information of each sample, indicating whether this sample is a "responder" (R) or
#' a "non-responder" (NR)
#'
#' @format a data.frame with 48 rows and 2 columns
#' \describe{
#'   \item{label}{sample label, matching with labels in X.label.melanoma}
#'   \item{condition}{condition of the sample label, either R or NR}
#' }
#'
#' @source \href{https://www.sciencedirect.com/science/article/pii/S0092867418313941}
#' {Sade-Feldman, Moshe, et al. (Cell. 2018)}
"X.label.info"

#' t-SNE embedding of the melanoma dataset
#'
#' A matrix containing 2D t-SNE embedding of the melanoma dataset
#'
#' @source \href{https://www.sciencedirect.com/science/article/pii/S0092867418313941}
#' {Sade-Feldman, Moshe, et al. (Cell. 2018)}
"X.2d.melanoma"

