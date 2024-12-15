#' Functional time series of curves
#'
#' FAR(1) process generating data.
#' Every row (200) represents one point of the domain (0,1) for which the curve evaluation is available.
#' Every column (100) represents the curve at one time instant. 
#'
#' @format A numeric matrix with 200 rows and 100 cols.
#' @usage data(data_1d)
#' @examples
#' data(data_1d)
#' head(data_1d)
"data_1d"



#' Functional time series of surfaces
#' 
#' FAR(1) process generating data.
#' Every element of the list represents the surface at one time instant.
#' Every element of the list is a numeric matrix: it represents the grid over which the 
#' evaluations of the surface are available. 
#'
#' @format A list with 20 10x10 matrices.
#' @usage data(data_2d)
#' @examples
#' data(data_2d)
#' head(data_2d)
"data_2d"
