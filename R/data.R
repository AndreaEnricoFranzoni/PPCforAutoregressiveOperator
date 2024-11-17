#' Dataset for 1dim domain data
#'
#' A dataset with 200 rows (discrete evaluation of the functional data) and 100 cols (time instants).
#'
#' @format A dataset with 200 rows and 100 cols:
#' \describe{
#'   \item{x}{a numeric vector (1:10)}
#'   \item{y}{a numeric vector with random values}
#' }
#' @usage data(data_1d)
#' @examples
#' data(data_1d)
#' head(data_1d)
"data_1d"



#' Dataset for 2dim domain data
#'
#' A list with 20 10x10 matrices: each matrix are the evaluation of the surface in a specific time instant (number of element of the list).
#'
#' @format A list with 20 10x10 matrices:
#' \describe{
#'   \item{x}{a numeric vector (1:10)}
#'   \item{y}{a numeric vector with random values}
#' }
#' @usage data(data_2d)
#' @examples
#' data(data_2d)
#' head(data_2d)
"data_2d"
