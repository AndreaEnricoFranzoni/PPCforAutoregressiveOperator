#' @title data_2d_wrapper_from_list
#' 
#' @name data_2d_wrapper_from_list
#' @description
#' Wrapping a R-list of matrices of doubles into an R matrix of double. Each element of the list represents a time instant,
#' in which the discrete evaluations of the surface is represented by the matrix. Each column of the output matrix represents a time instants,
#' while each rows a discrete evaluation of the surface. The grid is stored column-wise.
#' @param Xt R-list of matrices of doubles. Length of the list: n: number of time instants. Each matrix
#'           needs to have the same dimensions: [dim_x1,dim_x2], where dim_x1 is the number of discrete evaluations of the 
#'           surface along dimension one, dim_x2 the same along dimension two
#' @return an R matrix of double, as described above
#' @examples
#' library(PPCKO)
#' Xt = list(matrix(c(1,2,NaN,3),nrow=2,ncol=2),matrix(c(4,5,NaN,7),nrow=2,ncol=2))
#' data_2d_wrapper_from_list(Xt)
#' # return  [  1,   4]
#' #         [  2,   5]
#' #         [NaN, NaN]
#' #         [  3,   7]
#' @details NaNs values have to be passed anyway, and they are handled by PPCKO algorithm coherently.
#' @export
#' @author Andrea Enrico Franzoni
NULL


#' @title data_2d_wrapper_from_array
#'
#' 
#' @name data_2d_wrapper_from_array
#' @description
#' Wrap an R-array of double into an R matrix of double. Each column of the output matrix represents a time instants,
#' while each rows a discrete evaluation of the surface. The grid is stored column-wise.
#' 
#' @param Xt R-array of doubles, dimensions: [dim_x1,dim_x2,n], where dim_x1 is the number of discrete evaluations of the 
#'           surface along dimension one, dim_x2 the same along dimension two, n the number of time instants
#' @return an R matrix of double, as described above
#' @examples
#' Xt = array(c(1,2,NaN,3,4,5,NaN,7),dim=c(2,2,2))
#' PPCKO::data_2d_wrapper_from_array(Xt)
#' # return  [  1,   4]
#' #         [  2,   5]
#' #         [NaN, NaN]
#' #         [  3,   7]
#' @details
#' NaNs values have to be passed anyway, and they are handled by PPCKO algorithm coherently.
#' @export
#' @author Andrea Enrico Franzoni
NULL