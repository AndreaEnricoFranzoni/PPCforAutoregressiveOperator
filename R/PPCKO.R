#'@title Principal Predictive Components for Estimating an Autoregressive operator
#'@description PPCKO provides the Principal Predictive Components Kargin-Onatski algorithm
#' for performing one-step ahead prediction of Functional time series, forecasting curves
#' and surfaces.
#'
#'@references
#'\itemize{
#'\item V. Kargin and A. Onatski, \emph{"Curve forecasting by functional autoregression"}. Journal of Multivariate Analysis, 99, 2508-2526, 2008.
#'}
#'@seealso
#'\itemize{
#'\item Functional time series of curves:
#'\itemize{
#'\item algorithm: \code{\link{PPC_KO}}
#'\item assumptions check: \code{\link{KO_check_hps}}
#'\item results visualization: \code{\link{KO_show_results}}
#'\item example data: \code{\link{data_1d}}}
#'\item Functional time series of surfaces:
#'\itemize{
#'\item algorithm: \code{\link{PPC_KO_2d}}
#'\item assumptions check: \code{\link{KO_check_hps_2d}}
#'\item results visualization: \code{\link{KO_show_results_2d}}
#'\item example data: \code{\link{data_2d}}
#'\item data wrapper: \code{\link{data_2d_wrapper_from_list}}, \code{\link{data_2d_wrapper_from_array}}}}
#'
#'
#'@examples
#' library(PPCKO)
#' ############# CURVE CASE ##########################
#' 
#' # upload data
#' data(data_1d)
#' # check hp
#' check_hp_1d <- KO_check_hps(data_1d)
#' # forecasting algorithm
#' res_1d <- PPC_KO( X = data_1d,
#'                   id_CV = "NoCV",
#'                   alpha = 0.1,
#'                   k = 3)
#' # visualize results
#' KO_show_results(res_1d,check_hp_1d)
#'
#'
#'
#' ############# SURFACE CASE ##########################
#' 
#' # upload data
#' data(data_2d)
#' data_2d_wrap = data_2d_wrapper_from_list(data_2d)
#' # check hp
#' check_hp_2d <- KO_check_hps_2d(data_2d_wrap,dim_x1=10,dim_x2=10)
#' # forecasting algorithm
#' res_2d <- PPC_KO_2d( X = data_2d_wrap,
#'                      id_CV = "NoCV",
#'                      alpha = 0.1,
#'                      k = 3,
#'                      num_disc_ev_x1 = 10,
#'                      num_disc_ev_x2 = 10)
#' # visualize results
#' KO_show_results_2d(res_2d,check_hp_2d)
#'
#'@docType package
#'@name PPCKO
NULL