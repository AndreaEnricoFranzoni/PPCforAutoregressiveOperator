#' @title PPC_KO
#' @name PPC_KO
#' @description
#' Performs Principal Components Analysis Kargin-Onatski algorithm to compute one-step
#' ahead prediction of Functional time series of curves. 
#' CV is eventually performed taking an initial training set, and as validation set the next time instant. 
#' Then, the training set is shifted incrementally by one instant at each iteration, and so the validation set.
#' The validation error is the average of the estimate of the squared error between prediction and validation set.
#' @param X **`numeric matrix`**. Each row (m) represents a point of the domain in which the curve is evaluated.
#'          Each column (n) represents a time instant.
#' @param id_CV **`string`** (default: **`"NoCV"`**). Which version of PPCKO is performed. 
#'              \itemize{
#'              \item "NoCV": PPCKO is performed with the parameters passed as input;
#'              \item "CV_alpha": cv for regularization parameter is performed;
#'              \item "CV_k": cv for the number of retained PPCs is performed;
#'              \item "CV": cv for both the regularization parameter and the number of retained PPCs is performed.
#'              }
#' @param alpha **`double`** (default: **`0.75`**). Strictly positive. Regularization parameter. Will be ignored in "CV_alpha" and "CV" versions.
#' @param k **`integer`** (default: **`0`**). Between 0 and the number of available discrete evaluations of the curve (m).
#'          Number of retained PPCs. Will be ignored in "CV_k" and "CV" versions. If "NoCV" and "CV_alpha" versions:
#'          if it is 0: the number of PPCs retained is chosen through the level of explanatory power criterion (see next parameter);
#'          if it is greater than 0: the number of PPCs retained is k.
#' @param threshold_ppc **`double`** (default: **`0.95`**). Between 0 and 1. Threshold of explanatory power retained. Will be ignored in "CV_k" and "CV" versions,
#'                      and in "NoCV" and "CV_alpha" if k>0. Otherwise, determines the number of PPCs retained depending on data information.
#' @param alpha_vec **`numeric vector`** (default: **`NULL`**). The input space for the regularization parameter in "CV_alpha" and "CV"
#'                  versions. If NULL: logarithmic scale with increasing exponent from 1e-10 up to 1e11 is the input space.
#' @param k_vec **`integer vector`** (default: **`NULL`**). The input space for the number of retained PPCs in "CV_k" and "CV" versions.
#'              If NULL: input space are the integer from 1 up to m.
#' @param toll **`double`** (default: **`1e-4`**). If doing cv for the number of PPCs, by PPCKO construction, adding a PPC can not improve anymore
#'             the predictor. Tolerance for stopping adding PPCs in the cv, evaluated as toll*trace(covariance).
#' @param disc_ev **`numeric vector`** (default: **`NULL`**). Has to have size m. The point of the domain for which the evaluation is available.
#'                 If NULL: a discrete equally spaced grid with m points is assumed.
#' @param left_extreme **`double`** (default: **`0`**). Left extreme of the domain of the functional object.
#' @param right_extreme **`double`** (default: **`1`**). Right extreme of the domain of the functional object.
#' @param min_size_ts **`integer`** (default: **`NULL`**). Between 2 and max_size_ts. The dimension (number of time instants) of the first training set.
#'                    If NULL: is half of n if n even, ceil of half of n if n odd.
#' @param max_size_ts **`integer`** (default: **`NULL`**). Between min_size_ts and n-1. The dimension (number of time instants) of the last training set.
#'                    If NULL: n-1.
#' @param err_ret **`integer`** (default: **`0`**).
#'              \itemize{
#'              \item 0: validation errors are not returned (and not stored during the algorithm);
#'              \item 1: validation errors are returned;
#'              }
#' @param num_threads **`integer`** (default: **`NULL`**). Number of threads for going parallel multithreading.
#'                    Using 1 is equivalent to run the algorithm sequentially (not recommended if doing cv).
#'                    If NULL, or a wrong integer is passed, by default the number of threads used will be equal to the maximum number of threads available for the machine.
#' @param id_rem_nan **`string`** (default: **`NULL`**). Strategy for handling NaNs values. The NaNs are substitute
#'                   only if the row contains other value that are not NaNs, if not the evaluation is interpreted as not belonging
#'                   to the domain (strategy for handling more complex domains).         
#'                   \itemize{
#'                   \item "NO": NaNs are not replaced (**N.B.:  DO NOT USE IT**);
#'                   \item "MR": NaNs are replaced by the avergage of the non-NaNs values of the row (default);
#'                   \item "ZR": NaNs are replaced by 0.
#'                   }
#' @return **`list`** whose items are:
#'                   \itemize{
#'                   \item 'One-step ahead prediction': **`numeric vector`**: numeric vector with the predicted curve;
#'                   \item 'Alpha': **`double`**: regularization parameter used;
#'                   \item 'Number of PPCs retained': **`integer`**: number of retained PPCs;
#'                   \item 'Scores along PPCs': **`numeric vector`**: scores along every PPC. Projection of the last instant over the direction of the PPC;
#'                   \item 'Explanatory power PPCs': **`numeric vector`**: the cumulative explanatory power up to the PPC i-th;
#'                   \item 'Directions of PPCs': **`numeric matrix`**: matrix whose columns are the direction of each PPC;
#'                   \item 'Weights of PPCs': **`numeric matrix`**: matrix whose columns are the weights of each PPC;
#'                   \item 'Validation errors': **`numeric vector`** or **`numeric matrix`**: available only if err_ret==1. For "CV_alpha" and "CV_k"
#'                                              is a vector containing the validation errors for every parameter (for number of PPCs, it is truncated 
#'                                              to the number of PPCs actually tested in the cv process). For "CV" is a matrix, for each pair alpha (row) - k (col);
#'                  \item 'Function discrete evaluations points': the points of the domain for which the evaluations are available;
#'                  \item 'Left extreme domain': left extreme domain;
#'                  \item 'Right extreme domain': right extreme domain;
#'                  \item 'f_n': curve at the last instant;
#'                  \item 'CV': which algorithm version has been performed;
#'                  \item 'Alphas': input space for the regularization parameter;
#'                  \item 'K_s': input space for the number of PPCs retained.
#'                   }
#' @details
#' If more complex domains have to represented, a row of NaNs is put for specific points that do not belong 
#' to the domain but is important to report anyway (for example, in the center of the interval the function is not defined)
#' @references
#' - Paper: \href{https://core.ac.uk/download/pdf/82625156.pdf}{Principal Predictive Components Kargin-Onatski algorithm}
#' - Source code: \href{https://github.com/AndreaEnricoFranzoni/PPCforAutoregressiveOperator}{PPCKO implementation}
#' @export
#' @author Andrea Enrico Franzoni
NULL



#' @title PPC_KO_2d
#' @name PPC_KO_2d
#' @description
#' Performs Principal Components Analysis Kargin-Onatski algorithm to compute one-step
#' ahead prediction of Functional time series of surfaces. 
#' CV is eventually performed taking an initial training set, and as validation set the next time instant. 
#' Then, the training set is shifted incrementally by one instant at each iteration, and so the validation set.
#' The validation error is the average of the estimate of the squared error between prediction and validation set.
#' @param X **`numeric matrix`**. Each row (m) represents a point of the domain in which the surface is evaluated.
#'          The surface is represented by a grid, of dimensions (m1,m2), such that their product is m. 
#'          Then, the grid is encapsulated into a vector, column after column. Each column (n) represents a time instant.
#'          Some auxiliary functions ([data_2d_wrapper_from_list], [data_2d_wrapper_from_array]) are available for wrapping data into a coherent data structure for the algorithm.
#' @param id_CV **`string`** (default: **`"NoCV"`**). Which version of PPCKO is performed. 
#'              \itemize{
#'              \item "NoCV": PPCKO is performed with the parameters passed as input;
#'              \item "CV_alpha": cv for regularization parameter is performed;
#'              \item "CV_k": cv for the number of retained PPCs is performed;
#'              \item "CV": cv for both the regularization parameter and the number of retained PPCs is performed.
#'              }
#' @param alpha **`double`** (default: **`0.75`**). Strictly positive. Regularization parameter. Will be ignored in "CV_alpha" and "CV" versions.
#' @param k **`integer`** (default: **`0`**). Between 0 and the number of available discrete evaluations of the curve (m).
#'          Number of retained PPCs. Will be ignored in "CV_k" and "CV" versions. If "NoCV" and "CV_alpha" versions:
#'          if it is 0: the number of PPCs retained is chosen through the level of explanatory power criterion (see next parameter);
#'          if it is greater than 0: the number of PPCs retained is k.
#' @param threshold_ppc **`double`** (default: **`0.95`**). Between 0 and 1. Threshold of explanatory power retained. Will be ignored in "CV_k" and "CV" versions,
#'                      and in "NoCV" and "CV_alpha" if k>0. Otherwise, determines the number of PPCs retained depending on data information.
#' @param alpha_vec **`numeric vector`** (default: **`NULL`**). The input space for the regularization parameter in "CV_alpha" and "CV"
#'                  versions. If NULL: logarithmic scale with increasing exponent from 1e-10 up to 1e11 is the input space.
#' @param k_vec **`integer vector`** (default: **`NULL`**). The input space for the number of retained PPCs in "CV_k" and "CV" versions.
#'              If NULL: input space are the integer from 1 up to m.
#' @param toll **`double`** (default: **`1e-4`**). If doing cv for the number of PPCs, by PPCKO construction, adding a PPC can not improve anymore
#'             the predictor. Tolerance for stopping adding PPCs in the cv, evaluated as toll*trace(covariance).
#' @param disc_ev_x1 **`numeric vector`** (default: **`NULL`**). Has to have size m1. The points of the domain for which the evaluation is available along dimension one.
#'                   If NULL: a discrete equally spaced grid with m1 points is assumed.
#' @param num_disc_ev_x1 **`integer`** (default: **`10`**). The number of discrete evaluations along dimension one (has to be m1). **`IMPORTANT TO PASS IT CORRECTLY`**.
#' @param disc_ev_x2 **`numeric vector`** (default: **`NULL`**). Has to have size m2. The points of the domain for which the evaluation is available along dimension two.
#'                   If NULL: a discrete equally spaced grid with m2 points is assumed.
#' @param num_disc_ev_x2 **`integer`** (default: **`10`**). The number of discrete evaluations along dimension two (has to be m2). **`IMPORTANT TO PASS IT CORRECTLY`**.
#' @param left_extreme_x1 **`double`** (default: **`0`**). Left extreme of the domain of the functional object along dimension one.
#' @param right_extreme_x1 **`double`** (default: **`1`**). Right extreme of the domain of the functional object along dimension one.
#' @param left_extreme_x2 **`double`** (default: **`0`**). Left extreme of the domain of the functional object along dimension two.
#' @param right_extreme_x2 **`double`** (default: **`1`**). Right extreme of the domain of the functional object along dimension two.
#' @param min_size_ts **`integer`** (default: **`NULL`**). Between 2 and max_size_ts. The dimension (number of time instants) of the first training set.
#'                    If NULL: is half of n if n even, ceil of half of n if n odd.
#' @param max_size_ts **`integer`** (default: **`NULL`**). Between min_size_ts and n-1. The dimension (number of time instants) of the last training set.
#'                    If NULL: n-1.
#' @param err_ret **`integer`** (default: **`0`**).
#'              \itemize{
#'              \item 0: validation errors are not returned (and not stored during the algorithm);
#'              \item 1: validation errors are returned;
#'              }
#' @param num_threads **`integer`** (default: **`NULL`**). Number of threads for going parallel multithreading.
#'                    Using 1 is equivalent to run the algorithm sequentially (not recommended if doing cv).
#'                    If NULL, or a wrong integer is passed, by default the number of threads used will be equal to the maximum number of threads available for the machine.
#' @param id_rem_nan **`string`** (default: **`NULL`**). Strategy for handling NaNs values. The NaNs are substitute
#'                   only if the row contains other value that are not NaNs, if not the evaluation is interpreted as not belonging
#'                   to the domain (strategy for handling more complex domains).         
#'                   \itemize{
#'                   \item "NO": NaNs are not replaced (**N.B.:  DO NOT USE IT**);
#'                   \item "MR": NaNs are replaced by the avergage of the non-NaNs values of the row (default);
#'                   \item "ZR": NaNs are replaced by 0.
#'                   }
#' @return **`list`** whose items are:
#'                   \itemize{
#'                   \item 'One-step ahead prediction': **`numeric matrix`**: numeric matrix with the predicted surface;
#'                   \item 'Alpha': **`double`**: regularization parameter used;
#'                   \item 'Number of PPCs retained': **`integer`**: number of retained PPCs;
#'                   \item 'Scores along PPCs': **`numeric vector`**: scores along every PPC. Projection of the last instant over the direction of the PPC;
#'                   \item 'Explanatory power PPCs': **`numeric vector`**: the cumulative explanatory power up to the PPC i-th;
#'                   \item 'Directions of PPCs': **`list of numeric matrix`**: each element of the list is i-th PPC's direction;
#'                   \item 'Weights of PPCs': **`list of numeric matrix`**: each element of the list is i-th PPC's weight;
#'                   \item 'Validation errors': **`numeric vector`** or **`numeric matrix`**: available only if err_ret==1. For "CV_alpha" and "CV_k"
#'                                              is a vector containing the validation errors for every parameter (for number of PPCs, it is truncated 
#'                                              to the number of PPCs actually tested in the cv process). For "CV" is a matrix, for each pair alpha (row) - k (col);
#'                  \item 'Function discrete evaluations points dim1': the points of the domain for which the evaluations are available along dimension one;
#'                  \item 'Left extreme domain dim1': left extreme domain for dimension one;
#'                  \item 'Right extreme domain dim1': right extreme domain for dimension one;
#'                  \item 'Function discrete evaluations points dim2': the points of the domain for which the evaluations are available along dimension two;
#'                  \item 'Left extreme domain dim2': left extreme domain for dimension two;
#'                  \item 'Right extreme domain dim2': right extreme domain for dimension two;
#'                  \item 'f_n': surface at the last instant;
#'                  \item 'CV': which algorithm version has been performed;
#'                  \item 'Alphas': input space for the regularization parameter;
#'                  \item 'K_s': input space for the number of PPCs retained.
#'                   }
#' @details
#' If more complex domains have to represented, a row of NaNs is put for specific points that do not belong 
#' to the domain but is important to report anyway (for example, to represent a complicated geographical area)
#' @seealso [data_2d_wrapper_from_list], [data_2d_wrapper_from_array]
#' @references
#' - Paper: \href{https://core.ac.uk/download/pdf/82625156.pdf}{Principal Predictive Components Kargin-Onatski algorithm}
#' - Source code: \href{https://github.com/AndreaEnricoFranzoni/PPCforAutoregressiveOperator}{PPCKO implementation}
#' @export
#' @author Andrea Enrico Franzoni
NULL



#' @title KO_check_hps
#' @name KO_check_hps
#' @description
#' Evaluate the pointwise Augmented DIckey-Fuller (ADF) test p-values for the available evaluations of the curve.
#' @param Xt **`numeric matrix`**. Each row (m) represents a point of the domain in which the curve is evaluated.
#'          Each column (n) represents a time instant.
#' @return **`list`** whose items are:
#'         \itemize{
#'         \item 'P-values ADF': vector containing the pointwise ADF test p-values.
#'         }
#' @references 
#' Source code: \href{https://github.com/AndreaEnricoFranzoni/PPCforAutoregressiveOperator}{PPCKO implementation}
#' @export
#' @author Andrea Enrico Franzoni
NULL


#' @title KO_check_hps_2d
#' @name KO_check_hps_2d
#' @description
#' Evaluate the pointwise Augmented DIckey-Fuller (ADF) test p-values for the available evaluations of the surface
#' @param X **`numeric matrix`**. Each row (m) represents a point of the domain in which the surface is evaluated.
#'          The surface is represented by a grid, of dimensions (m1,m2), such that their product is m. 
#'          Then, the grid is encapsulated into a vector, column after column. Each column (n) represents a time instant.
#'          Some auxiliary functions ([data_2d_wrapper_from_list], [data_2d_wrapper_from_array]) are available for wrapping data into a coherent data structure for the algorithm.
#' @param dim_x1 **`integer`**. The number of discrete evaluations along dimension one (has to be m1).
#' @param dim_x2 **`integer`**. The number of discrete evaluations along dimension one (has to be m2).
#' @return **`list`** whose items are:
#'         \itemize{
#'         \item 'P-values ADF': matrix containing the pointwise ADF test p-values.
#'         }
#' @references 
#' Source code: \href{https://github.com/AndreaEnricoFranzoni/PPCforAutoregressiveOperator}{PPCKO implementation}
#' @export
#' @author Andrea Enrico Franzoni
NULL


#' @title KO_show_results
#' @name KO_show_results
#' @description
#' Print on the R console regularization parameter used, number of PPCs retained and cumulative explanatory power of the PPCs retained.
#' Plot pointwise p-value (eventually), prediction at the next instant comparing to the last one available,
#' the validation errors (eventually, if saved), direction and weight for each PPC retained.
#' @param results_ko **`list`** as output from [PPC_KO].
#' @param hp_ko **`list`** (default: **`NULL`**) as output from [KO_check_hps]. If NULL: p-values not plotted.
#' @param x_lab **`string`** (default: **`"x"`**). Name of the x-axis in the plots.
#' @param y_lab **`string`** (default: **`"y"`**). Name of the y-axis in the plots.
#' @param true_alphas **`bool`** (default: **`FALSE`**). If plotting validation errors, if FALSE, the values of alpha are put equally spaced. If TRUE, as their real value
#' @return No return. It simply prints and plots.
#' @seealso [PPC_KO], [KO_check_hps]
#' @references 
#' Source code: \href{https://github.com/AndreaEnricoFranzoni/PPCforAutoregressiveOperator}{PPCKO implementation}
#' @export
#' @author Andrea Enrico Franzoni
NULL



#' @title KO_show_results_2d
#' @name KO_show_results_2d
#' @description
#' Print on the R console regularization parameter used, number of PPCs retained and cumulative explanatory power of the PPCs retained.
#' Plot pointwise p-value (eventually), prediction at the next instant comparing to the last one available, the increment,
#' the validation errors (eventually, if saved), direction and weight for each PPC retained.
#' @param results_ko **`list`** as output from [PPC_KO_2d].
#' @param hp_ko **`list`** (default: **`NULL`**) as output from [KO_check_hps_2d]. If NULL: p-values not plotted.
#' @param x1_lab **`string`** (default: **`"x1"`**). Name of the ax of the first dimension in the plots.
#' @param x2_lab **`string`** (default: **`"x2"`**). Name of the ax of the second dimension in the plots.
#' @param z_lab **`string`** (default: **`"value"`**). Name of the z-ax in the plots.
#' @param true_alphas **`bool`** (default: **`FALSE`**). If plotting validation errors, if FALSE, the values of alpha are put equally spaced. If TRUE, as their real value
#' @return No return. It simply prints and plots.
#' @seealso [PPC_KO_2d], [KO_check_hps_2d]
#' @references 
#' Source code: \href{https://github.com/AndreaEnricoFranzoni/PPCforAutoregressiveOperator}{PPCKO implementation}
#' @export
#' @author Andrea Enrico Franzoni
NULL



#' @title data_2d_wrapper_from_list
#' @name data_2d_wrapper_from_list
#' @description
#' Wrapping a R-list of numeric matrices into an numeric matrix as suitable input for [PPC_KO_2d] and [KO_check_hps_2d]. Each element of the list represents a time instant,
#' in which the discrete evaluations of the surface is represented by the matrix. Each column of the output matrix represents a time instants,
#' while each rows a discrete evaluation of the surface.
#' @param Xt **`list of numeric matrices`**. Length of the list: n: number of time instants. Each matrix
#'           needs to have the same dimensions (dim_x1,dim_x2), where dim_x1 is the number of discrete evaluations of the 
#'           surface along dimension one, dim_x2 the same along dimension two.
#' @return **`numeric matrix`**, as described above.
#' @examples
#' library(PPCKO)
#' Xt = list(matrix(c(1,2,NaN,3),nrow=2,ncol=2),matrix(c(4,5,NaN,7),nrow=2,ncol=2))
#' data_2d_wrapper_from_list(Xt)
#' # return  [  1,   4]
#' #         [  2,   5]
#' #         [NaN, NaN]
#' #         [  3,   7]
#' @details NaNs values have to be passed anyway as explained in [PPC_KO_2d].
#' @seealso [PPC_KO_2d], [KO_check_hps_2d]
#' @references 
#' Source code: \href{https://github.com/AndreaEnricoFranzoni/PPCforAutoregressiveOperator}{PPCKO implementation}
#' @export
#' @author Andrea Enrico Franzoni
NULL



#' @title data_2d_wrapper_from_array
#' @name data_2d_wrapper_from_array
#' @description
#' Wrap an numeric array into an numeric matrix as suitable input for [PPC_KO_2d] and [KO_check_hps_2d]. Each column of the output matrix represents a time instants,
#' while each rows a discrete evaluation of the surface. The grid is stored column-wise.
#' @param Xt **`numeric vector`**, dimensions (dim_x1,dim_x2,n), where dim_x1 is the number of discrete evaluations of the 
#'           surface along dimension one, dim_x2 the same along dimension two, n the number of time instants.
#' @return **`numeric matrix`**, as described above.
#' @examples
#' Xt = array(c(1,2,NaN,3,4,5,NaN,7),dim=c(2,2,2))
#' PPCKO::data_2d_wrapper_from_array(Xt)
#' # return  [  1,   4]
#' #         [  2,   5]
#' #         [NaN, NaN]
#' #         [  3,   7]
#' @details NaNs values have to be passed anyway as explained in [PPC_KO_2d].
#' @seealso [PPC_KO_2d], [KO_check_hps_2d]
#' @references 
#' Source code: \href{https://github.com/AndreaEnricoFranzoni/PPCforAutoregressiveOperator}{PPCKO implementation}
#' @export
#' @author Andrea Enrico Franzoni
NULL