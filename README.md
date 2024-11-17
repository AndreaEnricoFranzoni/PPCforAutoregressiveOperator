# Principal Predictive Components for Estimating an Autoregressive Operator

R package for estimating an autoregressive operator using Principal Predictive Components (PPC) according to [Kargin-Onatski algorithm](https://core.ac.uk/download/pdf/82625156.pdf), to compute one-step ahead prediction of functional time series. The domain of the functional element can be unidimensional (time serie of functions) or bidimensional (time serie of surfaces).

Assumptions: data a stationary Functional AutoRegressive of order 1 (*FAR(1)*) process. 

The point predictor assumes the form of:

$\hat{f}_{n+1} = \sum _{i=1}^{k} \langle f_n, b_i \rangle  a_i$

where: 
- k is the number of PPCs;
- $f_n$ is the last instant functional element;
- $a_i$ are the directions along which the prediction is approximated (directions that have the most predictive power) (predictive loadings);
- $b_i$ weights the prediction along the $i$-th PPC: $\langle f_n, b_i \rangle$ is the projection along it (predictive factors).



# Prerequisites
~~~
library(Rcpp)
library(RcppEigen)
library(devtools)
~~~
If MacOS is used, having Fortran installed is mandatory. In case of error during the installation, follow the instructions in this [link](https://cran.r-project.org/bin/macosx/tools/).

If others errors occur, could be that the other needed dependencies have not been actually imported: can be done manually from Rstudio. If already installed:
~~~
library(ggplot2)
library(grid)
library(gridExtra)
library(viridis)
library(patchwork)
library(latex2exp)
~~~



# Installation

To install the package:
~~~
devtools::install_github("AndreaEnricoFranzoni/PPCforAutoregressiveOperator", force = TRUE)
~~~



# Test

To automatically test the package:
~~~
devtools::test()
~~~

The repository contains only the development of the algorithm through C++ and the interface with R through Rcpp. More details about examples, tests on the package and statistical properties of the predictor can be found [here](https://github.com/AndreaEnricoFranzoni/Functional_time_series_forecasting).



# Documentation

TODO



# Usage: unidimensional domain

~~~
PPCKO::PPC_KO( Rcpp::NumericMatrix           X,
               std::string                   id_CV         = "NoCV",
               double                        alpha         = 0.75,
               int                           k             = 0, 
               double                        threshold_ppc = 0.95,
               Rcpp::Nullable<NumericVector> alpha_vec     = R_NilValue,
               Rcpp::Nullable<IntegerVector> k_vec         = R_NilValue,
               double                        toll          = 1e-4,
               Rcpp::Nullable<NumericVector> disc_ev       = R_NilValue,
               double                        left_extreme  = 0,
               double                        right_extreme = 1,
               Rcpp::Nullable<int>           min_size_ts   = R_NilValue,
               Rcpp::Nullable<int>           max_size_ts   = R_NilValue,
               int                           err_ret       = 0,
               Rcpp::Nullable<int>           num_threads   = R_NilValue,
               Rcpp::Nullable<std::string>   id_rem_nan    = R_NilValue)
~~~
-**`X`**: matrix of numeric (real) values: each row (*m*) is the evaluation of the functional element in a point of its domain, each column (*n*) is a specific time instant (equispaced) at which the evaluation occurs.

-**`id_CV`**: string: indicates which version of KO algorithm is used:
- `NoCV (default)`: algorithm is performed with fixed input parameters; 
- `CV_alpha`: CV is performed to select the best value of `alpha`;
- `CV_k`: CV is performed to select the best `k`;
- `CV`: CV is performed to select the best pair {`alpha`-`k`}.
  
CV is performed as follows: the first train set will comprehend data from the beginning of the dataset to a specific time instant, while the respective validation set will be the next time instant. Then, for every CV iteration, the window is augmented: another time instant will be added to the previous train set, and validation set will be the next one, up until another specific time instant. For each validation set, the validation error is evaluated taking the squared mean error of the difference between the discrete evaluation of prediction and actual value, to estimate the $L^2$ norm of the loss. The average of these losses will be the validation error for a specific parameter/pair of parameters.

-**`alpha`**: double, $\in (0,+\infty)$: regularization parameter.

-**`k`**: integer, $\in$ { $0, \dots, m$ }: number of PPCs retained: if `0`, it will be the smallest number that satisfies `threshold_k` explanatory power retained; otherwise, `k` PPCs are retained. 

-**`threshold_k`**: double in $(0,1)$: how much explanatory power the PPCs retain: ignored if `k` is passed or selected through CV.

-**`alpha_vec`**: vector of values of `alpha` tried in the CV. Default value: looking for its magnitude from $10^{-10}$ to $10^{11}$.

-**`k_vec`**: vector of values of `k` tried in the CV. Default value: { $\{1, \dots, m\}$ }.

-**`toll`**: if CV for `k`: adding another PPC increases the predictive power. But after having reached a good amount of explanatory predictive power, adding others PPCs could be meaningless: if the absolute difference of two consecutive validation errors is smaller than `trace(cov)`*`toll`: stop adding others PPCs. Default is $10^{-4}$.

-**`disc_ev`**: discrete points of the domain for which the evaluations are available. Default value: equispaced grid between left and right extreme, of length *m*.

-**`left_extreme`**: left extreme of the domain of the data.

-**`right_extreme`**: right extreme of the domain of the data.

-**`min_size_ts`**: the dimension of the first train set, $\in$ { $2, \dots,$ `max_size_ts`}. If not passed: default value is $\lceil (n/2) \rceil$.

-**`max_size_ts`**: the dimension of the last train set, $\in$ {`min_size_ts` $, \dots,n-1$}. If not passed: default value is $n-1$.

-**`err_ret`**: `1`: validation errors stored and returned. `0`: validation errors not stored and not returned.

-**`num_threads`**: number of threads for parallelization through OMP. Default: the max number of available threads in the machine. If used the serialized version, will be automatically put equal to $1$.

-**`id_rem_nan`**: string that indicates which strategy to remove eventual NaNs is used:
- `MR` (default): NaNs are substituted by the mean of the temporal serie in that point of the domain;
- `ZR`: NaNs are substituted by zeros.

**RETURN**: list containing:

-**`One-step ahead prediction`**: vector containing the one step ahead prediction (for each available point of the domain).

-**`Alpha`**: value of the regularization parameter `alpha` used.

-**`Number of PPCs retained`**: number of PPCs retained (value of `k`).

-**`Scores along PPCs`**: vector containing the scores along the PPCs.

-**`Explanatory power PPCs`**: vector containing the cumulative explanatory power up to PPC $i$.

-**`Directions of PPCs`**: matrix in which each column is the direction $a_i$ (discrete evaluated function) along PPC $i$-th.

-**`Weights of PPCs`**: matrix in which each column is the weight $b_i$ (discrete evaluated function) for PPC $i$-th.

-**`Validation errors`**: if `err_ret==1`: if `id_CV==NoCV`: empty vector; if `id_CV==CV_alpha`: vector containing the validation errors for each tried value of `alpha`; if `id_CV==CV_k`: vector containing the validation errors for each value of `k` actually tried; if `id_CV==CV`: matrix containing the validation errors for the pairs `alpha`-`k` (only `k` actually tried, if not put the validation error for the greatest `k` actually tried given that `alpha`). If `err_ret==0`: element not returned.

-**`Function discrete evaluations points`**: points of the domain for which evaluations are available.

-**`Left extreme domain`**: left extreme of the domain.

-**`Right extreme domain`**: right extreme of the domain.

-**`f_n`**: evaluation of the functional data in last time instant available.

-**`CV`**: which algorithm has been performed.

-**`Alphas`**: vector of values of `alpha` tried in the CV.

-**`K_s`**: vector of values of `k` tried in the CV.


  
---
~~~
PPCKO::KO_check_hps(Rcpp::NumericMatrix X)
~~~
-**`X`**: matrix of numeric (real) values: each row (*m*) is the evaluation of the functional element in a point of its domain, each column (*n*) is a specific time instant (equispaced) at which the evaluation occurs.

**RETURN**: list containing:

-**`Pvalues ADF`**: vector containing the p-values of the pointwise (for each point of the domain for which the evaluation is available) ADF (Augmented Dickey-Fuller) test.



---
~~~
KO_show_results( results_ko,
                 hp_ko       = NULL,
                 x_lab       = "x",
                 y_lab       = "y",
                 true_alphas = FALSE)
~~~
-**`results_ko`**: the output of `PPCKO::PPC_KO`.

-**`hp_ko`**: the output of `PPCKO::KO_check_hps` (optional), done on the same data of `results_ko`.

-**`x_lab`**: name of the $x$-axis.

-**`y_lab`**: name of the $y$-axis.

-**`true_alphas`**: if CV for `alpha` is performed and its validation errors returned: `FALSE`: values of `alpha` are plotted equispaced. `TRUE`: values of `alpha` are plotted with their actual values. (For visualization purposes). 

**RETURN**: void function:

-printing regularization parameter used, number of PPCs retained and cumulative explanatory power of the PPCs retained.

-plotting the pointwise p-values of the ADF test (if `hp_ko` given).

-plotting the comparison between (actual) $f_n$ and (predicted) $f_{n+1}$.

-plotting the validation errors (if validation errors retained) depending on the CV performed.

-plotting the direction and the weight for each PPC.




# Usage: bidimensional domain

~~~
PPCKO::PPC_KO_2d(Rcpp::NumericMatrix           X,
                 std::string                   id_CV            = "NoCV",
                 double                        alpha            = 0.75,
                 int                           k                = 0, 
                 double                        threshold_ppc    = 0.95,
                 Rcpp::Nullable<NumericVector> alpha_vec        = R_NilValue,
                 Rcpp::Nullable<IntegerVector> k_vec            = R_NilValue,
                 double                        toll             = 1e-4,
                 Rcpp::Nullable<NumericVector> disc_ev_x1       = R_NilValue,
                 int                           num_disc_ev_x1   = 10,
                 Rcpp::Nullable<NumericVector> disc_ev_x2       = R_NilValue,
                 int                           num_disc_ev_x2   = 10,
                 double                        left_extreme_x1  = 0,
                 double                        right_extreme_x1 = 1,
                 double                        left_extreme_x2  = 0,
                 double                        right_extreme_x2 = 1,
                 Rcpp::Nullable<int>           min_size_ts      = R_NilValue,
                 Rcpp::Nullable<int>           max_size_ts      = R_NilValue,
                 int                           err_ret          = 0,
                 Rcpp::Nullable<std::string>   id_rem_nan       = R_NilValue)
~~~
Inputs have the same meaning of the unidimensional domain case: here are explained only the ones that differ.

-**`X`**: matrix of numeric (real) values: each row (m) is the evaluation of the functional element in a point of its domain, each column (n) is a specific time instant (equispaced) at which the evaluation occurs. In this case, originally, each time instants is a matrix, in which each entries contains the evaluation of the functional data in that specific point. Data have to be mapped to obtain an input equal to the unidimensional domain case: a vector in which all the columns are lined up represents a single time instant (m is the total number of evaluations in the grid): are then put next to each other sequentially. To represent more complex domains: NaNs are put in all the points of the grid that do not belong to the domain. [Below](#usage-utilities-to-map-bidimensional-domain-data), some functions are used to map the data (to understand in which format data can be stored). 

-**`disc_ev_xi`**: discrete points of the domain for which the evaluations are available along dimension $i$. Default value: equispaced grid between left and right extreme of dimension $i$. 

-**`num_disc_ev_xi`**: number of discrete evaluations available for dimension $i$. Their product is equal to m.

-**`left_extreme_xi`**: left extreme of the domain of the data for dimension $i$.

-**`right_extreme_xi`**: right extreme of the domain of the data for dimension $i$.

**RETURN**: list containing:

-**`One-step ahead prediction`**: matrix containing the one step ahead prediction (for each available point of the domain).

-**`Alpha`**: values of the regularization parameter `alpha` used.

-**`Number of PPCs retained`**: number of PPCs retained (value of `k`).

-**`Scores along PPCs`**: vector containing the scores along the PPCs.

-**`Explanatory power PPCs`**: vector containing the cumulative explanatory power up to PPC $i$.

-**`Directions of PPCs`**: list in which each element is the direction (discrete evaluated function) along PPC $i$.

-**`Weights of PPCs`**: list in which each element is the weight (discrete evaluated function) along PPC $i$.

-**`Validation errors`**: if `err_ret==1`: if `id_CV==NoCV`: empty vector; if `id_CV==CV_alpha`: vector containing the validation errors for each tried value of `alpha`; if `id_CV==CV_k`: vector containing the validation errors for each value of `k` actually tried; if `id_CV==CV`: matrix containing the validation errors for the pairs `alpha`-`k` (only `k` actually tried, if not put the validation error for the greatest `k` actually tried given that `alpha`). If `err_ret==0`: element not returned.

-**`"Function discrete evaluations points dimi`**: discretization of the domain of the functional data for which evaluations are available along dimension $i$.

-**`Left extreme domain dimi`**: left extreme of the domain of the data for dimension $i$.

-**`Right extreme domain dimi`**: right extreme of the domain of the data for dimension $i$.

-**`f_n`**: evaluation of the functional data in last time instant available.

-**`CV`**: which algorithm has been performed.

-**`Alphas`**: vector of values of `alpha` tried in the CV.

-**`K_s`**: vector of values of `k` tried in the CV.



---
~~~
PPCKO::KO_check_hps_2d( Rcpp::NumericMatrix X,
                        int dim_x1,
                        int dim_x2 )
~~~
-**`X`**: matrix of numeric (real) values: each row (m) is the evaluation of the functional element in a point of its domain, each column (n) is a specific time instant (equispaced) at which the evaluation occurs. In this case, originally, each time instants is a matrix, in which each entries contains the evaluation of the functional data in that specific point. Data have to be mapped to obtain an input equal to the unidimensional domain case: a vector in which all the columns are lined up represents a single time instant (m is the total number of evaluations in the grid): are then put next to each other sequentially. To represent more complex domains: NaNs are put in all the points of the grid that do not belong to the domain. Below, some functions are used to map the data (to understand in which format data can be stored).

-**`dim_xi`**: number of discrete evaluations (actual evaluations AND NaNs) along dimension $i$.

**RETURN**: list containing:

-**`Pvalues ADF`**: matrix containing the p-values of the pointwise (for each point of the domain for which the evaluation is available) ADF (Augmented Dickey-Fuller) test.



---
~~~
KO_show_results( results_ko,
                 hp_ko       = NULL,
                 x1_lab      = "x1",
                 x2_lab      = "x2",
                 z_lab       = "value",
                 true_alphas = FALSE)
~~~
-**`results_ko`**: the output of `PPCKO::PPC_KO_2d`.

-**`hp_ko`**: the output of `PPCKO::KO_check_hps_2d` (optional), done on the same data of `results_ko`.

-**`x1_lab`**: name of the $x1$-axis.

-**`x2_lab`**: name of the $x1$-axis.

-**`z_lab`**: name of the $z$-axis.

-**`true_alphas`**: if CV for `alpha` is performed and its validation errors returned: `FALSE`: values of `alpha` are plotted equispaced. `TRUE`: values of `alpha` are plotted with their actual values. (For visualization purposes). 

**RETURN**: void function:

-printing regularization parameter used, number of PPCs retained and cumulative explanatory power of the PPCs retained.

-plotting the pointwise p-values of the ADF test (if `hp_ko` given).

-plotting the comparison between (actual) $f_n$ and (predicted) $f_{n+1}$, other than the increment between the last time instant and the next one in each point of the domain.

-plotting the validation errors (if validation errors retained) depending on the CV performed.

-plotting the direction and the weight for each PPC.



# Usage: utilities to map bidimensional domain data
## From list of matrices
~~~
PPCKO::data_2d_wrapper_from_list(Rcpp::List Xt)
~~~
-**`Xt`**: list in which each element is a matrix containing the evaluation of the functional data. Each matrix is a time instants (temporally equispaced). Put NaNs in each matrix for the points that actually do not belong the data domain (to represent more complex domains).

**RETURN**: matrix with the data to be used as `X` for PPCKO algorithm and PPCKO check hps.



# Bibliography 
TODO
