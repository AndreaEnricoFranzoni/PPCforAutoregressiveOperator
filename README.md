# PPRforEstimatingAutoregressiveOperator

R package for estimating autoregressive operator using Principal Predictive Components (PPC) according to [Kargin-Onatski algorithm](https://core.ac.uk/download/pdf/82625156.pdf), in order to compute one step ahead prediction for time series. Assumption on data: Functional AutoRegressive of order 1 process.

# Installation

To install the package:
~~~
library(Rcpp)
library(RcppEigen)
library(devtools)
devtools::install_github("AndreaEnricoFranzoni/PPCforAutoregressiveOperator", force = TRUE)
~~~
If MacOS is used, having Fortran installed is mandatory. In case of error during the installation, follow the instructions in this [link](https://cran.r-project.org/bin/macosx/tools/).


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
               int                           err_ret       = 0,
               Rcpp::Nullable<std::string>   id_rem_nan    = R_NilValue)
~~~
-**`X`**: R matrix of numeric values: each row (m) is the evaluation of the functional element in a point of its domain, each column (n) is a specific time instant (equispaced) at which the evaluation occurs.

-**`id_CV`**: string: indicates which version of KO algorithm is used:
- `NoCV (default)`: algorithm is performed with input parameters; 
- `CV_alpha`: CV is performed to select the best value of `alpha`;
- `CV_k`: CV is performed to select the best `k`;
- `CV`: CV is performed to select the best pair `alpha`-`k`.

-**`alpha`**: double in $(0,+\infty)$: regularization parameter.

-**`k`**: integer in $\{1, \dots, m\}$: number of PPCs retained: if `0`, it will be the smallest number that satisfies `threshold_k` explanatory power retained; otherwise, `k` PPCs are retained. 

-**`threshold_k`**: double in $(0,1)$: how much explanatory power the PPCs retain: ignored if `k` is passed or select through CV.

-**`alpha_vec`**: vector of values of `alpha` tried in the CV. Default value: looking for its magnitude from $10^{-10}$ to $10^{11}$.

-**`k_vec`**: vector of values of `k` tried in the CV. Default value: $\{1, \dots, m\}$.

-**`toll`**: if CV for `k`: how much difference in the validation error to stop adding another PPC.

-**`disc_ev`**: discrete points of the domain for which the evaluations are available. Default value: equispaced grid between left and right extreme, of length m.

-**`left_extreme`**: left extreme of the domain of the data.

-**`right_extreme`**: right extreme of the domain of the data.

-**`err_ret`**: `1`: validation errors stored and returned. `0`: validation errors not stored and not returned.

-**`id_rem_nan`**: string that indicates which strategy to remove eventual NaNs is used:
- `MR` (default): NaNs are substituted by the mean of the temporal serie in that point of the domain;
- `ZR`: NaNs are substituted by zeros.
  
---
~~~
PPCKO::KO_check_hps(Rcpp::NumericMatrix X)
~~~
-**`X`**: R matrix of numeric values: each row (m) is the evaluation of the functional element in a point of its domain, each column (n) is a specific time instant (equispaced) at which the evaluation occurs.

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
                 int                           err_ret          = 0,
                 Rcpp::Nullable<std::string>   id_rem_nan       = R_NilValue)
~~~
Inputs have the sane meaning of the unidimensional domain case: here are explained only the ones that differ.

-**`X`**: R matrix of numeric values: each row (m) is the evaluation of the functional element in a point of its domain, each column (n) is a specific time instant (equispaced) at which the evaluation occurs. In this case, originally, each time instants is a matrix, in which each entries contains the evaluation of the functional data in that specific point. Data have to be mapped to obtain an input equal to the unidimensional domain case: a vector in which all the columns are lined up represents a single time instant (m is the total number of evaluations in the grid): are then put next to each other sequentially. To represent more complex domains: NaNs are put in all the points of the grid that do not belong to the domain. Below, some functions are used to map the data (to understand in which format data can be stored).

-**`disc_ev_xi`**: discrete points of the domain for which the evaluations are available along dimension $i$. Default value: equispaced grid between left and right extreme of dimension $i$. 

-**`num_disc_ev_xi`**: number of discrete evaluations available for dimension $i$. Their product is equal to m.

-**`left_extreme_xi`**: left extreme of the domain of the data for dimension $i$.

-**`right_extreme_xi`**: right extreme of the domain of the data for dimension $i$.

---
~~~
PPCKO::KO_check_hps_2d( Rcpp::NumericMatrix X,
                        int dim_x1,
                        int dim_x2 )
~~~
-**`X`**: R matrix of numeric values: each row (m) is the evaluation of the functional element in a point of its domain, each column (n) is a specific time instant (equispaced) at which the evaluation occurs. In this case, originally, each time instants is a matrix, in which each entries contains the evaluation of the functional data in that specific point. Data have to be mapped to obtain an input equal to the unidimensional domain case: a vector in which all the columns are lined up represents a single time instant (m is the total number of evaluations in the grid): are then put next to each other sequentially. To represent more complex domains: NaNs are put in all the points of the grid that do not belong to the domain. Below, some functions are used to map the data (to understand in which format data can be stored).

-**`dim_xi`**: number of discrete evaluations (actual evaluations AND NaNs) along dimension $i$.

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


# Testing
In the folder `tests`, the implementation of KO is used to try to reconstruct the result obtained by [Didericksen, Kokoszka & Zhang](https://www.semanticscholar.org/paper/Empirical-properties-of-forecasts-with-the-model-Didericksen-Kokoszka/c1fae9f292c2b42beffe4e4146a2bf9ca005f060), trying to demonstrate the goodness of KO.
