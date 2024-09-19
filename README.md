# PPRforEstimatingAutoregressiveOperator

R package for estimating autoregressive operator using Principal Predictive Components (PPC) according to [Kargin-Onatski algorithm](https://core.ac.uk/download/pdf/82625156.pdf), in order to compute one step ahead prediction for time series. Assumption on data: Functional AutoRegressive process.

# Installation

To install the package:
~~~
library(Rcpp)
library(RcppEigen)
library(devtools)
devtools::install_github("AndreaEnricoFranzoni/PPCforAutoregressiveOperator", force = TRUE)
~~~
If MacOS is used, having Fortran installed is mandatory. In case of error during the installation, follow the instructions in this [link](https://cran.r-project.org/bin/macosx/tools/).


# Usage

The regularization parameter is `alpha`, a positive real number.

The number of PPCs retained is `k`.

Example of usage:
~~~
PPCKO::PPC_KO( X = data,
               id_CV = "CV_alpha",
               id_p_for_k = "Yes",
               threshold_k = 0.95,
               alpha = 0.75,
               n_disc = 100,
               k = 0,
               alpha_min = 0.00000000001,
               alpha_max = 10,
               id_p_imposed = "No",
               id_rem_nan = "EMA"
              )
~~~
`X` is the training dataset: each row (total: m) is the evaaluation of the functional element at a specific instant in a point of its doamin, each column (n) is a time instant.

`id_CV` is a string that indicates which algorithm has to be used: can assume values `NoCV` (`alpha` given, `k` can be imposed as parameter (`k` in range {1,...,m}) or searched imposing the variance explained by the PPCs (`k`=0)), `CV_alpha` (the optimal alpha is searched using CV: as training set is used a subset of the time instants, and as validation set is used the next time instants. The mse of the prediction is evaluated. The window is moved up until having the last time instant as validation set; the average of the mse is computed for each value of the parameter into [`alpha_min`,`alpha_max`], discretized using `n_disc` subintervals), `CV_k` (the optimal value of `k` is searched using the same strategy used for `alpha`. `alpha` in this case is given) and `CV` (where a brute force is used to find the best pair `alpha`-`k`). For `alpha`, the algorithm looks for its magnitude.

`id_p_for_k` can assume values `Yes` or `No`: in the first case, the variance explained by the PPCs is already retained from the regularized covariance, while in the latter is done only on the estimate of the `phi` operator.

`threshold_k` indicates how much variance is retained (used in `NoCV` and `CV_alpha` versions).

`id_p_imposed` can assume values `Yes` or `No`: in the first case, `k` has to be passed, and only this components are retained all along the algorithm, while in the latter case, only the `k` components will be retained only after the estimate of `phi`.

`id_rem_nan` indicates which strategy is used to remove eventually NaNs values: `EMA` (default, exponential moving average), `WMA` (weighted moving average), `SMA` (simple moving average), `MR` (mean replacement) and `ZR` (replacement with `0`).

# Testing
In the folder `tests`, the implementation of KO is used to try to reconstruct the result obtained by [Didericksen, Kokoszka & Zhang](https://www.semanticscholar.org/paper/Empirical-properties-of-forecasts-with-the-model-Didericksen-Kokoszka/c1fae9f292c2b42beffe4e4146a2bf9ca005f060), trying to demonstrate the goodness of KO.
