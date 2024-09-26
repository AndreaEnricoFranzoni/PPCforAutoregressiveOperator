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


# Usage

~~~
PPCKO::PPC_KO( X = data,
               id_CV = "CV_alpha",
               threshold_k = 0.95,
               alpha = 0.75,
               k = 0,
               id_rem_nan = "EMA")
~~~
-**`X`**: R matrix of numeric values: training dataset: each row (m) is the evaluation of the functional element at a specific instant in a point of its doamin, each column (n) is a time instant.

-**`id_CV`**: string: indicates which version of KO algorithm is used:
- `NoCV`: algorithm is performed with input parameters; 
- `CV_alpha`: CV is performed to select the best magnitude of `alpha`. `k` is passed as below; 
- `CV_k`: CV is performed to select the best `k`. `alpha` is passed as below;
- `CV`: CV is performed to select the best pair `alpha`-`k`. 

-**`threshold_k`**: double in $(0,1)$: how much explanatory power the PPCs retain.

-**`alpha`**: double in $(0,+\infty)$: regularization parameter.

-**`k`**: integer in $\{1, \dots, m\}$: number of PPCs retained: if `0`, it will be the smallest number that satisfies `threshold_k` explanatory power retained; otherwise, `k` are retained. 

-**`id_rem_nan`**: string that indicates which strategy to remove eventual NaNs is used:
- `MR` (default): NaNs are substituted by the mean of the temporal serie;
- `ZR`: NaNs are substituted by zeros.


# Testing
In the folder `tests`, the implementation of KO is used to try to reconstruct the result obtained by [Didericksen, Kokoszka & Zhang](https://www.semanticscholar.org/paper/Empirical-properties-of-forecasts-with-the-model-Didericksen-Kokoszka/c1fae9f292c2b42beffe4e4146a2bf9ca005f060), trying to demonstrate the goodness of KO.
