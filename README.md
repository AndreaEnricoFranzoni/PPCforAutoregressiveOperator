# PPRforEstimatingAutoregressiveOperator

Creating an R package for estimating autoregressive operator using Principal Predictive Components (PPC) according to Kargin-Onatski algorithm.

# Installation
If you have MacOS, be sure to have Fortran installed. In case of error, follow the instructions in this [link](https://cran.r-project.org/bin/macosx/tools/).

To install the package:
~~~
library(Rcpp)
library(RcppEigen)
library(devtools)
devtools::install_github("AndreaEnricoFranzoni/PPCforAutoregressiveOperator", force = TRUE)
~~~


# Testing
In the folder `tests`, there is a comparison between the C++ based package and R function based, with respect to the version with and without cv.
Follow the comments in order to run it properly and being able to notice the difference in terms yet of goodness of forecasting yet on computational time.
In `Test.R`, the test is done on `CanadianWeather` dataset, on `monthlyPrecip`. The last month is used as test set.
In `Test2.R`, the test is done on some toy examples.
