# Principal Predictive Components for Estimating an Autoregressive Operator

**`PPCKO`**  is a C++ based R package for estimating an autoregressive operator. It relays on Principal Predictive Components (PPC) [Kargin-Onatski algorithm](https://core.ac.uk/download/pdf/82625156.pdf) [1](#ref-PPCKO) to compute one-step ahead prediction of time series of curves and surfaces.

Briefly: if data a stationary Functional AutoRegressive of order 1 process (*FAR(1)*), the point predictor assumes the form of:

$\hat{f}_{n+1} = \sum _{i=1}^{k} \langle f_n, b_i \rangle  a_i$

where: 
- k is the number of retained PPCs;
- $f_n$ is the curve/surface at the last available instant;
- $a_i$ are the directions along which the prediction is approximated (directions that have the most predictive power) (predictive loadings);
- $b_i$ weights the prediction along the $i$-th PPC: $\langle f_n, b_i \rangle$ is the projection along it (predictive factors).

> ❗️ **N.B.:** The repository contains only the development of the algorithm through C++ its interface on R using RcppEigen and documentation. More detailed analysis about statistical properties of the predictor and implementation performances can be found [here](https://github.com/AndreaEnricoFranzoni/Functional_time_series_forecasting).



# Prerequisites

R has to be updated at least to 4.0.0 version. If Windows is used, R version has to be at least 4.4.0.

On R console:
~~~
library(devtools)
~~~

Or, alternatively, if not installed:
~~~
install.packages("devtools")
library(devtools)
~~~

**`PPCKO`** depends also on having Fortran, Lapack, BLAS and, if parallel version is used, OpenMP installed.
Depending on the operative system, the instructions to set up everything can be found [here below](#prerequisites-depending-on-operative-system).



# Installation

To install the package:

- Parallel version (highly recommended)
~~~
devtools::install_github("AndreaEnricoFranzoni/PPCforAutoregressiveOperator")
~~~

- Serial version (not available at the moment)
~~~
devtools::install_github("AndreaEnricoFranzoni/PPCforAutoregressiveOperator@v2.0.0.no_omp") 
~~~

And upload the library in the R environment
~~~
library(PPCKO)
~~~



# Documentation 

~~~
?PPCKO
~~~
opens up the documentation. The links to the documentation of all **`PPCKO`** functions can be found here.



## Tutorials

~~~
?PPCKO
~~~
opens up the **`PPCKO`** documentation. In the section `Examples`, there is one tutorial for unidimensional case and one for bidimensional case. Although they do not run directly from here, can be easily reproduced locally. 



## Test

To automatically test the package installation:
~~~
devtools::test()
~~~



## Data

A time series of functions (synthetic data) can be loaded on the global environment using
~~~
data(data_1d)
~~~

A time series of surfaces (synthetic data) can be loaded on the global environment using
~~~
data(data_2d)
~~~



# Usage: unidimensional domain
## KO algorithm
~~~
PPC_KO( Rcpp::NumericMatrix           X,
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

-**`num_threads`**: number of threads for parallelization through OMP. Default: the max number of available threads in the machine. If the serialized version is used, will be automatically put equal to $1$.

-**`id_rem_nan`**: string that indicates which strategy to remove eventual NaNs is used:
- `MR` (default): NaNs are substituted by the mean of the temporal serie in that point of the domain;
- `ZR`: NaNs are substituted by zeros.

**RETURN**: list containing:

-**`One-step ahead prediction`**: vector containing the one step ahead prediction (for each available point of the domain).

-**`Alpha`**: value of the regularization parameter `alpha` used.

-**`Number of PPCs retained`**: number of PPCs retained (value of `k`).

-**`Scores along PPCs`**: vector containing the scores along the PPCs.

-**`Explanatory power PPCs`**: vector containing the cumulative explanatory power up to PPC $i$-th.

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


  
## KO algorithm hypothesiss check
~~~
KO_check_hps(Rcpp::NumericMatrix X)
~~~
-**`X`**: matrix of numeric (real) values: each row (*m*) is the evaluation of the functional element in a point of its domain, each column (*n*) is a specific time instant (equispaced) at which the evaluation occurs.

**RETURN**: list containing:

-**`Pvalues ADF`**: vector containing the p-values of the pointwise (for each point of the domain for which the evaluation is available) ADF (Augmented Dickey-Fuller) test.



## KO algorithm results visualization
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
## KO algorithm
~~~
PPC_KO_2d(Rcpp::NumericMatrix           X,
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
          Rcpp::Nullable<int>           num_threads      = R_NilValue,
          Rcpp::Nullable<std::string>   id_rem_nan       = R_NilValue)
~~~
Inputs have the same meaning of the unidimensional domain case: here are explained only the ones that differ.

-**`X`**: matrix of numeric (real) values: each row (*m*) is the evaluation of the functional element in a point of its domain, each column (*n*) is a specific time instant (equispaced) at which the evaluation occurs. In this case, originally, each time instants is a matrix, in which each entries contains the evaluation of the functional data in that specific point. Data have to be mapped to obtain an input equal to the unidimensional domain case: a vector in which all the columns are lined up represents a single time instant (*m* is the total number of actual evaluations in the grid): are then put next to each other sequentially. To represent more complex domains: NaNs are put in all the points of the grid that do not belong to the domain, in EACH ONE of the time instants, and they concurr in the count of *m*. In this way, the data reader understands that these points do not belong to the domain. [Below](#usage-utilities-to-map-bidimensional-domain-data), some functions are used to map the data from different storaging strategies to the ones needed from the package main function (to understand in which format data can be stored). 

-**`disc_ev_xi`**: discrete points of the domain for which the evaluations are available along dimension $i$. Default value: equispaced grid between left and right extreme of domain $i$-th. 

-**`num_disc_ev_xi`**: number of discrete evaluations (actual AND NaNs) available for dimension $i$.

-**`left_extreme_xi`**: left extreme of the domain of the data for dimension $i$. 

-**`right_extreme_xi`**: right extreme of the domain of the data for dimension $i$.

**RETURN**: list containing:

-**`One-step ahead prediction`**: matrix containing the one step ahead prediction (for each available point of the domain). NaNs are put in the points of the rectangle that do not belong to the domain.

-**`Alpha`**: values of the regularization parameter `alpha` used.

-**`Number of PPCs retained`**: number of PPCs retained (value of `k`).

-**`Scores along PPCs`**: vector containing the scores along the PPCs.

-**`Explanatory power PPCs`**: vector containing the cumulative explanatory power up to PPC $i$-th.

-**`Directions of PPCs`**: list in which each element is the direction (discrete evaluated function) along PPC $i$-th.

-**`Weights of PPCs`**: list in which each element is the weight (discrete evaluated function) along PPC $i$-th.

-**`Validation errors`**: if `err_ret==1`: if `id_CV==NoCV`: empty vector; if `id_CV==CV_alpha`: vector containing the validation errors for each tried value of `alpha`; if `id_CV==CV_k`: vector containing the validation errors for each value of `k` actually tried; if `id_CV==CV`: matrix containing the validation errors for the pairs `alpha`-`k` (only `k` actually tried, if not put the validation error for the greatest `k` actually tried given that `alpha`). If `err_ret==0`: element not returned.

-**`"Function discrete evaluations points dimi`**: discretization of the domain of the functional data for which evaluations are available along dimension $i$.

-**`Left extreme domain dimi`**: left extreme of the domain of the data for dimension $i$.

-**`Right extreme domain dimi`**: right extreme of the domain of the data for dimension $i$.

-**`f_n`**: evaluation of the functional data in last time instant available.

-**`CV`**: which algorithm has been performed.

-**`Alphas`**: vector of values of `alpha` tried in the CV.

-**`K_s`**: vector of values of `k` tried in the CV.



## KO algorithm hypothesis check
~~~
KO_check_hps_2d( Rcpp::NumericMatrix X,
                 int                 dim_x1,
                 int                 dim_x2 )
~~~
-**`X`**: matrix of numeric (real) values: each row (*m*) is the evaluation of the functional element in a point of its domain, each column (*n*) is a specific time instant (equispaced) at which the evaluation occurs. In this case, originally, each time instants is a matrix, in which each entries contains the evaluation of the functional data in that specific point. Data have to be mapped to obtain an input equal to the unidimensional domain case: a vector in which all the columns are lined up represents a single time instant (*m* is the total number of actual evaluations in the grid): are then put next to each other sequentially. To represent more complex domains: NaNs are put in all the points of the grid that do not belong to the domain, in EACH ONE of the time instants, and they concurr in the count of *m*. In this way, the data reader understands that these points do not belong to the domain. [Below](#usage-utilities-to-map-bidimensional-domain-data), some functions are used to map the data from different storaging strategies to the ones needed from the package main function (to understand in which format data can be stored). 

-**`dim_xi`**: number of discrete evaluations (actual evaluations AND NaNs) along dimension $i$.

**RETURN**: list containing:

-**`Pvalues ADF`**: matrix containing the p-values of the pointwise (for each point of the domain for which the evaluation is available) ADF (Augmented Dickey-Fuller) test.



## KO algorithm results visualization
~~~
KO_show_results_2d( results_ko,
                    hp_ko       = NULL,
                    x1_lab      = "x1",
                    x2_lab      = "x2",
                    z_lab       = "value",
                    true_alphas = FALSE)
~~~
-**`results_ko`**: the output of `PPCKO::PPC_KO_2d`.

-**`hp_ko`**: the output of `PPCKO::KO_check_hps_2d` (optional), done on the same data of `results_ko`.

-**`x1_lab`**: name of the $x1$-axis.

-**`x2_lab`**: name of the $x2$-axis.

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
data_2d_wrapper_from_list(Rcpp::List Xt)
~~~
-**`Xt`**: list in which each element is a matrix containing the evaluation of the functional data. Each matrix is a time instants (temporally equispaced). Put NaNs in each matrix for the points that actually do not belong the data domain (to represent more complex domains).

**RETURN**: matrix with the data to be used as `X` for PPCKO algorithm and PPCKO check hps.

## From array
~~~
data_2d_wrapper_from_array(Rcpp::NumericVector Xt)
~~~
-**`Xt`**: array that stores sequentially the surfaces in time. Element [i,j,k] means evaluation in the point (i,j) of the surface at instant k. Put NaNs in (i,j) at each instant for the points that actually do not belong the data domain (to represent more complex domains).

**RETURN**: matrix with the data to be used as `X` for PPCKO algorithm and PPCKO check hps.



# Prerequisites: depending on operative system

More detailed documentation can be found in [this section](https://cran.r-project.org) of `The R Manuals`.
Although installing **`PPCKO`** shpuld automatically install all the R dependecies, could be worth trying install them manaully if an error occurs.
~~~
install.packages("Rcpp")
install.packages("RcppEigen")
library(Rcpp)
library(RcppEigen)
~~~

## macOS

1. **Fortran**:  unlike Linux and Windows, Fortran has to be installed on macOS: instructions in this [link](https://cran.r-project.org/bin/macosx/tools/). Lapack and BLAS will be consequently installed.

2. **OpenMP**: unlike Linux and Windows, OMP is not installed by default on macOS. Open the terminal and digit the following commands.

- **Homebrew**
  - Check the presence of Homebrew
    ~~~
    brew --version
    ~~~
    - If this command does not give back the version of Homebrew, install it according to the macOS version
    ~~~
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
    ~~~
      1. *M1* or *M2*
      ~~~
      echo 'eval "$(/opt/homebrew/bin/brew shellenv)"' >> ~/.zprofile
      eval "$(/opt/homebrew/bin/brew shellenv)"
      ~~~
      2. *Intel*
      ~~~
      echo 'eval "$(/usr/local/bin/brew shellenv)"' >> ~/.zprofile
      eval "$(/usr/local/bin/brew shellenv)"
      ~~~

- **OMP**
  - Once Homebrew is set, check the presence of OMP
    ~~~
    brew list libomp
    ~~~
  - Install it in case of negative output
    ~~~
    brew install libomp
    ~~~

- **LLVM**
   
  [LLVM](https://llvm.org) is needed to configure clang on macOS in order to use an external OMP version
  - Check its presence
  ~~~
  llvm-config --version
  ~~~
  - Eventually, download it
  ~~~
  brew install llvm
  ~~~

## Windows

1. **Rtools**: can be installed from [here](https://cran.r-project.org/bin/windows/Rtools/). Version 4.4 is needed to install parallel version.

## Linux

Linux installation depends on its distributor. All the commands here reported will refer to Ubuntu and Debian distributors.

Standard developement packages have to be installed. In Ubuntu and Debian, for example, all the packages have been collected into a single one, that is possible to install digiting into the terminal:

   ~~~
sudo apt install r-base-dev
sudo apt install build-essential
   ~~~

## Linux image
Before being able to run the commands explained above for Linux, R has to be downloaded. Afterward, the installation of Fortran, Lapack, BLAS, devtools and its dependecies can be done by digiting into the terminal:
   ~~~
sudo apt-get update
sudo apt install gfortran
sudo apt install liblapack-dev libblas-dev
   ~~~
   ~~~
sudo apt-get install libcurl4-openssl-dev
sudo apt-get install libssl-dev
sudo apt-get install libz-dev
sudo apt-get install -y libcurl4-openssl-dev libssl-dev libxml2-dev
sudo apt install zlib1g-dev
sudo apt install -y libfreetype6-dev libfontconfig1-dev libharfbuzz-dev libcairo2-dev libpango1.0-dev pandoc
   ~~~

                


# Bibliography 
1. <a id="ref-PPCKO"></a> **V. Kargin, A. Onatski**, *Curve forecasting by functional autoregression*, Journal of Multivariate Analysis, 99, 2508-2526, 2008, [paper](https://www.sciencedirect.com/science/article/pii/S0047259X08000961)


2. **Autore**, "Titolo dell'Articolo", *Nome della Rivista*, Volume(Número), pp. xx-yy, Anno.
   - DOI/URL: [DOI o Collegamento](https://doi.org/esempio)
   - Nota: Breve riassunto del contenuto e motivo della citazione.

3. **Autore**, *Titolo del Documento Tecnico*, Istituzione, Anno.
   - URL: [Collegamento alla fonte](https://esempio.com)
   - Nota: Dettagli sulla fonte.
