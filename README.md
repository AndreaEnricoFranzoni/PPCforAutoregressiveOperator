# Principal Predictive Components for Estimating an Autoregressive Operator

**`PPCKO`**  is a C++ based R package for estimating an autoregressive operator, in the scenario of Functional Time Series (FTS). It relies on Principal Predictive Components (PPCs) [Kargin-Onatski algorithm](#ref-PPCKO) to compute one-step ahead prediction of time series of curves and surfaces. The package provides functions to compute the forecasting, retaining the PPCs, to test the pointwise stationarity of the FTS and to visualize the results.

Assuming data coming from a stationary Functional AutoRegressive of order 1 process (*FAR(1)*)

$f_{n+1} = \rho f_n + \epsilon_{n+1}$


the predictor is shaped as

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



Due to the high number of warnings, can be useful adding as argument of `install_github`
~~~
quiet=TRUE
~~~
to disable them

If problem arises when installing, due to problems in the dependencies, also the argument 
~~~
dependencies = TRUE
~~~
could be useful






# Documentation 

~~~
?PPCKO
~~~
opens up the documentation of R-interfaced functions, if digited on the R console. 

The documentation of all **`PPCKO`** C++ internal code can be found [here](https://andreaenricofranzoni.github.io/PPCforAutoregressiveOperator/).





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



# Prerequisites: depending on operative system

More detailed documentation can be found in [this section](https://cran.r-project.org) of `The R Manuals`.
Although installing **`PPCKO`** should automatically install all the R dependecies, could be worth trying to install them manaully if an error occurs.
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

- Check that the command that opens the graphic window used by the visualization results functions
  ~~~
  quartz()
  ~~~
  works on R, by digiting it on R console. If not, look [here](https://www.xquartz.org).

## Windows

- **Rtools**: can be installed from [here](https://cran.r-project.org/bin/windows/Rtools/). Version 4.4 is needed to install parallel version.

- Check that the command that opens the graphic window used by the visualization results functions
  ~~~
  windows()
  ~~~
  works on R, by digiting it on R console.

## Linux

- Linux installation depends on its distributor. All the commands here reported will refer to Ubuntu and Debian distributors. Standard developement packages have to be installed.   In Ubuntu and Debian, for example, all the packages have been collected into a single one, that is possible to install digiting into the terminal:

   ~~~
  sudo apt install r-base-dev
  sudo apt install build-essential
   ~~~

- Check that the command that opens the graphic window used by the visualization results functions
  ~~~
  x11()
  ~~~
  works on R, by digiting it on R console.

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
1. <a id="ref-PPCKO"></a> **Kargin, V., Onatski, A.**, `Curve forecasting by functional autoregression`, *Journal of Multivariate Analysis*, 99, 2508-2526, 2008, https://www.sciencedirect.com/science/article/pii/S0047259X08000961

2. <a id="ref-Rcpp"></a> **Eddelbuettel D., Francois R., Allaire J., Ushey K., Kou Q., Russell N., Ucar I., Bates D., Chambers J.**, `Rcpp: Seamless R and C++ Integration`, *R package version 1.0.13-1*, 2024, https://cran.r-project.org/web/packages/Rcpp/citation.html

3. <a id="ref-eigen"></a> **Guennebaud G., Jacob B., et Al.**, `Eigen v3.4`, 2021, https://eigen.tuxfamily.org/index.php?title=Main_Page

4. <a id="ref-spectra"></a> **Qui X., et Al.**, `Spectra - C++ Library For Large Scale Eigenvalue Problems`, 2022, https://spectralib.org

5. <a id="ref-mlpack"></a> **Curtin R. R., Edel M., Shrit O., Agrawal S., Basak S., Balamuta J. J., Birmingham R., Dutt K., Eddelbuettel D., Garg R., Jaiswal S., Kaushik A., Kim S., Mukherjee A., Gnana Sai N., Sharma N., Singh Parihar Y., Swain R., Sanderson C.**, `mlpack 4: a fast, header-only C++ machine learning library`, *Journal of Open Source Software*, 8, 82, 5026, 2023, https://www.mlpack.org

6. <a id="ref-armadillo"></a> **Sanderson C., Curtin R. R.**,  `Armadillo: a template-based C++ library for linear algebra`, *Journal of Open Source Software*, 1, 2, 26, 2016, https://arma.sourceforge.net

7. <a id="ref-armadillo2"></a> **Sanderson C., Curtin R. R.**,   `Practical Sparse Matrices in C++ with Hybrid Storage and Template-Based Expression Optimisation `, *Mathematical and Computational Applications*, 24, 3, 2019, https://arxiv.org/abs/1811.08768

8. <a id="ref-ensmallen"></a> **Curtin R. R., Edel M., Ganesh Prabhu R., Basak S., Lou Z., Sanderson C.**, `The ensmallen library for flexible numerical optimization`, *Journal of Machine Learning Research*, 22, 166, 2021, https://ensmallen.org

9. <a id="ref-cereal"></a> **Huebner L., et Al.**, `cereal - A C++11 library for serialization`, 2022, https://github.com/USCiLab/cereal

10. <a id="ref-pacsexamples"></a> **Formaggia L., et Al.**, `EXAMPLES AND EXERCISES FOR AMSC and APSC (PACS) COURSES`, 2024, https://github.com/pacs-course/pacs-examples

11. <a id="ref-tseries"></a> **Trapletti A., Hornik K., LeBaron B.**, `tseries: Time Series Analysis and Computational Finance`, *CRAN repository*, 2024, https://cran.r-project.org/web/packages/tseries/index.html
