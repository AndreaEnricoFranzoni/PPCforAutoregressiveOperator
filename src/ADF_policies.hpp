// Copyright (c) 2024 Andrea Enrico Franzoni (andreaenrico.franzoni@gmail.com)
//
// This file is part of PPCKO
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of PPCKO and associated documentation files (the PPCKO software), to deal
// PPCKO without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of PPCKO, and to permit persons to whom PPCKO is
// furnished to do so, subject to the following conditions:
//
// PPCKO IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH PPCKO OR THE USE OR OTHER DEALINGS IN
// PPCKO.

#ifndef ADF_PPC_POLICIES_HPP
#define ADF_PPC_POLICIES_HPP

#include <algorithm>
#include <numeric>
#include <array>

#include "traits_ko.hpp"
#include "ADF_lr.hpp"


/*!
* @file ADF_policies.hpp
* @brief Implement the policies, through functors and template parameter, for considering or not bigger lag orders to perform pointwisely Augmented Dickey-Fuller (ADF) test on the fts. 
* @note The test implementation is as in (https://cran.r-project.org/web/packages/tseries/index.html)
* @author Andrea Enrico Franzoni
*/



/*!
* @struct CaseNoLagOrderADF
* @brief Functor returning the values for which computing ratio for computing ADF-test statistic. No lag orders bigger than one considered.
*/
struct CaseNoLagOrderADF
{ 
  /*!
  * @brief Retains the coefficients for the ADF-test statistic
  * @param x time series (dimensions: (n+1) x 1, n number of time instants)
  * @param z embedded time series (dimensions: n x 1, n number of time instants)
  * @return array containing the number for which the ratio gives back the ADF-test statistic
  */
  std::array<double,2> 
  operator()(const KO_Traits::StoringVector &x, const KO_Traits::StoringMatrix &z) const
  {
    //design matrix: nx3:
    // 3 covariates: intercept, all the time instants (1:n) and time series 
    KO_Traits::StoringMatrix covariates(x.size()-1,2);
    covariates.col(0) = x.head(x.size()-1);
    std::iota(covariates.col(1).begin(),covariates.col(1).end(),static_cast<double>(1.0));
    
    //linear regression 
    //responses: the one step difference
    lr_adf lr(std::move(covariates),std::move(z.col(0)));
    lr.solve();
    
    //returning the values necessary for the test statistics (coefficients on time series and its standard error)
    std::array<double,2> a = {lr.coeff()(1),lr.se_coeff()(1)};
    return a;
  }
};



/*!
* @struct CaseLagOrderADF
* @brief Functor returning the values for which computing ratio for computing ADF-test statistic. Lag orders bigger than one are considered.
*/
struct CaseLagOrderADF
{
  /*!
  * @brief Retains the coefficients for the ADF-test statistic, using lag orders
  * @param x time series (dimensions: (n+1) x 1, n number of time instants)
  * @param z embedded time series (dimensions: (n+1-k_used) x k_used, n number of time instants)
  * @return array containing the number for which the ratio gives back the ADF-test statistic
  */
  std::array<double,2>
  operator()(const KO_Traits::StoringVector &x, const KO_Traits::StoringMatrix &z) const
  {
    //dimension
    auto k = z.cols();                    //dimension of the embedding space (k_used)
    auto n_row_units = x.size() - k;      //the actual number of statistical units in this regression
    
    // k is the lag order
    //design matrix: n - k x k + 1
    // covariates: intercept, time instants from k to n, time series from instant k to n, embedded space dimensions from 2 to k
    KO_Traits::StoringMatrix covariates(n_row_units,k+1);
    covariates.col(0) = x.segment(k-1,n_row_units);
    std::iota(covariates.col(1).begin(),covariates.col(1).end(),static_cast<double>(k));
    covariates.rightCols(k-1) = z.rightCols(k-1);
    
    //linear regression over the first col of the embedded space
    lr_adf lr(std::move(covariates),std::move(z.col(0)));
    lr.solve();
    
    //returning the values necessary for the test statistics (coefficients on time series and its standard error)
    std::array<double,2> a = {lr.coeff()(1),lr.se_coeff()(1)};
    return a;
  }
};

#endif  //ADF_PPC_POLICIES_HPP