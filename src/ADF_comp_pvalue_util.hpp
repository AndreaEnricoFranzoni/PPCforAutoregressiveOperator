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

#ifndef ADF_PPC_PVAL_HPP
#define ADF_PPC_PVAL_HPP

#include <vector>
#include <functional>

#include "interp_func.hpp"


/*!
* @file ADF_comp_pvalue_util.hpp
* @brief Helper function and constant values to evaluate ADF-test p-value from the value of the statistic: how much the statistic has extreme negative values
* @note The test implementation is as in (https://cran.r-project.org/web/packages/tseries/index.html)
* @author Andrea Enrico Franzoni
*/



/*!Critical values for Dickey–Fuller t-distribution*/ 
std::vector<std::vector<double>> table = { {-4.38, -4.15, -4.04, -3.99, -3.98, -3.96},
                                           {-3.95, -3.80, -3.73, -3.69, -3.68, -3.66},
                                           {-3.60, -3.50, -3.45, -3.43, -3.42, -3.41},
                                           {-3.24, -3.18, -3.15, -3.13, -3.13, -3.12},
                                           {-1.14, -1.19, -1.22, -1.23, -1.24, -1.25},
                                           {-0.80, -0.87, -0.90, -0.92, -0.93, -0.94},
                                           {-0.50, -0.58, -0.62, -0.64, -0.65, -0.66},
                                           {-0.15, -0.24, -0.28, -0.31, -0.32, -0.33} };

/*!Time instants of the critical values for Dickey–Fuller t-distribution*/
std::vector<double> tableT = {25.0, 50.0, 100.0, 250.0, 500.0, 100000.0};

/*!Honorable cases for pvalue evaluation*/
std::vector<double> tablep = {0.01, 0.025, 0.05, 0.10, 0.90, 0.95, 0.975, 0.99};


/*!
* @brief Auxiliary function for having benchmark value for p-value evaluation, depending on time series length
* @param n time series time instants
* @return a vector containing some benchmark values to be compared with the statistic
*/
std::vector<double>
tableipl(const double &n)
{
  //functor to interpolate
  interp_func interp_f{tableT};
  
  std::vector<double> tableipl_;
  tableipl_.resize(table.size());
  
  //if n is too big or too low with respect to benchmarck values: use only most extreme ones
  if(n<tableT.front()){std::transform(table.begin(),table.end(),tableipl_.begin(),[](auto el){return el.front();});}
  if(n>tableT.back()) {std::transform(table.begin(),table.end(),tableipl_.begin(),[](auto el){return el.back();});}
  
  //interpolate the number of time instants for retaining benchmark values for p-values for the time series
  std::transform(table.begin(),table.end(),tableipl_.begin(),[&n,&interp_f](auto el){return interp_f(n,el);});
  
  return tableipl_;
}

#endif  //ADF_PPC_PVAL_HPP