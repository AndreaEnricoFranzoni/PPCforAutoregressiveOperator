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

#ifndef ADF_PPC_STAT_HPP
#define ADF_PPC_STAT_HPP

#include <array>

/*!
* @file ADF_comp_stat.hpp
* @brief Implement the class to compute Augmented Dickey-Fuller (ADF) test statistic. 
* @note The test implementation is as in (https://cran.r-project.org/web/packages/tseries/index.html)
* @author Andrea Enrico Franzoni
*/



/*!
* @struct adf_stat
* @brief Contains the functor to evaluate the test statistic, according to LAG_policy
* @tparam LAG_policy policy to select how to evaluate the test statistic (if considering lag order bigger than one)
*/
template <class LAG_policy> 
class adf_stat
{
public:
  
  /*!
  * @brief Computing ADF-test statistic
  * @param x time series
  * @param z embedded time series (embedding to be coherent with the policy used)
  * @return ADF-test statistic
  */
  double
  operator()(const KO_Traits::StoringVector &x, const KO_Traits::StoringMatrix &z) const
  {
    //test statistic of ADF, computed according to lag order policy via template
    std::array<double,2> values = lag_order(x,z);
    return values[0]/values[1];
  }
  
private:
  /*!Policy to correctly computing test statistic*/
  LAG_policy lag_order;
};


#endif  //ADF_PPC_STAT_HPP