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

#ifndef ADF_PPC_HPP
#define ADF_PPC_HPP

#include <utility>
#include <vector>
#include <algorithm>

#include "traits_ko.hpp"
#include "ADF_comp_stat.hpp"
#include "ADF_comp_pvalue_util.hpp"


/*!
* @file ADF_test.hpp
* @brief Contains a class to perform pointwisely Augmented Dickey-Fuller (ADF) test on the fts. 
* @note The test implementation is as in (https://cran.r-project.org/web/packages/tseries/index.html)
* @author Andrea Enrico Franzoni
*/



/*!
* @class adf
* @brief Template class for performing pointwise ADF for a fts.
* @tparam LAG_policy indicates if lag orders bigger than one has to be taken into account while computing the test statistic
* @details The ADF is performed for all the rows of the matrix storing the fts as: every row indicates a domain point for which an evaluation of 
*          the functional object is available, every column a time instant.
*/
template<class LAG_policy>
class adf
{

private:
  /*!Matrix containing the fts: each row can be seen as a time series*/
  KO_Traits::StoringMatrix m_x; 
  /*!Lag order*/              
  int m_k;                                 
  /*!Lag order actually used for statistic computation: it is 'm_k'+1*/
  int m_k_used;
  /*!Number of time instants of the time series*/
  int m_tot_time_instants;                    
  /*!P-values of the pointwise ADF test*/
  std::vector<double> m_p_values;             
 
public:

  /*!
  * @brief Constructor taking the matrix containing the fts and the lag order
  * @param x matrix containing the fts
  * @param k lag order
  * @details Universal constructor: move semantic used to optimazing handling big size objects
  */
  template<typename STOR_OBJ>
  adf(STOR_OBJ&& x, int k)      //constructor
   : m_x{std::forward<STOR_OBJ>(x)}, m_k(k), m_k_used(k+static_cast<int>(1)) 
   { m_tot_time_instants = m_x.cols(); } 
  
  /*!
  * @brief Getter for the poitwise ADF test p-values
  * @return the private m_p_values
  */
  inline std::vector<double> p_values() const {return m_p_values;};
  
  /*!
  * @brief Computing the one time step difference for the time series
  * @param ts the time series
  * @return an std::vector<double> containing the one time step differences of the time series
  */
  std::vector<double>      one_step_diff(const KO_Traits::StoringVector &ts)        const;

  /*!
  * @brief Embeds the time series into a low-dimensional Euclidean space
  * @param ts time series
  * @param dimension dimension of the Euclidean space in which embedding the time series
  * @return a matrix containing the embedded time series
  * @details Each row of the resulting matrix consists of sequences x[t], x[t-1], …, x[t-dimension+1], where t is the original index of x.
  * @note Implementation taken from (https://rdocumentation.org/packages/stats/versions/3.6.2)
  */
  KO_Traits::StoringMatrix embed(const KO_Traits::StoringVector &ts, int dimension) const;

  /*!
  * @brief Evaluating the ADF-test statistic for a time series
  * @param ts time series
  * @return statistic evaluation
  */
  double                   statistic_eval(const KO_Traits::StoringVector &ts)       const;      
  
  /*!
  * @brief Evaluating thr ADF-test p-value for a time series
  * @param ts time series
  * @param tableipl table containing extreme values for the statistic for p-value computation
  * @param i number of the discrete evaluation corresponding to 'ts'
  * @return ADF-test p-value
  */
  double                   p_value_eval(const KO_Traits::StoringVector &ts, const std::vector<double> &tableipl, int i)  const;   

  /*!
  * @brief computing the ADF-test pointwisely for each point of the domain for the fts
  */
  void                     test();
  
};
  

#include "ADF_test_imp.hpp"

#endif  //ADF_PPC_HPP