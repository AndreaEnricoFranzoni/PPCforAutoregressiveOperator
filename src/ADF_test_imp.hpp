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

#include "ADF_test.hpp"


/*!
* @file ADF_test_imp.hpp
* @brief Implement the class to perform pointwisely Augmented Dickey-Fuller (ADF) test on the fts. 
* @note The test implementation is as in (https://cran.r-project.org/web/packages/tseries/index.html)
* @author Andrea Enrico Franzoni
*/



/*!
* @brief Computing the one time step difference for the time series
* @param ts the time series
* @return an std::vector<double> containing the one time step differences of the time series
*/
template<class LAG_policy>
std::vector<double>
adf<LAG_policy>::one_step_diff(const KO_Traits::StoringVector &ts)
const
{
  //compute the differences within adjacent values of the time series
  std::vector<double> diff;
  diff.resize(m_tot_time_instants-1);
  std::adjacent_difference(ts.cbegin()+1,ts.cend(),diff.begin());
  diff[0] = ts[1] - ts[0];
  
  return diff;
}



/*!
* @brief Embeds the time series into a low-dimensional Euclidean space
* @param ts time series
* @param dimension dimension of the Euclidean space in which embedding the time series
* @return a matrix containing the embedded time series
* @details Each row of the resulting matrix consists of sequences x[t], x[t-1], …, x[t-dimension+1], where t is the original index of x.
* @note Implementation taken from (https://rdocumentation.org/packages/stats/versions/3.6.2)
*/
template<class LAG_policy>
KO_Traits::StoringMatrix
adf<LAG_policy>::embed(const KO_Traits::StoringVector &ts, int dimension)
const
{ 
  //evaluate one step differences for the time serie
  std::vector<double> y(this->one_step_diff(ts));

  int n = y.size();
  int m = n - dimension + 1;    //number of rows of embedded ts
  
  KO_Traits::StoringMatrix embedded(m,dimension);
  
  //embedding
  std::vector<std::size_t> start_ind;
  start_ind.resize(static_cast<std::size_t>(m));
  std::iota(start_ind.begin(),start_ind.end(),static_cast<std::size_t>(0));
  
  for(std::size_t i = 0; i < static_cast<std::size_t>(dimension); ++i)
  {
    std::size_t gain = static_cast<std::size_t>(dimension) - i - 1;
    
    std::transform(start_ind.begin(),start_ind.end(),
                   embedded.col(i).begin(),
                   [&y,&gain](std::size_t el){ return y[el+gain];});
    
  }
  
  start_ind.clear();
  
  return embedded;
}



/*!
* @brief Evaluating the ADF-test statistic for a time series
* @param ts time series
* @return statistic evaluation
*/
template<class LAG_policy>
double
adf<LAG_policy>::statistic_eval(const KO_Traits::StoringVector &ts)
const
{
  //embedding of the time serie in a low-dimensional Euclidean space
  KO_Traits::StoringMatrix z = this->embed(ts,m_k_used);
  
  //the staistics is computed evaluating coefficients standard errors of a fitted linear regression
  //the design matrix depends on the time lag: Lag_policy policy accounts for it
  adf_stat<LAG_policy> statistic_adf;

  //evaluate the statistic
  return statistic_adf(ts,z);
}
  


/*!
* @brief Evaluating thr ADF-test p-value for a time series
* @param ts time series
* @param tableipl table containing extreme values for the statistic for p-value computation
* @param i number of the discrete evaluation corresponding to 'ts'
* @return ADF-test p-value
*/
template<class LAG_policy>
double
adf<LAG_policy>::p_value_eval(const KO_Traits::StoringVector &ts, const std::vector<double> &tableipl, int i) 
const
{ 
  //evaluation of the test statistic
  double stat = this->statistic_eval(ts);
  
  //if statistic too extreme the pvalue is put to 0 or 1
  if(stat<tableipl.front()){return 0.0;}
  if(stat>tableipl.back()){return 1.0;}
  
  //pvalue evaluation
  interp_func p_val_int{tableipl};          //retain the p-value using interpolation: preparing the functor, with the discrete evaluations available
  double pval = p_val_int(stat,tablep);     //performing interpolation
  
  return pval;
}



/*!
* @brief computing the ADF-test pointwisely for each point of the domain for the fts
*/
template<class LAG_policy>
void
adf<LAG_policy>::test()
{
  
   //preparing for the p_value evaluation
   m_p_values.reserve(m_x.rows());
   //retaining the extreme values for statistic computation according to time series dimension
   std::vector<double> tableipl_ = tableipl(static_cast<double>(m_tot_time_instants-1));
   
   //computing the pvalues for each time serie
   for (int i = 0; i < m_x.rows(); ++i) 
   {
      m_p_values.emplace_back( this->p_value_eval(m_x.row(i),tableipl_,i+1)  );
   }
}