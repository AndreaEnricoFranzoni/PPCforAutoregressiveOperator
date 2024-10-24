#ifndef KO_READ_DATA_HPP
#define KO_READ_DATA_HPP

#include <RcppEigen.h>

#include <algorithm>
#include <iterator>
#include <numeric>
#include <cmath>
#include <iostream>

#include "traits_ko.hpp"
#include "removing_nan.hpp"
#include "parameters_wrapper.hpp"


//
// [[Rcpp::depends(RcppEigen)]]
template<typename T> 
const KO_Traits::StoringMatrix 
reader_data(Rcpp::NumericMatrix X,
              const REM_NAN MA_t)
{
  std::cout << "Local version 695" << std::endl;

  //taking the dimensions: n_row is the number of time series, n_col is the number of time istants
  int n_row = X.nrow();
  int n_col = X.ncol();
  
  //data are read into a Eigen::MatrixXd: Eigen::Map exploited (to be modified if the trait is not an Eigen object anymore)
  KO_Traits::StoringMatrix x = Eigen::Map<KO_Traits::StoringMatrix>(X.begin(),n_row,n_col);
  
  //check if there are nans
  auto check_nan = std::find_if(x.reshaped().cbegin(),x.reshaped().cend(),[](T el){return isnan(el);});
  
  //if there are nans: remove them
  if (check_nan!=x.reshaped().end())
  {
    if(MA_t == REM_NAN::MR)     //replacing nans with the mean
    {
      removing_nan<T,REM_NAN::MR> temp(std::move(x));
      temp.remove_nan();
      return temp.data();
    }
    if(MA_t == REM_NAN::ZR)     //replacing nans with 0s
    {
      removing_nan<T,REM_NAN::ZR> temp(std::move(x));
      temp.remove_nan();
      return temp.data();
    }
  }

  return x;
}

#endif /*KO_READ_DATA_HPP*/