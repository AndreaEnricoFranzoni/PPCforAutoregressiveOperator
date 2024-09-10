#ifndef KO_READ_DATA_HPP
#define KO_READ_DATA_HPP

#include <RcppEigen.h>

#include <algorithm>
#include <iterator>
#include <numeric>
#include <cmath>
#include <iostream>

#include "KO_Traits.hpp"
#include "default_parameters.hpp"
#include "removing_nan.hpp"
#include "wrapper_params.hpp"


//
// [[Rcpp::depends(RcppEigen)]]
template<typename T> 
const KO_Traits::StoringMatrix 
reading_data(Rcpp::NumericMatrix X,
              const DEF_PARAMS_PPC::MA_type MA_t)
{
  //std::cout << "Local version 162" << std::endl;

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

    if(MA_t == DEF_PARAMS_PPC::MA_type::EMA)    //replacing nans with exponential moving average
    {
      REM_NAN_PPC::removing_nan<T,DEF_PARAMS_PPC::MA_type::EMA> temp(std::move(x));
      temp.remove_nan();
      return temp.data();
    }
    if(MA_t == DEF_PARAMS_PPC::MA_type::WMA)    //replacing nans with weighted moving average
    {
      REM_NAN_PPC::removing_nan<T,DEF_PARAMS_PPC::MA_type::WMA> temp(std::move(x));
      temp.remove_nan();
      return temp.data();
    }
    if(MA_t == DEF_PARAMS_PPC::MA_type::SMA)    //replacing nans with simple moving average
    {
      REM_NAN_PPC::removing_nan<T,DEF_PARAMS_PPC::MA_type::SMA> temp(std::move(x));
      temp.remove_nan();
      return temp.data();
    }
    if(MA_t == DEF_PARAMS_PPC::MA_type::MR)     //replacing nans with the mean
    {
      REM_NAN_PPC::removing_nan<T,DEF_PARAMS_PPC::MA_type::MR> temp(std::move(x));
      temp.remove_nan();
      return temp.data();
    }
    if(MA_t == DEF_PARAMS_PPC::MA_type::ZR)     //replacing nans with 0s
    {
      REM_NAN_PPC::removing_nan<T,DEF_PARAMS_PPC::MA_type::ZR> temp(std::move(x));
      temp.remove_nan();
      return temp.data();
    }
  }


  return x;
}

#endif /*KO_READ_DATA_HPP*/