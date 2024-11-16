#ifndef KO_READ_DATA_HPP
#define KO_READ_DATA_HPP

#include <RcppEigen.h>

#include <algorithm>
#include <iterator>
#include <numeric>
#include <cmath>
#include <utility>
#include <vector>
#include <iostream>

#include "traits_ko.hpp"
#include "removing_nan.hpp"
#include "parameters_wrapper.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif


//
// [[Rcpp::depends(RcppEigen)]]
template<typename T> 
std::pair<KO_Traits::StoringMatrix,std::vector<int>>
reader_data(Rcpp::NumericMatrix X,
            REM_NAN MA_t)
{
  std::cout << "Remote version 2.1" << std::endl;

#ifdef _OPENMP
  std::cout << "Remote par version 2.1.9" << std::endl;
#endif

  //taking the dimensions: n_row is the number of time series, n_col is the number of time istants
  int n_row = X.nrow();
  int n_col = X.ncol();
  
  
  //needed if we are simply mapping the 2dim grid
  if(MA_t == REM_NAN::NR)
  {
    KO_Traits::StoringMatrix x = Eigen::Map<KO_Traits::StoringMatrix>(X.begin(),n_row,n_col);
    std::vector<int> row_removed;
    return std::make_pair(x,row_removed);
  }
  
  
  //firstly: rows of all NaNs has to be removed (taking which ones)
  //rows to be retained and removed
  std::set<int> rows_to_be_retained = rows_entire_NaNs(X);      //rows to be retained
  
  std::vector<int> temp;
  temp.resize(X.nrow());
  std::iota(temp.begin(),temp.end(),static_cast<int>(0));     
  std::set<int> all_rows(temp.begin(),temp.end());
  temp.clear();
  
  std::set<int> rows_to_be_removed;                             //rows to be removed
  std::set_difference(all_rows.begin(),all_rows.end(),
                      rows_to_be_retained.begin(),rows_to_be_retained.end(),
                      std::inserter(rows_to_be_removed,rows_to_be_removed.begin()));
  all_rows.clear();
  
  //data matrix in which there are not rows of all NaNs and rows removed
  Rcpp::NumericMatrix X_clean = removing_NaNS_entire_rows(X,rows_to_be_retained);
  std::vector<int> rows_retained(rows_to_be_retained.begin(),rows_to_be_retained.end());   //AS IF INDECES START FROM 0
  
  
  //data are read into a Eigen::MatrixXd: Eigen::Map exploited (to be modified if the trait is not an Eigen object anymore)
  KO_Traits::StoringMatrix x = Eigen::Map<KO_Traits::StoringMatrix>(X_clean.begin(),rows_to_be_retained.size(),n_col);
  
  //check if there are nans
  auto check_nan = std::find_if(x.reshaped().cbegin(),x.reshaped().cend(),[](T el){return isnan(el);});
  
  //if there are nans: remove them
  if (check_nan!=x.reshaped().end())
  {
    if(MA_t == REM_NAN::MR)     //replacing nans with the mean
    {
      removing_nan<T,REM_NAN::MR> data_clean(std::move(x));
      data_clean.remove_nan();
      return std::make_pair(data_clean.data(),rows_retained);
    }
    if(MA_t == REM_NAN::ZR)     //replacing nans with 0s
    {
      removing_nan<T,REM_NAN::ZR> data_clean(std::move(x));
      data_clean.remove_nan();
      return std::make_pair(data_clean.data(),rows_retained);
    }
  }
  
  return std::make_pair(x,rows_retained);
}

#endif /*KO_READ_DATA_HPP*/
