#ifndef KO_READ_DATA_HPP
#define KO_READ_DATA_HPP

#include <RcppEigen.h>
#include "KO_traits.hpp"

//
// [[Rcpp::depends(RcppEigen)]]
inline const KO_Traits::StoringMatrix reading_data(Rcpp::NumericMatrix X)
{
  
  //taking the dimensions: n_row is the number of time series, n_col is the number of time istants
  int n_row = X.nrow();
  int n_col = X.ncol();
  
  //data are read into a Eigen::MatrixXd: Eigen::Map exploited (to be modified if the trait is not an Eigen object anymore)
  const KO_Traits::StoringMatrix x = Eigen::Map<KO_Traits::StoringMatrix>(X.begin(),n_row,n_col);
  
  return x;
}

#endif /*KO_READ_DATA_HPP*/