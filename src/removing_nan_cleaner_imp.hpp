#include "removing_nan.hpp"
#include <RcppEigen.h>

//remove rows of nan

// Funzione che verifica se una riga è interamente NaN
inline
bool 
is_row_all_nan(const Rcpp::NumericMatrix::ConstRow& row, int tot_cols) 
{
  return std::count_if(row.cbegin(),row.cend(),[](auto el){return isnan(el);}) == tot_cols  ?  true  :  false;
}


//function to identify rows of only NaNs
//
//  [[Rcpp::depends(RcppEigen)]]
std::set<int>
rows_entire_NaNs(const Rcpp::NumericMatrix &x)
{
  if (x.nrow() == 0 || x.ncol() == 0) 
  {
    std::string error_message1 = "Empty data matrix";
    throw std::invalid_argument(error_message1);
  }
    
  int tot_cols = x.ncol();
  std::set<int> rows_to_be_kept;
  for (int i = 0; i < x.nrow(); ++i) { if (!is_row_all_nan(x.row(i),tot_cols)) {  rows_to_be_kept.insert(i);}}
    
  if (rows_to_be_kept.empty())
  {
    std::string error_message2 = "Only-NaNs data matrix";
    throw std::invalid_argument(error_message2);
  }
    
  return rows_to_be_kept;
}


// function to remove the all-NaNs rows
//
//  [[Rcpp::depends(RcppEigen)]]
Rcpp::NumericMatrix 
removing_NaNS_entire_rows(const Rcpp::NumericMatrix &x,const std::set<int> &rows_to_be_kept)
{
  Rcpp::NumericMatrix x_clean(rows_to_be_kept.size(),x.ncol());
  int counter_row_clean = 0;
  std::for_each(rows_to_be_kept.cbegin(),rows_to_be_kept.cend(),[&x,&x_clean,&counter_row_clean](auto el){x_clean.row(counter_row_clean)=x.row(el); counter_row_clean++;});
  
  return x_clean;
}