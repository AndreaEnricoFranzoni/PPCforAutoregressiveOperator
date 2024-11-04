#ifndef KO_UTILS_HPP
#define KO_UTILS_HPP

#include "traits_ko.hpp"
#include <limits>

//function to dispacth the error depending if is from CV single-param or double param
Rcpp::List 
valid_err_disp
(const valid_err_variant& var) 
{
  
  return std::visit([](auto&& arg) -> Rcpp::List 
  {
    using T = std::decay_t<decltype(arg)>;
    
    // CV single param: errors in a vector -> only one vector back
    if constexpr (std::is_same_v<T, std::vector<double>>) 
    {
      return Rcpp::List::create(Rcpp::Named("Errors") = Rcpp::wrap(arg));  // NumericVector
    }
    // CV double param: errors in a vector of vector -> list of vector(one for every alpha) back
    else if constexpr (std::is_same_v<T, std::vector<std::vector<double>>>) 
    {
      Rcpp::List res(arg.size());
      for (size_t i = 0; i < arg.size(); ++i) 
      {
        res[i] = Rcpp::wrap(arg[i]);  // Ogni std::vector<double> diventa un NumericVector
      }
      
      return Rcpp::List::create(Rcpp::Named("Errors") = res);  // Lista di NumericVector
    } 
    
    else 
    {
      Rcpp::stop("Wrong type in variant!");
    }
    }, var);
}


//function to map a matrix into a vector: for 2d case
KO_Traits::StoringVector
from_matrix_to_col(const KO_Traits::StoringMatrix &mat)
{
  return KO_Traits::StoringVector(mat.reshaped());
}

//function to map a vector into a matrix: for 2d case
KO_Traits::StoringMatrix
from_col_to_matrix(const KO_Traits::StoringVector &col, int rows, int cols)
{
  return Eigen::Map<const KO_Traits::StoringMatrix>(col.data(),rows,cols);
}


//function to add nans in the result
//pred is the vector with the prediction, row ret the vector containing which points actually have evals, complete size is the size of all the points
KO_Traits::StoringVector
add_nans_vec(const KO_Traits::StoringVector &pred, const std::vector<int> &row_ret, int complete_size)
{
  if(row_ret.size()==0){return pred;}
  KO_Traits::StoringVector pred_comp(complete_size);
  
  pred_comp.setConstant(std::numeric_limits<double>::quiet_NaN());
  
  int counter=0;
  std::for_each(row_ret.cbegin(),row_ret.cend(),[&pred_comp,&pred,&counter](int el){pred_comp[el]=pred[counter]; counter++;});
  
  return pred_comp;
}

#endif  //KO_UTILS_HPP