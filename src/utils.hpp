#ifndef KO_UTILS_HPP
#define KO_UTILS_HPP

#include "traits_ko.hpp"

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

#endif  //KO_UTILS_HPP