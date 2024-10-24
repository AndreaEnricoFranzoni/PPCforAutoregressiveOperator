#ifndef KO_WRAP_PARAMS_HPP
#define KO_WRAP_PARAMS_HPP

#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <string>
#include <stdexcept>


//#include "removing_nan.hpp"


//utilities to wrap the input parameters of the R function



//check that threshold_ppc is in the correct range
inline
void
check_threshold_ppc(const double &threshold_ppc)
{
  if(threshold_ppc<=0 || threshold_ppc>=1)
  {
    std::string error_message = "threshold_ppc has to be in (0,1)";
    throw std::invalid_argument(error_message);
  }
}



//check that alpha is in the correct range
inline
void
check_alpha(const double &alpha)
{
  if( alpha<= 0 )
  {
    std::string error_message = "alpha has to be a positive real number";
    throw std::invalid_argument(error_message);
  }
}



//check that k is in the correct range
inline
void
check_k(const int &k, const int &max_k)
{
  if( k < 0 )
  {
    std::string error_message = "k has to be a positive integer or 0";
    throw std::invalid_argument(error_message);
  }
  if( k > max_k )
  {
    std::string error_message = "k has to be lower than the maximum number of PPCs";
    throw std::invalid_argument(error_message);
  }
}



//to wrap the vector of alphas passed in input
inline
std::vector<double>
wrap_alpha_vec(Rcpp::Nullable<Rcpp::NumericVector> alpha_vec)
{
  //if no alpha is given, the default value for the alphas is given
  if(alpha_vec.isNull())
  {
    std::vector<double> alphas;
    alphas.resize(21);
    
    std::iota(alphas.begin(),alphas.end(),static_cast<double>(-10));
    std::transform(alphas.begin(),alphas.end(),alphas.begin(),[](double el){return(pow(static_cast<double>(10),el));});
    
    return alphas;
  }
  
  
  std::vector<double> alphas = Rcpp::as<std::vector<double>>(alpha_vec);
  
  //all the alphas have to be positive numbers
  if(*(std::min_element(alphas.cbegin(),alphas.cend())) <= 0)
  {
    std::string error_message = "Every alpha has to be a positive real number";
    throw std::invalid_argument(error_message);
  }
  
  //sorting into ascending order the alphas to be coherent during the algorithm
  std::sort(alphas.begin(), alphas.end());
  return alphas;
}



//to wrap the vector of k passed in input
inline
std::vector<int>
wrap_k_vec(Rcpp::Nullable<Rcpp::IntegerVector> k_vec, int k_max)
{
  //if no k is given, the default is looking for all the possible directions (number of discrete evaluations innthe domain of the functional object)
    if(k_vec.isNull())
    {
      std::vector<int> k_s;
      k_s.resize(k_max);
      
      std::iota(k_s.begin(),k_s.end(),static_cast<int>(1));
    
      return k_s;
    }
    
    
    std::vector<int> k_s = Rcpp::as<std::vector<int>>(k_vec);
    
    //all the k_s has to be in the range [1,k_max]
    if(*(std::min_element(k_s.cbegin(),k_s.cend())) < 1)
    {
      std::string error_message1 = "k has to be at least 1";
      throw std::invalid_argument(error_message1);
    }
    
    if(*(std::max_element(k_s.cbegin(),k_s.cend())) > k_max)
    {
      std::string error_message2 = "k cannot be greater than the number of discrete evaluation of the functional object in the domain (" + std::to_string(k_max) + ")";
      throw std::invalid_argument(error_message2);
    }
    
    //sorting into ascending order the alphas to be coherent during the algorithm
    std::sort(k_s.begin(), k_s.end());
    return k_s;
}



//removing NaNs
enum REM_NAN
{
  MR  = 0,      //replacing nans with mean (easily changes the mean of the distribution)
  ZR  = 1,      //replacing nans with 0s (easily changes the sd of the distribution)
};


//reads the input string an gives back the correct value of the enumerator for replacing nans
inline
REM_NAN
wrap_id_rem_nans(Rcpp::Nullable<std::string> id_rem_nan)
{
  if(id_rem_nan.isNull())
  { 
    return REM_NAN::MR;
  }
  if(Rcpp::as< std::string >(id_rem_nan) == "MR")
  {
    return REM_NAN::MR;
  }
  if(Rcpp::as< std::string >(id_rem_nan) == "ZR")
  {
    return REM_NAN::ZR;
  }
  else
  {
    std::string error_message = "Wrong input string for handling NANs";
    throw std::invalid_argument(error_message);
  }
  
};


#endif  /*KO_WRAP_PARAMS_HPP*/