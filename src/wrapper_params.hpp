#ifndef KO_WRAP_PARAMS_HPP
#define KO_WRAP_PARAMS_HPP

#include <iostream>
#include <string>
#include <stdexcept>

#include "default_parameters.hpp"


namespace WRAP_PARAMS_PPC         //utilities to wrap the input parameters of the R function
{


//check that threshold_ppc is in the correct range
inline
void
check_threshold_ppc(double threshold_ppc)
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
check_alpha(double alpha)
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
check_k(int k, int max_k)
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



//reads the input string an gives back the correct value of the enumerator for replacing nans
inline
const DEF_PARAMS_PPC::MA_type
wrap_id_rem_nans(Rcpp::Nullable<std::string> id_rem_nan)
{
  if(id_rem_nan.isNull())
  { 
    return DEF_PARAMS_PPC::MA_type::MR;
  }
  if(Rcpp::as< std::string >(id_rem_nan) == "MR")
  {
    return DEF_PARAMS_PPC::MA_type::MR;
  }
  if(Rcpp::as< std::string >(id_rem_nan) == "ZR")
  {
    return DEF_PARAMS_PPC::MA_type::ZR;
  }
  else
  {
    std::string error_message = "Wrong input string for handling NANs";
    throw std::invalid_argument(error_message);
  }
  
};

} //end namespace WRAP_PARAMS_PPC

#endif  /*KO_WRAP_PARAMS_HPP*/