#ifndef KO_WRAP_PARAMS_HPP
#define KO_WRAP_PARAMS_HPP

#include <iostream>
#include <string>
#include <stdexcept>

#include "default_parameters.hpp"


namespace WRAP_PARAMS_PPC         //utilities to wrap the input parameters of the R function
{

//reads the input string an gives back the correct value of the enumerator for replacing nans
inline
const DEF_PARAMS_PPC::MA_type
wrap_id_rem_nans(Rcpp::Nullable<std::string> id_rem_nan)
{
  if(id_rem_nan.isNull())
  { 
    return DEF_PARAMS_PPC::MA_type::EMA;
  }
  if(Rcpp::as< std::string >(id_rem_nan) == "EMA")
  {
    return DEF_PARAMS_PPC::MA_type::EMA;
  }
  if(Rcpp::as< std::string >(id_rem_nan) == "WMA")
  {
    return DEF_PARAMS_PPC::MA_type::WMA;
  }
  if(Rcpp::as< std::string >(id_rem_nan) == "SMA")
  {
    return DEF_PARAMS_PPC::MA_type::SMA;
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