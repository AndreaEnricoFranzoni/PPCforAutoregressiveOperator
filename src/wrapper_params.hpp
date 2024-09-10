#ifndef KO_WRAP_PARAMS_HPP
#define KO_WRAP_PARAMS_HPP

#include <iostream>
#include <string>
#include <stdexcept>

#include "default_parameters.hpp"


namespace WRAP_PARAMS_PPC         //utilities to wrap the input parameters of the R function
{



inline
const int
wrap_k(int k)
{

  if(k < 0)
  {
    std::string error_message = "k has to be a positive integer";
    throw std::invalid_argument(error_message);
  }
  else 
  {  
    return static_cast<int>(k);
  }
  
};


//reads the input string that says if you impone p during CV k
inline
const bool
wrap_id_p_imposed(const std::string &id_p_imposed)
{
  if(id_p_imposed=="Yes")     //I am imposing that we are using k components to invert cov reg
  {
    return true;
  }
  if(id_p_imposed=="No")    //all components are used to invert rg, and then k is imposed on phi
  {
    return false;
  }
  else
  {
    std::string error_message = "Wrong input string for if you impone p";
    throw std::invalid_argument(error_message);
  }
}







//reads the input string that says how to choose the number of PPC
inline
const bool
wrap_id_p_for_k(const std::string &id_p_for_k)
{
  if(id_p_for_k=="Yes")     //only the p biggest components are used to invert rg, and then k=p
  {
    return true;
  }
  if(id_p_for_k=="No")    //all components are used to invert rg, and then k is evaluated on phi
  {
    return false;
  }
  else
  {
    std::string error_message = "Wrong input string for handling number of PPC";
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