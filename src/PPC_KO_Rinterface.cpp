#include <RcppEigen.h>


#include "default_parameters.hpp"
#include "wrapper_params.hpp"

#include "reading_data.hpp"
#include "PPC_KO.hpp"
#include "KO_factory.hpp"




using namespace Rcpp;


//
// [[Rcpp::depends(RcppEigen)]]

//
// [[Rcpp::export]]
Rcpp::List PPC_KO(Rcpp::NumericMatrix X,
                  std::string id_CV,
                  std::string id_p_for_k,
                  double threshold_k,
                  double alpha = 0.75,
                  int n_disc = 1000,
                  int k = 0,
                  double alpha_min = 0.00000000001,
                  double alpha_max = 10,
                  std::string id_p_imposed = "No",
                  Rcpp::Nullable<std::string> id_rem_nan = R_NilValue
                  )
{ 
  using T = double;       //version for real-values time series
  
  //wrapping all the parameters
  const DEF_PARAMS_PPC::MA_type id_RN = WRAP_PARAMS_PPC::wrap_id_rem_nans(id_rem_nan);
  const int k_w = WRAP_PARAMS_PPC::wrap_k(k);
  
  
  const bool p_for_k = WRAP_PARAMS_PPC::wrap_id_p_for_k(id_p_for_k);
  const bool p_imposed = WRAP_PARAMS_PPC::wrap_id_p_imposed(id_p_imposed);
  
  
  //reading data, handling NANs
  KO_Traits::StoringMatrix x = reading_data<T>(X,id_RN);
  
  //object solver  !!using std::move() to move the data (better if there are big data)!!
  std::unique_ptr<PPC::PPC_KO_base> ko = KO_Factory::KO_solver(id_CV,std::move(x),threshold_k,p_for_k,p_imposed,alpha,n_disc,alpha_min,alpha_max,k_w);
  //solving
  ko->solve();
  
  //estimate of the predictions
  auto one_step_ahead_pred = ko->prediction();
  //alpha used
  double alpha_used = ko->alpha();
  //number of PPC retained
  int n_PPC = ko->k();
  
  auto valid_err = ko->ValidErr();
  

  //saving results in a list
  Rcpp::List l;
  l["predictions"] = one_step_ahead_pred;
  l["alpha"] = alpha_used;
  l["PPCs_retained"] = n_PPC;
  l["valid_errors"] = valid_err;
  
  return l;
}




//
// [[Rcpp::export]]
Rcpp::List read_data_na(Rcpp::NumericMatrix X)
{ 
  using T = double;       //version for double
  
  DEF_PARAMS_PPC::MA_type id_RN = DEF_PARAMS_PPC::MA_type::MR;
  KO_Traits::StoringMatrix x = reading_data<T>(X, id_RN);
  
  
  //saving results in a list
  Rcpp::List l;
  l["data_read"] = x;

  return l;
}