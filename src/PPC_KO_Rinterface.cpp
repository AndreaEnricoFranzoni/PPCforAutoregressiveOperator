#include <RcppEigen.h>


#include "default_parameters.hpp"
#include "wrapper_params.hpp"

#include "reading_data.hpp"
#include "PPC_KO.hpp"
#include "KO_factory.hpp"

#include "ztest_reg.hpp"


using namespace Rcpp;


//
// [[Rcpp::depends(RcppEigen)]]


//
// [[Rcpp::export]]
Rcpp::List PPC_KO(Rcpp::NumericMatrix X,
                  std::string id_CV = "NoCV",
                  double threshold_ppc = 0.95,
                  double alpha = 0.75,
                  int k = 0,
                  Rcpp::Nullable<std::string> id_rem_nan = R_NilValue
                  )
{ 
  using T = double;       //version for real-values time series
  
  //wrapping parameters
  const DEF_PARAMS_PPC::MA_type id_RN = WRAP_PARAMS_PPC::wrap_id_rem_nans(id_rem_nan);
  
  
  //reading data, handling NANs
  KO_Traits::StoringMatrix x = reading_data<T>(X,id_RN);
  
  
  //checking parameters
  WRAP_PARAMS_PPC::check_threshold_ppc(threshold_ppc);
  WRAP_PARAMS_PPC::check_alpha(alpha);
  WRAP_PARAMS_PPC::check_k(k,x.rows());
  
  
  //object solver  !!using std::move() to move the data (better if there are big data)!!
  std::unique_ptr<PPC::PPC_KO_base> ko = KO_Factory::KO_solver(id_CV,std::move(x),threshold_ppc,alpha,k);
  
  //solving using the requested KO algo
  ko->solve();

  //estimate of the predictions
  auto one_step_ahead_pred  = ko->prediction();
  //alpha used
  double alpha_used         = ko->alpha();
  //number of PPC retained
  int n_PPC                 = ko->k();
  //scores along the k PPCs
  auto scores_PPC           = ko->scores();
  //estimate of rho
  auto rho_estimate         = ko->rho();
  
  auto valid_err = ko->ValidErr();
  

  //saving results in a list
  Rcpp::List l;
  l["predictions"]    = one_step_ahead_pred;
  l["alpha"]          = alpha_used;
  l["PPCs_retained"]  = n_PPC;
  l["scores"]         = scores_PPC;
  l["rho_hat"]        = rho_estimate;
  l["valid_errors"]   = valid_err;
  
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


//
// [[Rcpp::export]]
Rcpp::List test_lr(Rcpp::NumericMatrix X, Rcpp::NumericMatrix Y)
{
  using T = double; 
  
  KO_Traits::StoringMatrix x = reading_data<T>(X, DEF_PARAMS_PPC::MA_type::MR);
  KO_Traits::StoringMatrix y = reading_data<T>(Y, DEF_PARAMS_PPC::MA_type::MR);
  
  reg_test rg(std::move(x),std::move(y));
  
  rg.solve();
  
  auto coeff = rg.coeff();
  
  Rcpp::List l;
  l["coeff"] = coeff;
  return l;
}