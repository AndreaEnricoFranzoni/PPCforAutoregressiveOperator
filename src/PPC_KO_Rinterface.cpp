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
                  double threshold_k,
                  double alpha,
                  int n_disc,
                  Rcpp::Nullable<std::string> id_rem_nan = R_NilValue
                  )
{ 
  using T = double;       //version for real-values time series
  
  //wrapping all the parameters
  const DEF_PARAMS_PPC::MA_type id_RN = WRAP_PARAMS_PPC::wrap_id_rem_nans(id_rem_nan);
  
  //reading data, handling NANs
  KO_Traits::StoringMatrix x = reading_data<T>(X,id_RN);
  
  //object solver  !!using std::move() to move the data (better if there are big data)!!
  std::unique_ptr<PPC::PPC_KO_base> ko = KO_Factory::KO_solver(id_CV,std::move(x),alpha);
  
  //solving
  ko->solve();
  
  //estimate of the predictions
  auto pred = ko->prediction();
  //alpha used
  double alpha_used = ko->alpha();
  //number of PPC retained
  int k = ko->k();
  
  auto valid_err = ko->ValidErr();
  
  Rcout << "The regularization parameter used is " << alpha_used << std::endl;
  
  //saving results in a list
  Rcpp::List l;
  l["predictions"] = pred;
  l["alpha"] = alpha_used;
  l["PPCs_retained"] = k;
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