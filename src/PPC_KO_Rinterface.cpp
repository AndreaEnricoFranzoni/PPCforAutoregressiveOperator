#include <RcppEigen.h>
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
                  int n_disc
                  )
{ 
  Rcout << "Version 565" << std::endl;
  KO_Traits::StoringMatrix x = reading_data(X);

  
  //!!using std::move() to move the data (better if there are big data)!!
  std::unique_ptr<PPC::PPC_KO_base> ko = KO_Factory::KO_solver(id_CV,std::move(x),alpha);

  ko->solve();
  //estimate of the predictions
  auto pred = ko->prediction();
  double alpha_used = ko->alpha();
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