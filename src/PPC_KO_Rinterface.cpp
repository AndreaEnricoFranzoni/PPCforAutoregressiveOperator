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
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> e9f642e91c5135471cf4581b5d36597c550b2d0d
=======
>>>>>>> 0cd8c8530975b720726228746d1d643ca4251c3c
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
