#include <RcppEigen.h>


#include "traits_ko.hpp"
#include "parameters_wrapper.hpp"
#include "utils.hpp"
#include "data_reader.hpp"
#include "Factory_ko.hpp"

#include "ADF_test.hpp"
#include "ADF_policies.hpp"


using namespace Rcpp;


//
// [[Rcpp::depends(RcppEigen)]]



//
// [[Rcpp::export]]
Rcpp::List PPC_KO(Rcpp::NumericMatrix           X,
                  std::string                   id_CV         = "NoCV",
                  double                        alpha         = 0.75,
                  int                           k             = 0, 
                  double                        threshold_ppc = 0.95,
                  Rcpp::Nullable<NumericVector> alpha_vec     = R_NilValue,
                  Rcpp::Nullable<IntegerVector> k_vec         = R_NilValue,
                  double                        toll          = 1e-4,
                  Rcpp::Nullable<NumericVector> disc_ev       = R_NilValue,
                  double                        left_extreme  = 0,
                  double                        right_extreme = 1,
                  int                           dom_dim_s     = 1,
                  int                           err_ret       = 0,
                  Rcpp::Nullable<std::string>   id_rem_nan    = R_NilValue
                  )
{ 
  using T = double;       //version for real-values time series
  
  //wrapping and checking parameters
  check_threshold_ppc(threshold_ppc);
  check_alpha(alpha);
  check_k(k,X.nrow());
  std::vector<double> alphas = wrap_alpha_vec(alpha_vec);
  std::vector<int> k_s       = wrap_k_vec(k_vec,X.nrow());
  const REM_NAN id_RN = wrap_id_rem_nans(id_rem_nan);
  std::vector<double> disc_ev_points = wrap_disc_ev(disc_ev,left_extreme,right_extreme,X.nrow());
  
  //reading data, handling NANs
  KO_Traits::StoringMatrix x = reader_data<T>(X,id_RN);
  
  //checking which 
  bool dom_dim = dom_dim_s==1 ? false : true;
  bool err_ret_b = err_ret==1 ? true : false;
  
  //returning element
  Rcpp::List l;

  if(!dom_dim)                                //DOMAIN 1D
  {
    
    if(err_ret_b)                                         //VALIDATION ERRORS STORED AND RETURNED
    {
      
      if(k>0)                                                                                   //K IMPOSED
      { 
        //1D domain, k imposed, returning errors
        //solver
        auto ko = KO_Factory< DOM_DIM::uni_dim, K_IMP::YES, VALID_ERR_RET::YES_err, CV_STRAT::AUGMENTING_WINDOW, CV_ERR_EVAL::MSE >::KO_solver(id_CV,std::move(x),alpha,k,threshold_ppc,alphas,k_s,toll);
        //solving
        ko->call_ko();
        //results
        auto one_step_ahead_pred  = std::get<0>(ko->results());   //estimate of the prediction
        double alpha_used         = std::get<1>(ko->results());   //alpha used
        int n_PPC                 = std::get<2>(ko->results());   //number of PPC retained
        auto scores_PPC           = std::get<3>(ko->results());   //scores along the k PPCs
        auto explanatory_power    = std::get<4>(ko->results());   //explanatory power
        auto directions           = std::get<5>(ko->results());   //directions
        auto weights              = std::get<6>(ko->results());   //weights
        auto valid_err            = std::get<7>(ko->results());   //valid errors
        Rcpp::List errors = valid_err_disp(valid_err);            //dispatching correctly valid errors
        
        //saving results in a list, that will be returned
        //Rcpp::List l;
        l["One-step ahead prediction"]  = one_step_ahead_pred;
        l["Alpha"]                     = alpha_used;
        l["Number of PPCs retained"]   = n_PPC;
        l["Scores along PPCs"]         = scores_PPC;
        l["Explanatory power PPCs"]    = explanatory_power;
        l["Directions of PPCs"]        = directions;
        l["Weights of PPCs"]           = weights;
        l["Validation errors"]         = errors["Errors"];
      }
      
      else                                                                                      //K NOT IMPOSED
      { 
        //1D domain, k not imposed, returning errors
        //solver
        auto ko = KO_Factory< DOM_DIM::uni_dim, K_IMP::NO, VALID_ERR_RET::YES_err, CV_STRAT::AUGMENTING_WINDOW, CV_ERR_EVAL::MSE >::KO_solver(id_CV,std::move(x),alpha,k,threshold_ppc,alphas,k_s,toll);
        //solving
        ko->call_ko();
        //results
        auto one_step_ahead_pred  = std::get<0>(ko->results());   //estimate of the prediction
        double alpha_used         = std::get<1>(ko->results());   //alpha used
        int n_PPC                 = std::get<2>(ko->results());   //number of PPC retained
        auto scores_PPC           = std::get<3>(ko->results());   //scores along the k PPCs
        auto explanatory_power    = std::get<4>(ko->results());   //explanatory power
        auto directions           = std::get<5>(ko->results());   //directions
        auto weights              = std::get<6>(ko->results());   //weights
        auto valid_err            = std::get<7>(ko->results());   //valid errors
        Rcpp::List errors = valid_err_disp(valid_err);            //dispatching correctly valid errors
        
        //saving results in a list, that will be returned
        //Rcpp::List l;
        l["One-step ahead prediction"]  = one_step_ahead_pred;
        l["Alpha"]                     = alpha_used;
        l["Number of PPCs retained"]   = n_PPC;
        l["Scores along PPCs"]         = scores_PPC;
        l["Explanatory power PPCs"]    = explanatory_power;
        l["Directions of PPCs"]        = directions;
        l["Weights of PPCs"]           = weights;
        l["Validation errors"]         = errors["Errors"];
        
      }
    }
    
    
    else                                                //VALIDATION ERRORS NOT STORED
    {
      
      if(k>0)                                                                                   //K IMPOSED
      {
        //1D domain, k imposed, not returning errors
        //solver
        auto ko = KO_Factory< DOM_DIM::uni_dim, K_IMP::YES, VALID_ERR_RET::NO_err, CV_STRAT::AUGMENTING_WINDOW, CV_ERR_EVAL::MSE >::KO_solver(id_CV,std::move(x),alpha,k,threshold_ppc,alphas,k_s,toll);
        //solving
        ko->call_ko();
        //results
        auto one_step_ahead_pred  = std::get<0>(ko->results());   //estimate of the prediction
        double alpha_used         = std::get<1>(ko->results());   //alpha used
        int n_PPC                 = std::get<2>(ko->results());   //number of PPC retained
        auto scores_PPC           = std::get<3>(ko->results());   //scores along the k PPCs
        auto explanatory_power    = std::get<4>(ko->results());   //explanatory power
        auto directions           = std::get<5>(ko->results());   //directions
        auto weights              = std::get<6>(ko->results());   //weights
        //auto valid_err            = std::get<7>(ko->results());   //valid errors
        //Rcpp::List errors = valid_err_disp(valid_err);            //dispatching correctly valid errors
        
        //saving results in a list, that will be returned
        //Rcpp::List l;
        l["One-step ahead prediction"]  = one_step_ahead_pred;
        l["Alpha"]                     = alpha_used;
        l["Number of PPCs retained"]   = n_PPC;
        l["Scores along PPCs"]         = scores_PPC;
        l["Explanatory power PPCs"]    = explanatory_power;
        l["Directions of PPCs"]        = directions;
        l["Weights of PPCs"]           = weights;
        //l["Validation errors"]         = errors["Errors"];
      }
      
      else                                                                                      //K NOT IMPOSED
      {
        //1D domain, k not imposed, not returning errors
        //solver
        auto ko = KO_Factory< DOM_DIM::uni_dim, K_IMP::NO, VALID_ERR_RET::NO_err, CV_STRAT::AUGMENTING_WINDOW, CV_ERR_EVAL::MSE >::KO_solver(id_CV,std::move(x),alpha,k,threshold_ppc,alphas,k_s,toll);
        //solving
        ko->call_ko();
        //results
        auto one_step_ahead_pred  = std::get<0>(ko->results());   //estimate of the prediction
        double alpha_used         = std::get<1>(ko->results());   //alpha used
        int n_PPC                 = std::get<2>(ko->results());   //number of PPC retained
        auto scores_PPC           = std::get<3>(ko->results());   //scores along the k PPCs
        auto explanatory_power    = std::get<4>(ko->results());   //explanatory power
        auto directions           = std::get<5>(ko->results());   //directions
        auto weights              = std::get<6>(ko->results());   //weights
        //auto valid_err            = std::get<7>(ko->results());   //valid errors
        //Rcpp::List errors = valid_err_disp(valid_err);            //dispatching correctly valid errors
        
        //saving results in a list, that will be returned
        //Rcpp::List l;
        l["One-step ahead prediction"]  = one_step_ahead_pred;
        l["Alpha"]                     = alpha_used;
        l["Number of PPCs retained"]   = n_PPC;
        l["Scores along PPCs"]         = scores_PPC;
        l["Explanatory power PPCs"]    = explanatory_power;
        l["Directions of PPCs"]        = directions;
        l["Weights of PPCs"]           = weights;
        //l["Validation errors"]         = errors["Errors"];
      }        
    }
  }
  else                                        //DOMAIN 2D
  {
    
    if(err_ret_b)                                         //VALIDATION ERRORS STORED AND RETURNED
    {
      
      if(k>0)                                                                                   //K IMPOSED
      {
        //2D domain, k imposed, returning errors
        //solver
        auto ko = KO_Factory< DOM_DIM::bi_dim, K_IMP::YES, VALID_ERR_RET::YES_err, CV_STRAT::AUGMENTING_WINDOW, CV_ERR_EVAL::MSE >::KO_solver(id_CV,std::move(x),alpha,k,threshold_ppc,alphas,k_s,toll);
        //solving
        ko->call_ko();
        //results
        auto one_step_ahead_pred  = std::get<0>(ko->results());   //estimate of the prediction
        double alpha_used         = std::get<1>(ko->results());   //alpha used
        int n_PPC                 = std::get<2>(ko->results());   //number of PPC retained
        auto scores_PPC           = std::get<3>(ko->results());   //scores along the k PPCs
        auto explanatory_power    = std::get<4>(ko->results());   //explanatory power
        auto directions           = std::get<5>(ko->results());   //directions
        auto weights              = std::get<6>(ko->results());   //weights
        auto valid_err            = std::get<7>(ko->results());   //valid errors
        Rcpp::List errors = valid_err_disp(valid_err);            //dispatching correctly valid errors
        
        //saving results in a list, that will be returned
        //Rcpp::List l;
        l["One-step ahead prediction"]  = one_step_ahead_pred;
        l["Alpha"]                     = alpha_used;
        l["Number of PPCs retained"]   = n_PPC;
        l["Scores along PPCs"]         = scores_PPC;
        l["Explanatory power PPCs"]    = explanatory_power;
        l["Directions of PPCs"]        = directions;
        l["Weights of PPCs"]           = weights;
        l["Validation errors"]         = errors["Errors"];
      }
      
      else                                                                                      //K NOT IMPOSED
      {
        //2D domain, k not imposed, returning errors
        //solver
        auto ko = KO_Factory< DOM_DIM::bi_dim, K_IMP::NO, VALID_ERR_RET::YES_err, CV_STRAT::AUGMENTING_WINDOW, CV_ERR_EVAL::MSE >::KO_solver(id_CV,std::move(x),alpha,k,threshold_ppc,alphas,k_s,toll);
        //solving
        ko->call_ko();
        //results
        auto one_step_ahead_pred  = std::get<0>(ko->results());   //estimate of the prediction
        double alpha_used         = std::get<1>(ko->results());   //alpha used
        int n_PPC                 = std::get<2>(ko->results());   //number of PPC retained
        auto scores_PPC           = std::get<3>(ko->results());   //scores along the k PPCs
        auto explanatory_power    = std::get<4>(ko->results());   //explanatory power
        auto directions           = std::get<5>(ko->results());   //directions
        auto weights              = std::get<6>(ko->results());   //weights
        auto valid_err            = std::get<7>(ko->results());   //valid errors
        Rcpp::List errors = valid_err_disp(valid_err);            //dispatching correctly valid errors
        
        //saving results in a list, that will be returned
        //Rcpp::List l;
        l["One-step ahead prediction"]  = one_step_ahead_pred;
        l["Alpha"]                     = alpha_used;
        l["Number of PPCs retained"]   = n_PPC;
        l["Scores along PPCs"]         = scores_PPC;
        l["Explanatory power PPCs"]    = explanatory_power;
        l["Directions of PPCs"]        = directions;
        l["Weights of PPCs"]           = weights;
        l["Validation errors"]         = errors["Errors"];
      }
    }
    
    
    else                                                //VALIDATION ERRORS NOT STORED
    {
      
      if(k>0)                                                                                   //K IMPOSED
      { 
        //2D domain, k imposed, not returning errors
        //solver
        auto ko = KO_Factory< DOM_DIM::bi_dim, K_IMP::YES, VALID_ERR_RET::NO_err, CV_STRAT::AUGMENTING_WINDOW, CV_ERR_EVAL::MSE >::KO_solver(id_CV,std::move(x),alpha,k,threshold_ppc,alphas,k_s,toll);
        //solving
        ko->call_ko();
        //results
        auto one_step_ahead_pred  = std::get<0>(ko->results());   //estimate of the prediction
        double alpha_used         = std::get<1>(ko->results());   //alpha used
        int n_PPC                 = std::get<2>(ko->results());   //number of PPC retained
        auto scores_PPC           = std::get<3>(ko->results());   //scores along the k PPCs
        auto explanatory_power    = std::get<4>(ko->results());   //explanatory power
        auto directions           = std::get<5>(ko->results());   //directions
        auto weights              = std::get<6>(ko->results());   //weights
        //auto valid_err            = std::get<7>(ko->results());   //valid errors
        //Rcpp::List errors = valid_err_disp(valid_err);            //dispatching correctly valid errors
        
        //saving results in a list, that will be returned
        l["One-step ahead prediction"]  = one_step_ahead_pred;
        l["Alpha"]                     = alpha_used;
        l["Number of PPCs retained"]   = n_PPC;
        l["Scores along PPCs"]         = scores_PPC;
        l["Explanatory power PPCs"]    = explanatory_power;
        l["Directions of PPCs"]        = directions;
        l["Weights of PPCs"]           = weights;
        //l["Validation errors"]         = errors["Errors"];
      }
      
      else                                                                                      //K NOT IMPOSED
      {
        //2D domain, k not imposed, not returning errors
        //solver
        auto ko = KO_Factory< DOM_DIM::bi_dim, K_IMP::NO, VALID_ERR_RET::NO_err, CV_STRAT::AUGMENTING_WINDOW, CV_ERR_EVAL::MSE >::KO_solver(id_CV,std::move(x),alpha,k,threshold_ppc,alphas,k_s,toll);
        //solving
        ko->call_ko();
        //results
        auto one_step_ahead_pred  = std::get<0>(ko->results());   //estimate of the prediction
        double alpha_used         = std::get<1>(ko->results());   //alpha used
        int n_PPC                 = std::get<2>(ko->results());   //number of PPC retained
        auto scores_PPC           = std::get<3>(ko->results());   //scores along the k PPCs
        auto explanatory_power    = std::get<4>(ko->results());   //explanatory power
        auto directions           = std::get<5>(ko->results());   //directions
        auto weights              = std::get<6>(ko->results());   //weights
        //auto valid_err            = std::get<7>(ko->results());   //valid errors
        //Rcpp::List errors = valid_err_disp(valid_err);            //dispatching correctly valid errors
        
        //saving results in a list, that will be returned
        l["One-step ahead prediction"]  = one_step_ahead_pred;
        l["Alpha"]                     = alpha_used;
        l["Number of PPCs retained"]   = n_PPC;
        l["Scores along PPCs"]         = scores_PPC;
        l["Explanatory power PPCs"]    = explanatory_power;
        l["Directions of PPCs"]        = directions;
        l["Weights of PPCs"]           = weights;
        //l["Validation errors"]         = errors["Errors"];
      }        
    }    
  }
  
  
  l["Function discrete evaluations points"] = disc_ev_points;
  l["Left extreme domain"]                  = left_extreme; 
  l["Right extreme domain"]                 = right_extreme;
  l["f_n"]                                  = X(_, X.ncol() - 1);
  l["CV"]                                   = id_CV;
  l["Alphas"]                               = alphas;
  l["K_s"]                                  = k_s;
  
  return l;
}



//
// [[Rcpp::export]]
Rcpp::List KO_check_hps(Rcpp::NumericMatrix X)
{
  using T = double; 
  
  KO_Traits::StoringMatrix x = reader_data<T>(X, REM_NAN::MR);
  
  int number_time_instants = x.cols();
  
  std::size_t k = static_cast<std::size_t>(std::trunc(std::cbrt(static_cast<double>(number_time_instants)-1)));
  
  Rcpp::List l;
  
  if(k > 1)
  {
    //ADF with lag order k
    adf<CaseLagOrderADF> adf_t(std::move(x),k);
    adf_t.test();
    auto pv = adf_t.p_values();
    
    l["Pvalues ADF"] = pv;
    return l;
  }
  
  //ADF without lag order
  adf<CaseLagOrderADF> adf_t(std::move(x),0);
  adf_t.test();
  auto pv = adf_t.p_values();
  
  l["Pvalues ADF"] = pv;
  return l;
}







//
// [[Rcpp::export]]
Rcpp::List read_data_na(Rcpp::NumericMatrix X,std::string id_CV)
{ 
  using T = double;       //version for double
  
  REM_NAN id_RN = REM_NAN::MR;
  KO_Traits::StoringMatrix x = reader_data<T>(X, id_RN);
  
  /*
   * TEST t(id_CV);
   t.testing();
   
   auto alpha = std::get<0>(t.test_res());
   auto err   = std::get<1>(t.test_res());
   
   //saving results in a list
   
   l["alpha"] = alpha;
   
   Rcpp::List errors = valid_err_disp(err);
   l["Err"] = errors["Errors"];
   */

  Rcpp::List l;
  l["data_read"] = x;
  return l;
}