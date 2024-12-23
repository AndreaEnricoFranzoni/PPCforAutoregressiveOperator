#include <RcppEigen.h>

#include <string>
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
                  Rcpp::Nullable<int>           min_size_ts   = R_NilValue,
                  Rcpp::Nullable<int>           max_size_ts   = R_NilValue,
                  int                           err_ret       = 0,
                  Rcpp::Nullable<int>           num_threads   = R_NilValue,
                  Rcpp::Nullable<std::string>   id_rem_nan    = R_NilValue
                  )
{ 
  //1D DOMAIN
  using T = double;       //version for real-values time series
  
  //wrapping and checking parameters
  check_threshold_ppc(threshold_ppc);
  check_alpha(alpha);
  check_k(k,X.nrow());
  std::vector<double> alphas         = wrap_alpha_vec(alpha_vec);
  std::vector<int> k_s               = wrap_k_vec(k_vec,X.nrow());
  const REM_NAN id_RN                = wrap_id_rem_nans(id_rem_nan);
  std::vector<double> disc_ev_points = wrap_disc_ev(disc_ev,left_extreme,right_extreme,X.nrow());
  auto sizes_CV_sets                 = wrap_sizes_set_CV(min_size_ts,max_size_ts,X.ncol());
  int min_dim_train_set              = sizes_CV_sets.first;
  int max_dim_train_set              = sizes_CV_sets.second;
  int number_threads                 = wrap_num_thread(num_threads);

  //reading data, handling NANs
  auto data_read = reader_data<T>(X,id_RN);
  KO_Traits::StoringMatrix x = data_read.first;
  
  //checking which 
  bool err_ret_b = err_ret==1 ? true : false;
  
  //returning element
  Rcpp::List l;
  
  
  
  Rcout << "--------------------------------------------------------------------------------------------" << std::endl;
  Rcout << "Running Kargin-Onatski algorithm, " << wrap_string_CV_to_be_printed(id_CV) << std::endl;
  Rcout << "Functional data defined over: [" << left_extreme << "," << right_extreme << "], with " << disc_ev_points.size() << " discrete evaluations" << std::endl;

  
  
  if(err_ret_b)                                         //VALIDATION ERRORS STORED AND RETURNED
  {
    
    if(k>0)                                                                                   //K IMPOSED
    { 
      //1D domain, k imposed, returning errors
      //solver
      auto ko = KO_Factory< DOM_DIM::uni_dim, K_IMP::YES, VALID_ERR_RET::YES_err, CV_STRAT::AUGMENTING_WINDOW, CV_ERR_EVAL::MSE >::KO_solver(id_CV,std::move(x),alpha,k,threshold_ppc,alphas,k_s,toll,min_dim_train_set,max_dim_train_set,number_threads);
      //solving
      ko->call_ko();
      //results
      auto one_step_ahead_pred  = add_nans_vec(std::get<0>(ko->results()),data_read.second,X.nrow());  //estimate of the prediction (NaN for the points in which you do not have measurements)
      double alpha_used         = std::get<1>(ko->results());   //alpha used
      int n_PPC                 = std::get<2>(ko->results());   //number of PPC retained
      auto scores_PPC           = std::get<3>(ko->results());   //scores along the k PPCs
      auto explanatory_power    = std::get<4>(ko->results());   //explanatory power
      auto directions           = std::get<5>(ko->results());
      auto weights              = std::get<6>(ko->results());
      auto valid_err            = std::get<7>(ko->results());   //valid errors
      Rcpp::List errors = valid_err_disp(valid_err);            //dispatching correctly valid errors
      Rcpp::List directions_wrapped;
      Rcpp::List weights_wrapped;
      for(std::size_t i = 0; i < n_PPC; ++i)
      {
        std::string name_d = "Direction PPC " + std::to_string(i+1);
        std::string name_w = "Weight PPC " + std::to_string(i+1);
        directions_wrapped[name_d] = add_nans_vec(directions.col(i),data_read.second,X.nrow());
        weights_wrapped[name_w] = add_nans_vec(weights.col(i),data_read.second,X.nrow());
      }
      
      //saving results in a list, that will be returned
      l["One-step ahead prediction"] = one_step_ahead_pred;
      l["Alpha"]                     = alpha_used;
      l["Number of PPCs retained"]   = n_PPC;
      l["Scores along PPCs"]         = scores_PPC;
      l["Explanatory power PPCs"]    = explanatory_power;
      l["Directions of PPCs"]        = directions_wrapped;
      l["Weights of PPCs"]           = weights_wrapped;
      l["Validation errors"]         = errors["Errors"];
    }
    
    else                                                                                      //K NOT IMPOSED
    { 
      //1D domain, k not imposed, returning errors
      //solver
      auto ko = KO_Factory< DOM_DIM::uni_dim, K_IMP::NO, VALID_ERR_RET::YES_err, CV_STRAT::AUGMENTING_WINDOW, CV_ERR_EVAL::MSE >::KO_solver(id_CV,std::move(x),alpha,k,threshold_ppc,alphas,k_s,toll,min_dim_train_set,max_dim_train_set,number_threads);
      //solving
      ko->call_ko();
      //results
      auto one_step_ahead_pred  = add_nans_vec(std::get<0>(ko->results()),data_read.second,X.nrow());  //estimate of the prediction (NaN for the points in which you do not have measurements)
      double alpha_used         = std::get<1>(ko->results());   //alpha used
      int n_PPC                 = std::get<2>(ko->results());   //number of PPC retained
      auto scores_PPC           = std::get<3>(ko->results());   //scores along the k PPCs
      auto explanatory_power    = std::get<4>(ko->results());   //explanatory power
      auto directions           = std::get<5>(ko->results());
      auto weights              = std::get<6>(ko->results());
      auto valid_err            = std::get<7>(ko->results());   //valid errors
      Rcpp::List errors = valid_err_disp(valid_err);            //dispatching correctly valid errors
      Rcpp::List directions_wrapped;
      Rcpp::List weights_wrapped;
      for(std::size_t i = 0; i < n_PPC; ++i)
      {
        std::string name_d = "Direction PPC " + std::to_string(i+1);
        std::string name_w = "Weight PPC " + std::to_string(i+1);
        directions_wrapped[name_d] = add_nans_vec(directions.col(i),data_read.second,X.nrow());
        weights_wrapped[name_w] = add_nans_vec(weights.col(i),data_read.second,X.nrow());
      }
      
      //saving results in a list, that will be returned
      l["One-step ahead prediction"] = one_step_ahead_pred;
      l["Alpha"]                     = alpha_used;
      l["Number of PPCs retained"]   = n_PPC;
      l["Scores along PPCs"]         = scores_PPC;
      l["Explanatory power PPCs"]    = explanatory_power;
      l["Directions of PPCs"]        = directions_wrapped;
      l["Weights of PPCs"]           = weights_wrapped;
      l["Validation errors"]         = errors["Errors"];
      
    }
  }
  
  
  else                                                //VALIDATION ERRORS NOT STORED
  {
    
    if(k>0)                                                                                   //K IMPOSED
    {
      //1D domain, k imposed, not returning errors
      //solver
      auto ko = KO_Factory< DOM_DIM::uni_dim, K_IMP::YES, VALID_ERR_RET::NO_err, CV_STRAT::AUGMENTING_WINDOW, CV_ERR_EVAL::MSE >::KO_solver(id_CV,std::move(x),alpha,k,threshold_ppc,alphas,k_s,toll,min_dim_train_set,max_dim_train_set,number_threads);
      //solving
      ko->call_ko();
      //results
      auto one_step_ahead_pred  = add_nans_vec(std::get<0>(ko->results()),data_read.second,X.nrow());  //estimate of the prediction (NaN for the points in which you do not have measurements)
      double alpha_used         = std::get<1>(ko->results());   //alpha used
      int n_PPC                 = std::get<2>(ko->results());   //number of PPC retained
      auto scores_PPC           = std::get<3>(ko->results());   //scores along the k PPCs
      auto explanatory_power    = std::get<4>(ko->results());   //explanatory power
      auto directions           = std::get<5>(ko->results());
      auto weights              = std::get<6>(ko->results());
      Rcpp::List directions_wrapped;
      Rcpp::List weights_wrapped;
      for(std::size_t i = 0; i < n_PPC; ++i)
      {
        std::string name_d = "Direction PPC " + std::to_string(i+1);
        std::string name_w = "Weight PPC " + std::to_string(i+1);
        directions_wrapped[name_d] = add_nans_vec(directions.col(i),data_read.second,X.nrow());
        weights_wrapped[name_w] = add_nans_vec(weights.col(i),data_read.second,X.nrow());
      }

      //saving results in a list, that will be returned
      l["One-step ahead prediction"] = one_step_ahead_pred;
      l["Alpha"]                     = alpha_used;
      l["Number of PPCs retained"]   = n_PPC;
      l["Scores along PPCs"]         = scores_PPC;
      l["Explanatory power PPCs"]    = explanatory_power;
      l["Directions of PPCs"]        = directions_wrapped;
      l["Weights of PPCs"]           = weights_wrapped;
    }
    
    else                                                                                      //K NOT IMPOSED
    {
      //1D domain, k not imposed, not returning errors
      //solver
      auto ko = KO_Factory< DOM_DIM::uni_dim, K_IMP::NO, VALID_ERR_RET::NO_err, CV_STRAT::AUGMENTING_WINDOW, CV_ERR_EVAL::MSE >::KO_solver(id_CV,std::move(x),alpha,k,threshold_ppc,alphas,k_s,toll,min_dim_train_set,max_dim_train_set,number_threads);
      //solving
      ko->call_ko();
      //results
      auto one_step_ahead_pred  = add_nans_vec(std::get<0>(ko->results()),data_read.second,X.nrow());  //estimate of the prediction (NaN for the points in which you do not have measurements)
      double alpha_used         = std::get<1>(ko->results());   //alpha used
      int n_PPC                 = std::get<2>(ko->results());   //number of PPC retained
      auto scores_PPC           = std::get<3>(ko->results());   //scores along the k PPCs
      auto explanatory_power    = std::get<4>(ko->results());   //explanatory power
      auto directions           = std::get<5>(ko->results());
      auto weights              = std::get<6>(ko->results());
      Rcpp::List directions_wrapped;
      Rcpp::List weights_wrapped;
      for(std::size_t i = 0; i < n_PPC; ++i)
      {
        std::string name_d = "Direction PPC " + std::to_string(i+1);
        std::string name_w = "Weight PPC " + std::to_string(i+1);
        directions_wrapped[name_d] = add_nans_vec(directions.col(i),data_read.second,X.nrow());
        weights_wrapped[name_w] = add_nans_vec(weights.col(i),data_read.second,X.nrow());
      }

      //saving results in a list, that will be returned
      l["One-step ahead prediction"] = one_step_ahead_pred;
      l["Alpha"]                     = alpha_used;
      l["Number of PPCs retained"]   = n_PPC;
      l["Scores along PPCs"]         = scores_PPC;
      l["Explanatory power PPCs"]    = explanatory_power;
      l["Directions of PPCs"]        = directions_wrapped;
      l["Weights of PPCs"]           = weights_wrapped;
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
Rcpp::List PPC_KO_2d(Rcpp::NumericMatrix           X,
                     std::string                   id_CV            = "NoCV",
                     double                        alpha            = 0.75,
                     int                           k                = 0, 
                     double                        threshold_ppc    = 0.95,
                     Rcpp::Nullable<NumericVector> alpha_vec        = R_NilValue,
                     Rcpp::Nullable<IntegerVector> k_vec            = R_NilValue,
                     double                        toll             = 1e-4,
                     Rcpp::Nullable<NumericVector> disc_ev_x1       = R_NilValue,
                     int                           num_disc_ev_x1   = 10,
                     Rcpp::Nullable<NumericVector> disc_ev_x2       = R_NilValue,
                     int                           num_disc_ev_x2   = 10,
                     double                        left_extreme_x1  = 0,
                     double                        right_extreme_x1 = 1,
                     double                        left_extreme_x2  = 0,
                     double                        right_extreme_x2 = 1,
                     Rcpp::Nullable<int>           min_size_ts      = R_NilValue,
                     Rcpp::Nullable<int>           max_size_ts      = R_NilValue,
                     int                           err_ret          = 0,
                     Rcpp::Nullable<int>           num_threads      = R_NilValue,
                     Rcpp::Nullable<std::string>   id_rem_nan       = R_NilValue
)
{ 
  //2D DOMAIN
  using T = double;       //version for real-values time series
  
  //wrapping and checking parameters
  check_threshold_ppc(threshold_ppc);
  check_alpha(alpha);
  check_k(k,X.nrow());
  std::vector<double> alphas = wrap_alpha_vec(alpha_vec);
  std::vector<int> k_s       = wrap_k_vec(k_vec,X.nrow());
  const REM_NAN id_RN = wrap_id_rem_nans(id_rem_nan);
  std::vector<double> disc_ev_points_x1 = wrap_disc_ev(disc_ev_x1,left_extreme_x1,right_extreme_x1,num_disc_ev_x1);
  std::vector<double> disc_ev_points_x2 = wrap_disc_ev(disc_ev_x2,left_extreme_x2,right_extreme_x2,num_disc_ev_x2);
  auto sizes_CV_sets                    = wrap_sizes_set_CV(min_size_ts,max_size_ts,X.ncol());
  int min_dim_train_set                 = sizes_CV_sets.first;
  int max_dim_train_set                 = sizes_CV_sets.second;
  int number_threads                    = wrap_num_thread(num_threads);

  //reading data, handling NANs
  auto data_read = reader_data<T>(X,id_RN);
  KO_Traits::StoringMatrix x = data_read.first;
  
  //checking which 
  bool err_ret_b = err_ret==1 ? true : false;
  
  //returning element
  Rcpp::List l;
  
  Rcout << "--------------------------------------------------------------------------------------------" << std::endl;
  Rcout << "Running Kargin-Onatski algorithm, " << wrap_string_CV_to_be_printed(id_CV) << std::endl;
  Rcout << "Functional data defined over: [" << left_extreme_x1 << "," << right_extreme_x1 << "] x [" << left_extreme_x2 << "," << right_extreme_x2 <<"], with " << disc_ev_points_x1.size() << " x " << disc_ev_points_x2.size() << " discrete evaluations" << std::endl;
  
  
  
  if(err_ret_b)                                         //VALIDATION ERRORS STORED AND RETURNED
  {
    
    if(k>0)                                                                                   //K IMPOSED
    {
      //2D domain, k imposed, returning errors
      //solver
      auto ko = KO_Factory< DOM_DIM::bi_dim, K_IMP::YES, VALID_ERR_RET::YES_err, CV_STRAT::AUGMENTING_WINDOW, CV_ERR_EVAL::MSE >::KO_solver(id_CV,std::move(x),alpha,k,threshold_ppc,alphas,k_s,toll,min_dim_train_set,max_dim_train_set,number_threads);
      //solving
      ko->call_ko();
      //results
      auto one_step_ahead_pred  = from_col_to_matrix(add_nans_vec(std::get<0>(ko->results()),data_read.second,X.nrow()),disc_ev_points_x1.size(),disc_ev_points_x2.size());  //estimate of the prediction (NaN for the points in which you do not have measurements)
      double alpha_used         = std::get<1>(ko->results());   //alpha used
      int n_PPC                 = std::get<2>(ko->results());   //number of PPC retained
      auto scores_PPC           = std::get<3>(ko->results());   //scores along the k PPCs
      auto explanatory_power    = std::get<4>(ko->results());   //explanatory power
      auto directions           = std::get<5>(ko->results());   //directions
      auto weights              = std::get<6>(ko->results());   //weights
      auto valid_err            = std::get<7>(ko->results());   //valid errors
      Rcpp::List errors = valid_err_disp(valid_err);            //dispatching correctly valid errors
      Rcpp::List directions_wrapped;
      Rcpp::List weights_wrapped;
      for(std::size_t i = 0; i < n_PPC; ++i)
      {
        std::string name_d = "Direction PPC " + std::to_string(i+1);
        std::string name_w = "Weight PPC " + std::to_string(i+1);
        directions_wrapped[name_d] = from_col_to_matrix(add_nans_vec(directions.col(i),data_read.second,X.nrow()),disc_ev_points_x1.size(),disc_ev_points_x2.size());
        weights_wrapped[name_w] = from_col_to_matrix(add_nans_vec(weights.col(i),data_read.second,X.nrow()),disc_ev_points_x1.size(),disc_ev_points_x2.size());
      }
      
      l["One-step ahead prediction"] = one_step_ahead_pred;
      l["Alpha"]                     = alpha_used;
      l["Number of PPCs retained"]   = n_PPC;
      l["Scores along PPCs"]         = scores_PPC;
      l["Explanatory power PPCs"]    = explanatory_power;
      l["Directions of PPCs"]        = directions_wrapped;
      l["Weights of PPCs"]           = weights_wrapped;
      l["Validation errors"]         = errors["Errors"];
    }
    
    else                                                                                      //K NOT IMPOSED
    {
      //2D domain, k not imposed, returning errors
      //solver
      auto ko = KO_Factory< DOM_DIM::bi_dim, K_IMP::NO, VALID_ERR_RET::YES_err, CV_STRAT::AUGMENTING_WINDOW, CV_ERR_EVAL::MSE >::KO_solver(id_CV,std::move(x),alpha,k,threshold_ppc,alphas,k_s,toll,min_dim_train_set,max_dim_train_set,number_threads);
      //solving
      ko->call_ko();
      //results
      auto one_step_ahead_pred  = from_col_to_matrix(add_nans_vec(std::get<0>(ko->results()),data_read.second,X.nrow()),disc_ev_points_x1.size(),disc_ev_points_x2.size());  //estimate of the prediction (NaN for the points in which you do not have measurements)
      double alpha_used         = std::get<1>(ko->results());   //alpha used
      int n_PPC                 = std::get<2>(ko->results());   //number of PPC retained
      auto scores_PPC           = std::get<3>(ko->results());   //scores along the k PPCs
      auto explanatory_power    = std::get<4>(ko->results());   //explanatory power
      auto directions           = std::get<5>(ko->results());   //directions
      auto weights              = std::get<6>(ko->results());   //weights
      auto valid_err            = std::get<7>(ko->results());   //valid errors
      Rcpp::List errors = valid_err_disp(valid_err);            //dispatching correctly valid errors
      Rcpp::List directions_wrapped;
      Rcpp::List weights_wrapped;
      for(std::size_t i = 0; i < n_PPC; ++i)
      {
        std::string name_d = "Direction PPC " + std::to_string(i+1);
        std::string name_w = "Weight PPC " + std::to_string(i+1);
        directions_wrapped[name_d] = from_col_to_matrix(add_nans_vec(directions.col(i),data_read.second,X.nrow()),disc_ev_points_x1.size(),disc_ev_points_x2.size());
        weights_wrapped[name_w] = from_col_to_matrix(add_nans_vec(weights.col(i),data_read.second,X.nrow()),disc_ev_points_x1.size(),disc_ev_points_x2.size());
      }
      
      //saving results in a list, that will be returned
      //Rcpp::List l;
      l["One-step ahead prediction"] = one_step_ahead_pred;
      l["Alpha"]                     = alpha_used;
      l["Number of PPCs retained"]   = n_PPC;
      l["Scores along PPCs"]         = scores_PPC;
      l["Explanatory power PPCs"]    = explanatory_power;
      l["Directions of PPCs"]        = directions_wrapped;
      l["Weights of PPCs"]           = weights_wrapped;
      l["Validation errors"]         = errors["Errors"];
    }
  }
  
  
  else                                                //VALIDATION ERRORS NOT STORED
  {
    
    if(k>0)                                                                                   //K IMPOSED
    { 
      //2D domain, k imposed, not returning errors
      //solver
      auto ko = KO_Factory< DOM_DIM::bi_dim, K_IMP::YES, VALID_ERR_RET::NO_err, CV_STRAT::AUGMENTING_WINDOW, CV_ERR_EVAL::MSE >::KO_solver(id_CV,std::move(x),alpha,k,threshold_ppc,alphas,k_s,toll,min_dim_train_set,max_dim_train_set,number_threads);
      //solving
      ko->call_ko();
      //results
      auto one_step_ahead_pred  = from_col_to_matrix(add_nans_vec(std::get<0>(ko->results()),data_read.second,X.nrow()),disc_ev_points_x1.size(),disc_ev_points_x2.size());  //estimate of the prediction (NaN for the points in which you do not have measurements)
      double alpha_used         = std::get<1>(ko->results());   //alpha used
      int n_PPC                 = std::get<2>(ko->results());   //number of PPC retained
      auto scores_PPC           = std::get<3>(ko->results());   //scores along the k PPCs
      auto explanatory_power    = std::get<4>(ko->results());   //explanatory power
      auto directions           = std::get<5>(ko->results());   //directions
      auto weights              = std::get<6>(ko->results());   //weights
      Rcpp::List directions_wrapped;
      Rcpp::List weights_wrapped;
      for(std::size_t i = 0; i < n_PPC; ++i)
      {
        std::string name_d = "Direction PPC " + std::to_string(i+1);
        std::string name_w = "Weight PPC " + std::to_string(i+1);
        directions_wrapped[name_d] = from_col_to_matrix(add_nans_vec(directions.col(i),data_read.second,X.nrow()),disc_ev_points_x1.size(),disc_ev_points_x2.size());
        weights_wrapped[name_w] = from_col_to_matrix(add_nans_vec(weights.col(i),data_read.second,X.nrow()),disc_ev_points_x1.size(),disc_ev_points_x2.size());
      }

      //saving results in a list, that will be returned
      l["One-step ahead prediction"] = one_step_ahead_pred;
      l["Alpha"]                     = alpha_used;
      l["Number of PPCs retained"]   = n_PPC;
      l["Scores along PPCs"]         = scores_PPC;
      l["Explanatory power PPCs"]    = explanatory_power;
      l["Directions of PPCs"]        = directions_wrapped;
      l["Weights of PPCs"]           = weights_wrapped;
    }
    
    else                                                                                      //K NOT IMPOSED
    {
      //2D domain, k not imposed, not returning errors
      //solver
      auto ko = KO_Factory< DOM_DIM::bi_dim, K_IMP::NO, VALID_ERR_RET::NO_err, CV_STRAT::AUGMENTING_WINDOW, CV_ERR_EVAL::MSE >::KO_solver(id_CV,std::move(x),alpha,k,threshold_ppc,alphas,k_s,toll,min_dim_train_set,max_dim_train_set,number_threads);
      //solving
      ko->call_ko();
      //results
      auto one_step_ahead_pred  = from_col_to_matrix(add_nans_vec(std::get<0>(ko->results()),data_read.second,X.nrow()),disc_ev_points_x1.size(),disc_ev_points_x2.size());  //estimate of the prediction (NaN for the points in which you do not have measurements)
      double alpha_used         = std::get<1>(ko->results());   //alpha used
      int n_PPC                 = std::get<2>(ko->results());   //number of PPC retained
      auto scores_PPC           = std::get<3>(ko->results());   //scores along the k PPCs
      auto explanatory_power    = std::get<4>(ko->results());   //explanatory power
      auto directions           = std::get<5>(ko->results());   //directions
      auto weights              = std::get<6>(ko->results());   //weights
      Rcpp::List directions_wrapped;
      Rcpp::List weights_wrapped;
      for(std::size_t i = 0; i < n_PPC; ++i)
      {
        std::string name_d = "Direction PPC " + std::to_string(i+1);
        std::string name_w = "Weight PPC " + std::to_string(i+1);
        directions_wrapped[name_d] = from_col_to_matrix(add_nans_vec(directions.col(i),data_read.second,X.nrow()),disc_ev_points_x1.size(),disc_ev_points_x2.size());
        weights_wrapped[name_w] = from_col_to_matrix(add_nans_vec(weights.col(i),data_read.second,X.nrow()),disc_ev_points_x1.size(),disc_ev_points_x2.size());
        
      }
      
      //saving results in a list, that will be returned
      l["One-step ahead prediction"] = one_step_ahead_pred;
      l["Alpha"]                     = alpha_used;
      l["Number of PPCs retained"]   = n_PPC;
      l["Scores along PPCs"]         = scores_PPC;
      l["Explanatory power PPCs"]    = explanatory_power;
      l["Directions of PPCs"]        = directions_wrapped;
      l["Weights of PPCs"]           = weights_wrapped;
    }        
  }  
  
  NumericMatrix f_n(disc_ev_points_x1.size(),disc_ev_points_x2.size());
  std::copy(X(_, X.ncol() - 1).begin(), X(_, X.ncol() - 1).end(), f_n.begin());
  
  l["Function discrete evaluations points dim1"] = disc_ev_points_x1;
  l["Left extreme domain dim1"]                  = left_extreme_x1; 
  l["Right extreme domain dim1"]                 = right_extreme_x1;
  l["Function discrete evaluations points dim2"] = disc_ev_points_x2;
  l["Left extreme domain dim2"]                  = left_extreme_x2; 
  l["Right extreme domain dim2"]                 = right_extreme_x2;
  l["f_n"]                                       = f_n;
  l["CV"]                                        = id_CV;
  l["Alphas"]                                    = alphas;
  l["K_s"]                                       = k_s;
  
  return l;
}



//
// [[Rcpp::export]]
Rcpp::List KO_check_hps(Rcpp::NumericMatrix X)
{
  using T = double; 
  
  Rcout << "--------------------------------------" << std::endl;
  Rcout << "Pointwise ADF test p-values evaluation" << std::endl;
  
  auto data_read = reader_data<T>(X, REM_NAN::MR);
  KO_Traits::StoringMatrix x = data_read.first;
  
  int number_time_instants = x.cols();
  
  std::size_t k = static_cast<std::size_t>(std::trunc(std::cbrt(static_cast<double>(number_time_instants)-1)));
  
  Rcpp::List l;
  
  if(k > 1)
  {
    //ADF with lag order k
    adf<CaseLagOrderADF> adf_t(std::move(x),k);
    adf_t.test();
    auto pv = adf_t.p_values();
    auto pv_final = add_nans_vec(Eigen::Map<const KO_Traits::StoringVector>(pv.data(),pv.size()),data_read.second,X.nrow());
    
    l["P-values ADF"] = pv_final;
    return l;
  }
  
  //ADF without lag order
  adf<CaseLagOrderADF> adf_t(std::move(x),0);
  adf_t.test();
  auto pv = adf_t.p_values();
  auto pv_final = add_nans_vec(Eigen::Map<const KO_Traits::StoringVector>(pv.data(),pv.size()),data_read.second,X.nrow());
  
  
  l["P-values ADF"] = pv_final;
  return l;
}



//
// [[Rcpp::export]]
Rcpp::List KO_check_hps_2d(Rcpp::NumericMatrix X, int dim_x1, int dim_x2)
{ 
  //X is the matrix containing the 2d grid already dispatched
  using T = double; 
  
  Rcout << "--------------------------------------" << std::endl;
  Rcout << "Pointwise ADF test p-values evaluation" << std::endl;
  
  auto data_read = reader_data<T>(X, REM_NAN::MR);
  KO_Traits::StoringMatrix x = data_read.first;
  
  int number_time_instants = x.cols();
  
  std::size_t k = static_cast<std::size_t>(std::trunc(std::cbrt(static_cast<double>(number_time_instants)-1)));
  
  Rcpp::List l;
  
  if(k > 1)
  {
    //ADF with lag order k
    adf<CaseLagOrderADF> adf_t(std::move(x),k);
    adf_t.test();
    auto pv = adf_t.p_values();
    auto pv_final = from_col_to_matrix(add_nans_vec(Eigen::Map<const KO_Traits::StoringVector>(pv.data(),pv.size()),data_read.second,X.nrow()),dim_x1,dim_x2);
    
    l["P-values ADF"] = pv_final;
    return l;
  }
  
  //ADF without lag order
  adf<CaseLagOrderADF> adf_t(std::move(x),0);
  adf_t.test();
  auto pv = adf_t.p_values();
  auto pv_final = from_col_to_matrix(add_nans_vec(Eigen::Map<const KO_Traits::StoringVector>(pv.data(),pv.size()),data_read.second,X.nrow()),dim_x1,dim_x2);
  
  l["P-values ADF"] = pv_final;
  return l;
}










//
// [[Rcpp::export]]
Rcpp::NumericMatrix data_2d_wrapper_from_list(Rcpp::List Xt)
{
  using T = double;
  
  //this works only for 1-step time series
  int number_time_instants = Xt.size();
  if(number_time_instants==0)
  {
    std::string error_message1 = "Empty list";
    throw std::invalid_argument(error_message1);
  }
  
  int number_point_evaluations = as<NumericMatrix>(Xt[0]).nrow()*as<NumericMatrix>(Xt[0]).ncol();
  if(number_point_evaluations==0)
  {
    std::string error_message2 = "List of empty matrices";
    throw std::invalid_argument(error_message2);
  }
  
  Rcpp::NumericMatrix x(number_point_evaluations,number_time_instants);
  
  for(std::size_t i = 0; i < static_cast<std::size_t>(number_time_instants); ++i)
  {
    KO_Traits::StoringMatrix col = reader_data<T>(as<NumericMatrix>(Xt[i]),REM_NAN::NR).first;
    KO_Traits::StoringVector first_col = from_matrix_to_col(col);
    Rcpp::NumericVector col_wrapped = Rcpp::wrap(first_col);
    x(_,i) = col_wrapped;
  }
   
  return x;
}









//
// [[Rcpp::export]]
Rcpp::NumericMatrix data_2d_wrapper_from_array(Rcpp::NumericVector Xt)
{
  using T = double;
  
  IntegerVector dimensions = Xt.attr("dim");
  int dim1_size = dimensions[0];
  int dim2_size = dimensions[1];
  int number_point_evaluations = dim1_size*dim2_size;
  int number_time_instants = dimensions[2];   //this works only for 1-step time series
  
  Rcout << "Dim1: " << dim1_size << ", dim2: " << dim2_size << ", times: " << number_time_instants << std::endl;
  if(number_time_instants==0)
  {
    std::string error_message1 = "Empty array";
    throw std::invalid_argument(error_message1);
  }
  
  //object that will be returned
  Rcpp::NumericMatrix x(number_point_evaluations,number_time_instants);
  
  for(std::size_t i = 0; i < static_cast<std::size_t>(number_time_instants); ++i)
  { 
    //to save the result for each instant
    Rcpp::NumericMatrix instant_data(dim1_size, dim2_size);
    
    for (int i1 = 0; i1 < dim1_size; ++i1) {
      for (int i2 = 0; i2 < dim2_size; ++i2) {
        instant_data(i1, i2) = Xt[i1 + i2*dim1_size + i*dim1_size*dim2_size];  
      }
    }
    
    //wrapping into usual data format
    KO_Traits::StoringMatrix col = reader_data<T>(instant_data,REM_NAN::NR).first;
    KO_Traits::StoringVector first_col = from_matrix_to_col(col);
    Rcpp::NumericVector col_wrapped = Rcpp::wrap(first_col);
    x(_,i) = col_wrapped;
  }
  
  return x;
}



//
// [[Rcpp::export]]
Rcpp::List test(int sz)
{
  std::vector<std::array<double,2>> res;
  res.reserve(sz);
  
  for(int i = 0; i < sz; ++i){
    res.emplace_back(std::array<double,2>{0.2+i,0.3+i});
  }
  
  Rcpp::List l;
  l["Res"] = res;
  
  return l;
}