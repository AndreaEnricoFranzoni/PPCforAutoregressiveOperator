// Copyright (c) 2024 Andrea Enrico Franzoni (andreaenrico.franzoni@gmail.com)
//
// This file is part of PPCKO
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of PPCKO and associated documentation files (the PPCKO software), to deal
// PPCKO without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of PPCKO, and to permit persons to whom PPCKO is
// furnished to do so, subject to the following conditions:
//
// PPCKO IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH PPCKO OR THE USE OR OTHER DEALINGS IN
// PPCKO.

#include <RcppEigen.h>

#include <string>
#include "traits_ko.hpp"
#include "parameters_wrapper.hpp"
#include "utils.hpp"
#include "data_reader.hpp"
#include "Factory_ko.hpp"

#include "ADF_test.hpp"
#include "ADF_policies.hpp"


/*!
* @file PPC_KO_Rinterface.cpp
* @brief Contains the R-interfaced functions of the package 'PPCKO', which implement Principal Predictive Components (PPC)
*        Kargin-Onatski (KO) algorithm to perform one-step ahead prediction of Functional Time Series (FTS).
* @author Andrea Enrico Franzoni
*/


using namespace Rcpp;

//
// [[Rcpp::depends(RcppEigen)]]


/*!
* @brief Function to perform one-step ahead prediction of Functional Time Series of curves using PPCKO.
* @param X Rcpp::NumericMatrix (matrix of double) containing the curve time series: each row (m) is the evaluation of the curve in a point of its domain, each column (n) a time instant
* @param id_CV string denoting which version of the algorithm the user wants to use: 'NoCV', 'CV_alpha', 'CV_k' or 'CV' 
* @param alpha regularization parameter (positive real number)
* @param k number of PPCs: if 0, is selected through explanatory power criterion; if between 1 and m: k is imposed
* @param threshold_ppc minimum requested proportion of explanatory power: used only if k=0. Duble between 0 and 1
* @param alpha_vec vector of double containing the input space for cross-validation over the regularization parameter 
* @param k_vec vector of integers containing the input space for cross-validation over the number of retained PPCs
* @param toll tolerance over which cross-validation process stops looking for bigger number of PPCs if the validation errors between two consecutive iterations is smaller than it
* @param disc_ev vector of double containing the points of the domain for which the curve's evaluations are available
* @param left_extreme double indicating the left extreme of the domain of the curve
* @param right_extreme double indicating the right extreme of the domain of the curve
* @param min_size_ts time instants that define the smallest training set (between 2 and max_size_ts)
* @param max_size_ts time instants that define the biggest training set (between min_size_ts and n-1)
* @param err_ret true if validation errors are stored and returned, false if not
* @param ex_solver true if solving PPCKO inverting the regularized covariance matrix, false if relaying on GEP to avoid id
* @param num_threads number of threads to be used in OMP parallel directives
* @param id_rem_nan string that defines how to handle NaNs for some instant: 'MR': replacing them with the mean of the fts in that point, 'ZR' with 0s
* @return an R list containing:
* - one step ahead prediction of the fts
* - used regularization parameter
* - number of retained PPCs
* - scores along the PPCs (scalar product between the directions and fts at last available instant)
* - cumulative explanatory power of the PPCs (coherent only if ex_solver true)
* - directions of the PPCs
* - weights of the PPCs
* - standard deviation of the scores of each direction (evaluate on the scalar products between direction and the curve between instants 2 and n)
* - standard deviation of the scores of each weight (evaluate on the scalar products between weight and the curve between instants 1 and n-1)
* - mean function of the fts
* - validation errors (only if err_ret true)
* - points of the domain for which the evaluation of the curve are available
* - left extreme of the curve's domain
* - right extreme of the curve's domain
* - curve at the last available instant
* - which PPCKO version has been performed
* - input space for regularization parameter
* - input space for the number of PPCs
*/
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
                  bool                          err_ret       = false,
                  bool                          ex_solver     = true,
                  Rcpp::Nullable<int>           num_threads   = R_NilValue,
                  Rcpp::Nullable<std::string>   id_rem_nan    = R_NilValue
                  )
{ 
  using T = double;                   //real-values functional time series
  
  //wrapping and checking parameters
  check_threshold_ppc(threshold_ppc);
  check_alpha(alpha);
  check_k(k,X.nrow());
  check_solver(ex_solver,id_CV,k);
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
    
  //returning element
  Rcpp::List l;
  
  Rcout << "--------------------------------------------------------------------------------------------" << std::endl;
  Rcout << "Running Kargin-Onatski algorithm, " << wrap_string_CV_to_be_printed(id_CV) << std::endl;
  Rcout << "Functional data defined over: [" << left_extreme << "," << right_extreme << "], with " << disc_ev_points.size() << " discrete evaluations" << std::endl;

  if(ex_solver)               //EXACT ALGORITHM FOR PPCs
  {

    if(err_ret)                                         //VALIDATION ERRORS STORED AND RETURNED
    {

      if(k>0)                                                                                   //K IMPOSED
      { 
        //exact solver, k imposed, returning errors
        //solver
        auto ko = KO_Factory< SOLVER::ex_solver, K_IMP::YES, VALID_ERR_RET::YES_err, CV_STRAT::AUGMENTING_WINDOW, CV_ERR_EVAL::MSE >::KO_solver(id_CV,std::move(x),alpha,k,threshold_ppc,alphas,k_s,toll,min_dim_train_set,max_dim_train_set,number_threads);
        //solving
        ko->call_ko();
        //results
        auto one_step_ahead_pred  = add_nans_vec(std::get<0>(ko->results()),data_read.second,X.nrow());  //estimate of the prediction (NaN for the points in which you do not have measurements)
        double alpha_used         = std::get<1>(ko->results());                                          //alpha used
        int n_PPC                 = std::get<2>(ko->results());                                          //number of retained PPCs
        auto scores_PPC           = std::get<3>(ko->results());                                          //scores along the k PPCs
        auto explanatory_power    = std::get<4>(ko->results());                                          //explanatory power
        auto directions           = std::get<5>(ko->results());                                          //PPCs directions
        auto weights              = std::get<6>(ko->results());                                          //PPCs weights
        auto sd_scores_dir_wei    = std::get<7>(ko->results());                                          //sd of scores of directions and weights
        auto mean_func            = add_nans_vec(std::get<8>(ko->results()),data_read.second,X.nrow());  //mean function
        auto valid_err            = std::get<9>(ko->results());                                          //validation errors
        //dispatching correctly validation errors
        Rcpp::List errors = valid_err_disp(valid_err);
        //wrapping directions and weights to be returned properly            
        Rcpp::List directions_wrapped;
        Rcpp::List weights_wrapped;
        std::vector<double> scores_dir_sd;
        scores_dir_sd.reserve(n_PPC);
        std::vector<double> scores_wei_sd;
        scores_wei_sd.reserve(n_PPC);
        for(std::size_t i = 0; i < n_PPC; ++i)
        {
          std::string name_d = "Direction PPC " + std::to_string(i+1);
          std::string name_w = "Weight PPC " + std::to_string(i+1);
          directions_wrapped[name_d] = add_nans_vec(directions.col(i),data_read.second,X.nrow());
          weights_wrapped[name_w] = add_nans_vec(weights.col(i),data_read.second,X.nrow());
          scores_dir_sd.emplace_back(sd_scores_dir_wei[i][0]);
          scores_wei_sd.emplace_back(sd_scores_dir_wei[i][1]);
        }
      
        //saving results in a list, that will be returned
        l["One-step ahead prediction"] = one_step_ahead_pred;
        l["Alpha"]                     = alpha_used;
        l["Number of PPCs retained"]   = n_PPC;
        l["Scores along PPCs"]         = scores_PPC;
        l["Explanatory power PPCs"]    = explanatory_power;
        l["Directions of PPCs"]        = directions_wrapped;
        l["Weights of PPCs"]           = weights_wrapped;
        l["Sd scores directions"]      = scores_dir_sd;
        l["Sd scores weights"]         = scores_wei_sd;
        l["Mean function"]             = mean_func;
        l["Validation errors"]         = errors["Errors"];
      }
    
      else                                                                                      //K NOT IMPOSED
      { 
        //1D domain, k not imposed, returning errors
        //solver
        auto ko = KO_Factory< SOLVER::ex_solver, K_IMP::NO, VALID_ERR_RET::YES_err, CV_STRAT::AUGMENTING_WINDOW, CV_ERR_EVAL::MSE >::KO_solver(id_CV,std::move(x),alpha,k,threshold_ppc,alphas,k_s,toll,min_dim_train_set,max_dim_train_set,number_threads);
        //solving
        ko->call_ko();
        //results
        auto one_step_ahead_pred  = add_nans_vec(std::get<0>(ko->results()),data_read.second,X.nrow());  //estimate of the prediction (NaN for the points in which you do not have measurements)
        double alpha_used         = std::get<1>(ko->results());                                          //alpha used
        int n_PPC                 = std::get<2>(ko->results());                                          //number of retained PPCs
        auto scores_PPC           = std::get<3>(ko->results());                                          //scores along the k PPCs
        auto explanatory_power    = std::get<4>(ko->results());                                          //explanatory power
        auto directions           = std::get<5>(ko->results());                                          //PPCs directions
        auto weights              = std::get<6>(ko->results());                                          //PPcs weights
        auto sd_scores_dir_wei    = std::get<7>(ko->results());                                          //sd of scores of directions and weights
        auto mean_func            = add_nans_vec(std::get<8>(ko->results()),data_read.second,X.nrow());  //mean function
        auto valid_err            = std::get<9>(ko->results());                                          //validation errors
        //dispatching correctly validation errors
        Rcpp::List errors = valid_err_disp(valid_err);  
        //wrapping directions and weights to be returned properly              
        Rcpp::List directions_wrapped;
        Rcpp::List weights_wrapped;
        std::vector<double> scores_dir_sd;
        scores_dir_sd.reserve(n_PPC);
        std::vector<double> scores_wei_sd;
        scores_wei_sd.reserve(n_PPC);
        for(std::size_t i = 0; i < n_PPC; ++i)
        {
          std::string name_d = "Direction PPC " + std::to_string(i+1);
          std::string name_w = "Weight PPC " + std::to_string(i+1);
          directions_wrapped[name_d] = add_nans_vec(directions.col(i),data_read.second,X.nrow());
          weights_wrapped[name_w] = add_nans_vec(weights.col(i),data_read.second,X.nrow());
          scores_dir_sd.emplace_back(sd_scores_dir_wei[i][0]);
          scores_wei_sd.emplace_back(sd_scores_dir_wei[i][1]);
        }
      
        //saving results in a list, that will be returned
        l["One-step ahead prediction"] = one_step_ahead_pred;
        l["Alpha"]                     = alpha_used;
        l["Number of PPCs retained"]   = n_PPC;
        l["Scores along PPCs"]         = scores_PPC;
        l["Explanatory power PPCs"]    = explanatory_power;
        l["Directions of PPCs"]        = directions_wrapped;
        l["Weights of PPCs"]           = weights_wrapped;
        l["Sd scores directions"]      = scores_dir_sd;
        l["Sd scores weights"]         = scores_wei_sd;
        l["Mean function"]             = mean_func;
        l["Validation errors"]         = errors["Errors"];
      
    }
  }
  
  
  else                                                //VALIDATION ERRORS NOT STORED
  {
    
    if(k>0)                                                                                   //K IMPOSED
    {
      //1D domain, k imposed, not returning errors
      //solver
      auto ko = KO_Factory< SOLVER::ex_solver, K_IMP::YES, VALID_ERR_RET::NO_err, CV_STRAT::AUGMENTING_WINDOW, CV_ERR_EVAL::MSE >::KO_solver(id_CV,std::move(x),alpha,k,threshold_ppc,alphas,k_s,toll,min_dim_train_set,max_dim_train_set,number_threads);
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
      auto sd_scores_dir_wei    = std::get<7>(ko->results());
      auto mean_func            = add_nans_vec(std::get<8>(ko->results()),data_read.second,X.nrow());
      Rcpp::List directions_wrapped;
      Rcpp::List weights_wrapped;
      std::vector<double> scores_dir_sd;
      scores_dir_sd.reserve(n_PPC);
      std::vector<double> scores_wei_sd;
      scores_wei_sd.reserve(n_PPC);
      for(std::size_t i = 0; i < n_PPC; ++i)
      {
        std::string name_d = "Direction PPC " + std::to_string(i+1);
        std::string name_w = "Weight PPC " + std::to_string(i+1);
        directions_wrapped[name_d] = add_nans_vec(directions.col(i),data_read.second,X.nrow());
        weights_wrapped[name_w] = add_nans_vec(weights.col(i),data_read.second,X.nrow());
        scores_dir_sd.emplace_back(sd_scores_dir_wei[i][0]);
        scores_wei_sd.emplace_back(sd_scores_dir_wei[i][1]);
      }

      //saving results in a list, that will be returned
      l["One-step ahead prediction"] = one_step_ahead_pred;
      l["Alpha"]                     = alpha_used;
      l["Number of PPCs retained"]   = n_PPC;
      l["Scores along PPCs"]         = scores_PPC;
      l["Explanatory power PPCs"]    = explanatory_power;
      l["Directions of PPCs"]        = directions_wrapped;
      l["Weights of PPCs"]           = weights_wrapped;
      l["Sd scores directions"]      = scores_dir_sd;
      l["Sd scores weights"]         = scores_wei_sd;
      l["Mean function"]             = mean_func;
    }
    
    else                                                                                      //K NOT IMPOSED
    {
      //1D domain, k not imposed, not returning errors
      //solver
      auto ko = KO_Factory< SOLVER::ex_solver, K_IMP::NO, VALID_ERR_RET::NO_err, CV_STRAT::AUGMENTING_WINDOW, CV_ERR_EVAL::MSE >::KO_solver(id_CV,std::move(x),alpha,k,threshold_ppc,alphas,k_s,toll,min_dim_train_set,max_dim_train_set,number_threads);
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
      auto sd_scores_dir_wei    = std::get<7>(ko->results());
      auto mean_func            = add_nans_vec(std::get<8>(ko->results()),data_read.second,X.nrow());
      Rcpp::List directions_wrapped;
      Rcpp::List weights_wrapped;
      std::vector<double> scores_dir_sd;
      scores_dir_sd.reserve(n_PPC);
      std::vector<double> scores_wei_sd;
      scores_wei_sd.reserve(n_PPC);
      for(std::size_t i = 0; i < n_PPC; ++i)
      {
        std::string name_d = "Direction PPC " + std::to_string(i+1);
        std::string name_w = "Weight PPC " + std::to_string(i+1);
        directions_wrapped[name_d] = add_nans_vec(directions.col(i),data_read.second,X.nrow());
        weights_wrapped[name_w] = add_nans_vec(weights.col(i),data_read.second,X.nrow());
        scores_dir_sd.emplace_back(sd_scores_dir_wei[i][0]);
        scores_wei_sd.emplace_back(sd_scores_dir_wei[i][1]);
      }

      //saving results in a list, that will be returned
      l["One-step ahead prediction"] = one_step_ahead_pred;
      l["Alpha"]                     = alpha_used;
      l["Number of PPCs retained"]   = n_PPC;
      l["Scores along PPCs"]         = scores_PPC;
      l["Explanatory power PPCs"]    = explanatory_power;
      l["Directions of PPCs"]        = directions_wrapped;
      l["Weights of PPCs"]           = weights_wrapped;
      l["Sd scores directions"]      = scores_dir_sd;
      l["Sd scores weights"]         = scores_wei_sd;
      l["Mean function"]             = mean_func;
    }        
  }
  }
  else              //using generalized eigvls problem to retrieve PPCs
  {
  if(err_ret)                                         //VALIDATION ERRORS STORED AND RETURNED
  {
    
    if(k>0)                                                                                   //K IMPOSED
    { 
      //1D domain, k imposed, returning errors
      //solver
      auto ko = KO_Factory< SOLVER::gep_solver, K_IMP::YES, VALID_ERR_RET::YES_err, CV_STRAT::AUGMENTING_WINDOW, CV_ERR_EVAL::MSE >::KO_solver(id_CV,std::move(x),alpha,k,threshold_ppc,alphas,k_s,toll,min_dim_train_set,max_dim_train_set,number_threads);
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
      auto sd_scores_dir_wei    = std::get<7>(ko->results());
      auto mean_func            = add_nans_vec(std::get<8>(ko->results()),data_read.second,X.nrow());
      auto valid_err            = std::get<9>(ko->results());   //valid errors
      Rcpp::List errors = valid_err_disp(valid_err);            //dispatching correctly valid errors
      Rcpp::List directions_wrapped;
      Rcpp::List weights_wrapped;
      std::vector<double> scores_dir_sd;
      scores_dir_sd.reserve(n_PPC);
      std::vector<double> scores_wei_sd;
      scores_wei_sd.reserve(n_PPC);
      for(std::size_t i = 0; i < n_PPC; ++i)
      {
        std::string name_d = "Direction PPC " + std::to_string(i+1);
        std::string name_w = "Weight PPC " + std::to_string(i+1);
        directions_wrapped[name_d] = add_nans_vec(directions.col(i),data_read.second,X.nrow());
        weights_wrapped[name_w] = add_nans_vec(weights.col(i),data_read.second,X.nrow());
        scores_dir_sd.emplace_back(sd_scores_dir_wei[i][0]);
        scores_wei_sd.emplace_back(sd_scores_dir_wei[i][1]);
      }
      
      //saving results in a list, that will be returned
      l["One-step ahead prediction"] = one_step_ahead_pred;
      l["Alpha"]                     = alpha_used;
      l["Number of PPCs retained"]   = n_PPC;
      l["Scores along PPCs"]         = scores_PPC;
      l["Explanatory power PPCs"]    = explanatory_power;
      l["Directions of PPCs"]        = directions_wrapped;
      l["Weights of PPCs"]           = weights_wrapped;
      l["Sd scores directions"]      = scores_dir_sd;
      l["Sd scores weights"]         = scores_wei_sd;
      l["Mean function"]             = mean_func;
      l["Validation errors"]         = errors["Errors"];
    }
    
    else                                                                                      //K NOT IMPOSED
    { 
      //1D domain, k not imposed, returning errors
      //solver
      auto ko = KO_Factory< SOLVER::gep_solver, K_IMP::NO, VALID_ERR_RET::YES_err, CV_STRAT::AUGMENTING_WINDOW, CV_ERR_EVAL::MSE >::KO_solver(id_CV,std::move(x),alpha,k,threshold_ppc,alphas,k_s,toll,min_dim_train_set,max_dim_train_set,number_threads);
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
      auto sd_scores_dir_wei    = std::get<7>(ko->results());
      auto mean_func            = add_nans_vec(std::get<8>(ko->results()),data_read.second,X.nrow());
      auto valid_err            = std::get<9>(ko->results());   //valid errors
      Rcpp::List errors = valid_err_disp(valid_err);            //dispatching correctly valid errors
      Rcpp::List directions_wrapped;
      Rcpp::List weights_wrapped;
      std::vector<double> scores_dir_sd;
      scores_dir_sd.reserve(n_PPC);
      std::vector<double> scores_wei_sd;
      scores_wei_sd.reserve(n_PPC);
      for(std::size_t i = 0; i < n_PPC; ++i)
      {
        std::string name_d = "Direction PPC " + std::to_string(i+1);
        std::string name_w = "Weight PPC " + std::to_string(i+1);
        directions_wrapped[name_d] = add_nans_vec(directions.col(i),data_read.second,X.nrow());
        weights_wrapped[name_w] = add_nans_vec(weights.col(i),data_read.second,X.nrow());
        scores_dir_sd.emplace_back(sd_scores_dir_wei[i][0]);
        scores_wei_sd.emplace_back(sd_scores_dir_wei[i][1]);
      }
      
      //saving results in a list, that will be returned
      l["One-step ahead prediction"] = one_step_ahead_pred;
      l["Alpha"]                     = alpha_used;
      l["Number of PPCs retained"]   = n_PPC;
      l["Scores along PPCs"]         = scores_PPC;
      l["Explanatory power PPCs"]    = explanatory_power;
      l["Directions of PPCs"]        = directions_wrapped;
      l["Weights of PPCs"]           = weights_wrapped;
      l["Sd scores directions"]      = scores_dir_sd;
      l["Sd scores weights"]         = scores_wei_sd;
      l["Mean function"]             = mean_func;
      l["Validation errors"]         = errors["Errors"];
      
    }
  }
  
  
  else                                                //VALIDATION ERRORS NOT STORED
  {
    
    if(k>0)                                                                                   //K IMPOSED
    {
      //1D domain, k imposed, not returning errors
      //solver
      auto ko = KO_Factory< SOLVER::gep_solver, K_IMP::YES, VALID_ERR_RET::NO_err, CV_STRAT::AUGMENTING_WINDOW, CV_ERR_EVAL::MSE >::KO_solver(id_CV,std::move(x),alpha,k,threshold_ppc,alphas,k_s,toll,min_dim_train_set,max_dim_train_set,number_threads);
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
      auto sd_scores_dir_wei    = std::get<7>(ko->results());
      auto mean_func            = add_nans_vec(std::get<8>(ko->results()),data_read.second,X.nrow());
      Rcpp::List directions_wrapped;
      Rcpp::List weights_wrapped;
      std::vector<double> scores_dir_sd;
      scores_dir_sd.reserve(n_PPC);
      std::vector<double> scores_wei_sd;
      scores_wei_sd.reserve(n_PPC);
      for(std::size_t i = 0; i < n_PPC; ++i)
      {
        std::string name_d = "Direction PPC " + std::to_string(i+1);
        std::string name_w = "Weight PPC " + std::to_string(i+1);
        directions_wrapped[name_d] = add_nans_vec(directions.col(i),data_read.second,X.nrow());
        weights_wrapped[name_w] = add_nans_vec(weights.col(i),data_read.second,X.nrow());
        scores_dir_sd.emplace_back(sd_scores_dir_wei[i][0]);
        scores_wei_sd.emplace_back(sd_scores_dir_wei[i][1]);
      }

      //saving results in a list, that will be returned
      l["One-step ahead prediction"] = one_step_ahead_pred;
      l["Alpha"]                     = alpha_used;
      l["Number of PPCs retained"]   = n_PPC;
      l["Scores along PPCs"]         = scores_PPC;
      l["Explanatory power PPCs"]    = explanatory_power;
      l["Directions of PPCs"]        = directions_wrapped;
      l["Weights of PPCs"]           = weights_wrapped;
      l["Sd scores directions"]      = scores_dir_sd;
      l["Sd scores weights"]         = scores_wei_sd;
      l["Mean function"]             = mean_func;
    }
    
    else                                                                                      //K NOT IMPOSED
    {
      //1D domain, k not imposed, not returning errors
      //solver
      auto ko = KO_Factory< SOLVER::gep_solver, K_IMP::NO, VALID_ERR_RET::NO_err, CV_STRAT::AUGMENTING_WINDOW, CV_ERR_EVAL::MSE >::KO_solver(id_CV,std::move(x),alpha,k,threshold_ppc,alphas,k_s,toll,min_dim_train_set,max_dim_train_set,number_threads);
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
      auto sd_scores_dir_wei    = std::get<7>(ko->results());
      auto mean_func            = add_nans_vec(std::get<8>(ko->results()),data_read.second,X.nrow());
      Rcpp::List directions_wrapped;
      Rcpp::List weights_wrapped;
      std::vector<double> scores_dir_sd;
      scores_dir_sd.reserve(n_PPC);
      std::vector<double> scores_wei_sd;
      scores_wei_sd.reserve(n_PPC);
      for(std::size_t i = 0; i < n_PPC; ++i)
      {
        std::string name_d = "Direction PPC " + std::to_string(i+1);
        std::string name_w = "Weight PPC " + std::to_string(i+1);
        directions_wrapped[name_d] = add_nans_vec(directions.col(i),data_read.second,X.nrow());
        weights_wrapped[name_w] = add_nans_vec(weights.col(i),data_read.second,X.nrow());
        scores_dir_sd.emplace_back(sd_scores_dir_wei[i][0]);
        scores_wei_sd.emplace_back(sd_scores_dir_wei[i][1]);
      }

      //saving results in a list, that will be returned
      l["One-step ahead prediction"] = one_step_ahead_pred;
      l["Alpha"]                     = alpha_used;
      l["Number of PPCs retained"]   = n_PPC;
      l["Scores along PPCs"]         = scores_PPC;
      l["Explanatory power PPCs"]    = explanatory_power;
      l["Directions of PPCs"]        = directions_wrapped;
      l["Weights of PPCs"]           = weights_wrapped;
      l["Sd scores directions"]      = scores_dir_sd;
      l["Sd scores weights"]         = scores_wei_sd;
      l["Mean function"]             = mean_func;
    }        
  }
  }

  //return some information useful for plots
  l["Function discrete evaluations points"] = disc_ev_points;
  l["Left extreme domain"]                  = left_extreme; 
  l["Right extreme domain"]                 = right_extreme;
  l["f_n"]                                  = X(_, X.ncol() - 1);
  l["CV"]                                   = id_CV;
  l["Alphas"]                               = alphas;
  l["K_s"]                                  = k_s;
  
  return l;
}







/*!
* @brief Function to perform one-step ahead prediction of Functional Time Series of surfaces using PPCKO.
* @param X Rcpp::NumericMatrix (matrix of double) containing the surface time series: each row (m) is the evaluation of the curve in a point of its domain, each column (n) a time instant
* @param id_CV string denoting which version of the algorithm the user wants to use: 'NoCV', 'CV_alpha', 'CV_k' or 'CV' 
* @param alpha regularization parameter (positive real number)
* @param k number of PPCs: if 0, is selected through explanatory power criterion; if between 1 and m: k is imposed
* @param threshold_ppc minimum requested proportion of explanatory power: used only if k=0. Duble between 0 and 1
* @param alpha_vec vector of double containing the input space for cross-validation over the regularization parameter 
* @param k_vec vector of integers containing the input space for cross-validation over the number of retained PPCs
* @param toll tolerance over which cross-validation process stops looking for bigger number of PPCs if the validation errors between two consecutive iterations is smaller than it
* @param disc_ev_x1 vector of double containing the points of the domain for which the surface's evaluations are available along dimension 1
* @param num_disc_ev_x1 number of evaluation available along dimension 1
* @param disc_ev_x2 vector of double containing the points of the domain for which the surface's evaluations are available along dimension 2
* @param num_disc_ev_x2 number of evaluation available along dimension 2
* @param left_extreme_x1 double indicating the left extreme of the domain of the surface along dimension 1
* @param right_extreme_x1 double indicating the right extreme of the domain of the surface along dimension 1
* @param left_extreme_x2 double indicating the left extreme of the domain of the surface along dimension 2
* @param right_extreme_x2 double indicating the right extreme of the domain of the surface along dimension 2
* @param min_size_ts time instants that define the smallest training set (between 2 and max_size_ts)
* @param max_size_ts time instants that define the biggest training set (between min_size_ts and n-1)
* @param err_ret true if validation errors are stored and returned, false if not
* @param ex_solver true if solving PPCKO inverting the regularized covariance matrix, false if relaying on GEP to avoid id
* @param num_threads number of threads to be used in OMP parallel directives
* @param id_rem_nan string that defines how to handle NaNs for some instant: 'MR': replacing them with the mean of the fts in that point, 'ZR' with 0s
* @return an R list containing:
* - one step ahead prediction of the fts
* - used regularization parameter
* - number of retained PPCs
* - scores along the PPCs (scalar product between the directions and fts at last available instant)
* - cumulative explanatory power of the PPCs (coherent only if ex_solver true)
* - directions of the PPCs
* - weights of the PPCs
* - standard deviation of the scores of each direction (evaluate on the scalar products between direction and the curve between instants 2 and n)
* - standard deviation of the scores of each weight (evaluate on the scalar products between weight and the curve between instants 1 and n-1)
* - mean function of the fts
* - validation errors (only if err_ret true)
* - points of the domain for which the evaluation of the surface are available along dimension 1
* - left extreme of the surface's domain along dimension 1
* - right extreme of the surface's domain along dimension 1
* - points of the domain for which the evaluation of the surface are available along dimension 2
* - left extreme of the surface's domain along dimension 2
* - right extreme of the surface's domain along dimension 2
* - curve at the last available instant
* - which PPCKO version has been performed
* - input space for regularization parameter
* - input space for the number of PPCs
*/
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
                     bool                          err_ret          = false,
                     bool                          ex_solver        = true,
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
  check_solver(ex_solver,id_CV,k);
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
  
  
  //returning element
  Rcpp::List l;
  
  Rcout << "--------------------------------------------------------------------------------------------" << std::endl;
  Rcout << "Running Kargin-Onatski algorithm, " << wrap_string_CV_to_be_printed(id_CV) << std::endl;
  Rcout << "Functional data defined over: [" << left_extreme_x1 << "," << right_extreme_x1 << "] x [" << left_extreme_x2 << "," << right_extreme_x2 <<"], with " << disc_ev_points_x1.size() << " x " << disc_ev_points_x2.size() << " discrete evaluations" << std::endl;
  
  if(ex_solver)   //EX SOLVER
  {
  if(err_ret)                                         //VALIDATION ERRORS STORED AND RETURNED
  {
    
    if(k>0)                                                                                   //K IMPOSED
    {
      //2D domain, k imposed, returning errors
      //solver
      auto ko = KO_Factory< SOLVER::ex_solver, K_IMP::YES, VALID_ERR_RET::YES_err, CV_STRAT::AUGMENTING_WINDOW, CV_ERR_EVAL::MSE >::KO_solver(id_CV,std::move(x),alpha,k,threshold_ppc,alphas,k_s,toll,min_dim_train_set,max_dim_train_set,number_threads);
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
      auto sd_scores_dir_wei    = std::get<7>(ko->results());
      auto mean_func            = from_col_to_matrix(add_nans_vec(std::get<8>(ko->results()),data_read.second,X.nrow()),disc_ev_points_x1.size(),disc_ev_points_x2.size());
      auto valid_err            = std::get<9>(ko->results());   //valid errors
      Rcpp::List errors = valid_err_disp(valid_err);            //dispatching correctly valid errors
      Rcpp::List directions_wrapped;
      Rcpp::List weights_wrapped;
      std::vector<double> scores_dir_sd;
      scores_dir_sd.reserve(n_PPC);
      std::vector<double> scores_wei_sd;
      scores_wei_sd.reserve(n_PPC);
      for(std::size_t i = 0; i < n_PPC; ++i)
      {
        std::string name_d = "Direction PPC " + std::to_string(i+1);
        std::string name_w = "Weight PPC " + std::to_string(i+1);
        directions_wrapped[name_d] = from_col_to_matrix(add_nans_vec(directions.col(i),data_read.second,X.nrow()),disc_ev_points_x1.size(),disc_ev_points_x2.size());
        weights_wrapped[name_w] = from_col_to_matrix(add_nans_vec(weights.col(i),data_read.second,X.nrow()),disc_ev_points_x1.size(),disc_ev_points_x2.size());
        scores_dir_sd.emplace_back(sd_scores_dir_wei[i][0]);
        scores_wei_sd.emplace_back(sd_scores_dir_wei[i][1]);
      }
      
      l["One-step ahead prediction"] = one_step_ahead_pred;
      l["Alpha"]                     = alpha_used;
      l["Number of PPCs retained"]   = n_PPC;
      l["Scores along PPCs"]         = scores_PPC;
      l["Explanatory power PPCs"]    = explanatory_power;
      l["Directions of PPCs"]        = directions_wrapped;
      l["Weights of PPCs"]           = weights_wrapped;
      l["Sd scores directions"]      = scores_dir_sd;
      l["Sd scores weights"]         = scores_wei_sd;
      l["Mean function"]             = mean_func;
      l["Validation errors"]         = errors["Errors"];
    }
    
    else                                                                                      //K NOT IMPOSED
    {
      //2D domain, k not imposed, returning errors
      //solver
      auto ko = KO_Factory< SOLVER::ex_solver, K_IMP::NO, VALID_ERR_RET::YES_err, CV_STRAT::AUGMENTING_WINDOW, CV_ERR_EVAL::MSE >::KO_solver(id_CV,std::move(x),alpha,k,threshold_ppc,alphas,k_s,toll,min_dim_train_set,max_dim_train_set,number_threads);
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
      auto sd_scores_dir_wei    = std::get<7>(ko->results());
      auto mean_func            = from_col_to_matrix(add_nans_vec(std::get<8>(ko->results()),data_read.second,X.nrow()),disc_ev_points_x1.size(),disc_ev_points_x2.size());
      auto valid_err            = std::get<9>(ko->results());   //valid errors
      Rcpp::List errors = valid_err_disp(valid_err);            //dispatching correctly valid errors
      Rcpp::List directions_wrapped;
      Rcpp::List weights_wrapped;
      std::vector<double> scores_dir_sd;
      scores_dir_sd.reserve(n_PPC);
      std::vector<double> scores_wei_sd;
      scores_wei_sd.reserve(n_PPC);
      for(std::size_t i = 0; i < n_PPC; ++i)
      {
        std::string name_d = "Direction PPC " + std::to_string(i+1);
        std::string name_w = "Weight PPC " + std::to_string(i+1);
        directions_wrapped[name_d] = from_col_to_matrix(add_nans_vec(directions.col(i),data_read.second,X.nrow()),disc_ev_points_x1.size(),disc_ev_points_x2.size());
        weights_wrapped[name_w] = from_col_to_matrix(add_nans_vec(weights.col(i),data_read.second,X.nrow()),disc_ev_points_x1.size(),disc_ev_points_x2.size());
        scores_dir_sd.emplace_back(sd_scores_dir_wei[i][0]);
        scores_wei_sd.emplace_back(sd_scores_dir_wei[i][1]);
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
      l["Sd scores directions"]      = scores_dir_sd;
      l["Sd scores weights"]         = scores_wei_sd;
      l["Mean function"]             = mean_func;
      l["Validation errors"]         = errors["Errors"];
    }
  }
  
  
  else                                                //VALIDATION ERRORS NOT STORED
  {
    
    if(k>0)                                                                                   //K IMPOSED
    { 
      //2D domain, k imposed, not returning errors
      //solver
      auto ko = KO_Factory< SOLVER::ex_solver, K_IMP::YES, VALID_ERR_RET::NO_err, CV_STRAT::AUGMENTING_WINDOW, CV_ERR_EVAL::MSE >::KO_solver(id_CV,std::move(x),alpha,k,threshold_ppc,alphas,k_s,toll,min_dim_train_set,max_dim_train_set,number_threads);
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
      auto sd_scores_dir_wei    = std::get<7>(ko->results());
      auto mean_func            = from_col_to_matrix(add_nans_vec(std::get<8>(ko->results()),data_read.second,X.nrow()),disc_ev_points_x1.size(),disc_ev_points_x2.size());
      Rcpp::List directions_wrapped;
      Rcpp::List weights_wrapped;
      std::vector<double> scores_dir_sd;
      scores_dir_sd.reserve(n_PPC);
      std::vector<double> scores_wei_sd;
      scores_wei_sd.reserve(n_PPC);
      for(std::size_t i = 0; i < n_PPC; ++i)
      {
        std::string name_d = "Direction PPC " + std::to_string(i+1);
        std::string name_w = "Weight PPC " + std::to_string(i+1);
        directions_wrapped[name_d] = from_col_to_matrix(add_nans_vec(directions.col(i),data_read.second,X.nrow()),disc_ev_points_x1.size(),disc_ev_points_x2.size());
        weights_wrapped[name_w] = from_col_to_matrix(add_nans_vec(weights.col(i),data_read.second,X.nrow()),disc_ev_points_x1.size(),disc_ev_points_x2.size());
        scores_dir_sd.emplace_back(sd_scores_dir_wei[i][0]);
        scores_wei_sd.emplace_back(sd_scores_dir_wei[i][1]);
      }

      //saving results in a list, that will be returned
      l["One-step ahead prediction"] = one_step_ahead_pred;
      l["Alpha"]                     = alpha_used;
      l["Number of PPCs retained"]   = n_PPC;
      l["Scores along PPCs"]         = scores_PPC;
      l["Explanatory power PPCs"]    = explanatory_power;
      l["Directions of PPCs"]        = directions_wrapped;
      l["Weights of PPCs"]           = weights_wrapped;
      l["Sd scores directions"]      = scores_dir_sd;
      l["Sd scores weights"]         = scores_wei_sd;
      l["Mean function"]             = mean_func;
    }
    
    else                                                                                      //K NOT IMPOSED
    {
      //2D domain, k not imposed, not returning errors
      //solver
      auto ko = KO_Factory< SOLVER::ex_solver, K_IMP::NO, VALID_ERR_RET::NO_err, CV_STRAT::AUGMENTING_WINDOW, CV_ERR_EVAL::MSE >::KO_solver(id_CV,std::move(x),alpha,k,threshold_ppc,alphas,k_s,toll,min_dim_train_set,max_dim_train_set,number_threads);
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
      auto sd_scores_dir_wei    = std::get<7>(ko->results());
      auto mean_func            = from_col_to_matrix(add_nans_vec(std::get<8>(ko->results()),data_read.second,X.nrow()),disc_ev_points_x1.size(),disc_ev_points_x2.size());
      Rcpp::List directions_wrapped;
      Rcpp::List weights_wrapped;
      std::vector<double> scores_dir_sd;
      scores_dir_sd.reserve(n_PPC);
      std::vector<double> scores_wei_sd;
      scores_wei_sd.reserve(n_PPC);
      for(std::size_t i = 0; i < n_PPC; ++i)
      {
        std::string name_d = "Direction PPC " + std::to_string(i+1);
        std::string name_w = "Weight PPC " + std::to_string(i+1);
        directions_wrapped[name_d] = from_col_to_matrix(add_nans_vec(directions.col(i),data_read.second,X.nrow()),disc_ev_points_x1.size(),disc_ev_points_x2.size());
        weights_wrapped[name_w] = from_col_to_matrix(add_nans_vec(weights.col(i),data_read.second,X.nrow()),disc_ev_points_x1.size(),disc_ev_points_x2.size());
        scores_dir_sd.emplace_back(sd_scores_dir_wei[i][0]);
        scores_wei_sd.emplace_back(sd_scores_dir_wei[i][1]);
      }
      
      //saving results in a list, that will be returned
      l["One-step ahead prediction"] = one_step_ahead_pred;
      l["Alpha"]                     = alpha_used;
      l["Number of PPCs retained"]   = n_PPC;
      l["Scores along PPCs"]         = scores_PPC;
      l["Explanatory power PPCs"]    = explanatory_power;
      l["Directions of PPCs"]        = directions_wrapped;
      l["Weights of PPCs"]           = weights_wrapped;
      l["Sd scores directions"]      = scores_dir_sd;
      l["Sd scores weights"]         = scores_wei_sd;
      l["Mean function"]             = mean_func;
    }        
  } 
  }
  else  //GEP
  {
  if(err_ret)                                         //VALIDATION ERRORS STORED AND RETURNED
  {
    
    if(k>0)                                                                                   //K IMPOSED
    {
      //2D domain, k imposed, returning errors
      //solver
      auto ko = KO_Factory< SOLVER::gep_solver, K_IMP::YES, VALID_ERR_RET::YES_err, CV_STRAT::AUGMENTING_WINDOW, CV_ERR_EVAL::MSE >::KO_solver(id_CV,std::move(x),alpha,k,threshold_ppc,alphas,k_s,toll,min_dim_train_set,max_dim_train_set,number_threads);
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
      auto sd_scores_dir_wei    = std::get<7>(ko->results());
      auto mean_func            = from_col_to_matrix(add_nans_vec(std::get<8>(ko->results()),data_read.second,X.nrow()),disc_ev_points_x1.size(),disc_ev_points_x2.size());
      auto valid_err            = std::get<9>(ko->results());   //valid errors
      Rcpp::List errors = valid_err_disp(valid_err);            //dispatching correctly valid errors
      Rcpp::List directions_wrapped;
      Rcpp::List weights_wrapped;
      std::vector<double> scores_dir_sd;
      scores_dir_sd.reserve(n_PPC);
      std::vector<double> scores_wei_sd;
      scores_wei_sd.reserve(n_PPC);
      for(std::size_t i = 0; i < n_PPC; ++i)
      {
        std::string name_d = "Direction PPC " + std::to_string(i+1);
        std::string name_w = "Weight PPC " + std::to_string(i+1);
        directions_wrapped[name_d] = from_col_to_matrix(add_nans_vec(directions.col(i),data_read.second,X.nrow()),disc_ev_points_x1.size(),disc_ev_points_x2.size());
        weights_wrapped[name_w] = from_col_to_matrix(add_nans_vec(weights.col(i),data_read.second,X.nrow()),disc_ev_points_x1.size(),disc_ev_points_x2.size());
        scores_dir_sd.emplace_back(sd_scores_dir_wei[i][0]);
        scores_wei_sd.emplace_back(sd_scores_dir_wei[i][1]);
      }
      
      l["One-step ahead prediction"] = one_step_ahead_pred;
      l["Alpha"]                     = alpha_used;
      l["Number of PPCs retained"]   = n_PPC;
      l["Scores along PPCs"]         = scores_PPC;
      l["Explanatory power PPCs"]    = explanatory_power;
      l["Directions of PPCs"]        = directions_wrapped;
      l["Weights of PPCs"]           = weights_wrapped;
      l["Sd scores directions"]      = scores_dir_sd;
      l["Sd scores weights"]         = scores_wei_sd;
      l["Mean function"]             = mean_func;
      l["Validation errors"]         = errors["Errors"];
    }
    
    else                                                                                      //K NOT IMPOSED
    {
      //2D domain, k not imposed, returning errors
      //solver
      auto ko = KO_Factory< SOLVER::gep_solver, K_IMP::NO, VALID_ERR_RET::YES_err, CV_STRAT::AUGMENTING_WINDOW, CV_ERR_EVAL::MSE >::KO_solver(id_CV,std::move(x),alpha,k,threshold_ppc,alphas,k_s,toll,min_dim_train_set,max_dim_train_set,number_threads);
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
      auto sd_scores_dir_wei    = std::get<7>(ko->results());
      auto mean_func            = from_col_to_matrix(add_nans_vec(std::get<8>(ko->results()),data_read.second,X.nrow()),disc_ev_points_x1.size(),disc_ev_points_x2.size());
      auto valid_err            = std::get<9>(ko->results());   //valid errors
      Rcpp::List errors = valid_err_disp(valid_err);            //dispatching correctly valid errors
      Rcpp::List directions_wrapped;
      Rcpp::List weights_wrapped;
      std::vector<double> scores_dir_sd;
      scores_dir_sd.reserve(n_PPC);
      std::vector<double> scores_wei_sd;
      scores_wei_sd.reserve(n_PPC);
      for(std::size_t i = 0; i < n_PPC; ++i)
      {
        std::string name_d = "Direction PPC " + std::to_string(i+1);
        std::string name_w = "Weight PPC " + std::to_string(i+1);
        directions_wrapped[name_d] = from_col_to_matrix(add_nans_vec(directions.col(i),data_read.second,X.nrow()),disc_ev_points_x1.size(),disc_ev_points_x2.size());
        weights_wrapped[name_w] = from_col_to_matrix(add_nans_vec(weights.col(i),data_read.second,X.nrow()),disc_ev_points_x1.size(),disc_ev_points_x2.size());
        scores_dir_sd.emplace_back(sd_scores_dir_wei[i][0]);
        scores_wei_sd.emplace_back(sd_scores_dir_wei[i][1]);
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
      l["Sd scores directions"]      = scores_dir_sd;
      l["Sd scores weights"]         = scores_wei_sd;
      l["Mean function"]             = mean_func;
      l["Validation errors"]         = errors["Errors"];
    }
  }
  
  
  else                                                //VALIDATION ERRORS NOT STORED
  {
    
    if(k>0)                                                                                   //K IMPOSED
    { 
      //2D domain, k imposed, not returning errors
      //solver
      auto ko = KO_Factory< SOLVER::gep_solver, K_IMP::YES, VALID_ERR_RET::NO_err, CV_STRAT::AUGMENTING_WINDOW, CV_ERR_EVAL::MSE >::KO_solver(id_CV,std::move(x),alpha,k,threshold_ppc,alphas,k_s,toll,min_dim_train_set,max_dim_train_set,number_threads);
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
      auto sd_scores_dir_wei    = std::get<7>(ko->results());
      auto mean_func            = from_col_to_matrix(add_nans_vec(std::get<8>(ko->results()),data_read.second,X.nrow()),disc_ev_points_x1.size(),disc_ev_points_x2.size());
      Rcpp::List directions_wrapped;
      Rcpp::List weights_wrapped;
      std::vector<double> scores_dir_sd;
      scores_dir_sd.reserve(n_PPC);
      std::vector<double> scores_wei_sd;
      scores_wei_sd.reserve(n_PPC);
      for(std::size_t i = 0; i < n_PPC; ++i)
      {
        std::string name_d = "Direction PPC " + std::to_string(i+1);
        std::string name_w = "Weight PPC " + std::to_string(i+1);
        directions_wrapped[name_d] = from_col_to_matrix(add_nans_vec(directions.col(i),data_read.second,X.nrow()),disc_ev_points_x1.size(),disc_ev_points_x2.size());
        weights_wrapped[name_w] = from_col_to_matrix(add_nans_vec(weights.col(i),data_read.second,X.nrow()),disc_ev_points_x1.size(),disc_ev_points_x2.size());
        scores_dir_sd.emplace_back(sd_scores_dir_wei[i][0]);
        scores_wei_sd.emplace_back(sd_scores_dir_wei[i][1]);
      }

      //saving results in a list, that will be returned
      l["One-step ahead prediction"] = one_step_ahead_pred;
      l["Alpha"]                     = alpha_used;
      l["Number of PPCs retained"]   = n_PPC;
      l["Scores along PPCs"]         = scores_PPC;
      l["Explanatory power PPCs"]    = explanatory_power;
      l["Directions of PPCs"]        = directions_wrapped;
      l["Weights of PPCs"]           = weights_wrapped;
      l["Sd scores directions"]      = scores_dir_sd;
      l["Sd scores weights"]         = scores_wei_sd;
      l["Mean function"]             = mean_func;
    }
    
    else                                                                                      //K NOT IMPOSED
    {
      //2D domain, k not imposed, not returning errors
      //solver
      auto ko = KO_Factory< SOLVER::gep_solver, K_IMP::NO, VALID_ERR_RET::NO_err, CV_STRAT::AUGMENTING_WINDOW, CV_ERR_EVAL::MSE >::KO_solver(id_CV,std::move(x),alpha,k,threshold_ppc,alphas,k_s,toll,min_dim_train_set,max_dim_train_set,number_threads);
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
      auto sd_scores_dir_wei    = std::get<7>(ko->results());
      auto mean_func            = from_col_to_matrix(add_nans_vec(std::get<8>(ko->results()),data_read.second,X.nrow()),disc_ev_points_x1.size(),disc_ev_points_x2.size());
      Rcpp::List directions_wrapped;
      Rcpp::List weights_wrapped;
      std::vector<double> scores_dir_sd;
      scores_dir_sd.reserve(n_PPC);
      std::vector<double> scores_wei_sd;
      scores_wei_sd.reserve(n_PPC);
      for(std::size_t i = 0; i < n_PPC; ++i)
      {
        std::string name_d = "Direction PPC " + std::to_string(i+1);
        std::string name_w = "Weight PPC " + std::to_string(i+1);
        directions_wrapped[name_d] = from_col_to_matrix(add_nans_vec(directions.col(i),data_read.second,X.nrow()),disc_ev_points_x1.size(),disc_ev_points_x2.size());
        weights_wrapped[name_w] = from_col_to_matrix(add_nans_vec(weights.col(i),data_read.second,X.nrow()),disc_ev_points_x1.size(),disc_ev_points_x2.size());
        scores_dir_sd.emplace_back(sd_scores_dir_wei[i][0]);
        scores_wei_sd.emplace_back(sd_scores_dir_wei[i][1]);
      }
      
      //saving results in a list, that will be returned
      l["One-step ahead prediction"] = one_step_ahead_pred;
      l["Alpha"]                     = alpha_used;
      l["Number of PPCs retained"]   = n_PPC;
      l["Scores along PPCs"]         = scores_PPC;
      l["Explanatory power PPCs"]    = explanatory_power;
      l["Directions of PPCs"]        = directions_wrapped;
      l["Weights of PPCs"]           = weights_wrapped;
      l["Sd scores directions"]      = scores_dir_sd;
      l["Sd scores weights"]         = scores_wei_sd;
      l["Mean function"]             = mean_func;
    }        
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


/*!
* @brief Function to perform pointwise ADF-test p-values for curve fts
* @param X Rcpp::NumericMatrix (matrix of double) containing the curve time series: each row (m) is the evaluation of the curve in a point of its domain, each column (n) a time instant
* @return an R list containing the pointiwise p-values
*/
//
// [[Rcpp::export]]
Rcpp::List KO_check_hps(Rcpp::NumericMatrix X)
{
  using T = double; 
  
  Rcout << "--------------------------------------" << std::endl;
  Rcout << "Pointwise ADF test p-values evaluation" << std::endl;
  
  //read data, handle NaNs
  auto data_read = reader_data<T>(X, REM_NAN::MR);
  KO_Traits::StoringMatrix x = data_read.first;
  
  int number_time_instants = x.cols();
  
  //to check if lag orders bigger than ones have to be taken into account
  std::size_t k = static_cast<std::size_t>(std::trunc(std::cbrt(static_cast<double>(number_time_instants)-1)));
  
  //returning element
  Rcpp::List l;
  
  if(k > 1)
  {
    //ADF with lag order k
    adf<CaseLagOrderADF> adf_t(std::move(x),k);
    adf_t.test();                   //test
    auto pv = adf_t.p_values();     //pvalues
    auto pv_final = add_nans_vec(Eigen::Map<const KO_Traits::StoringVector>(pv.data(),pv.size()),data_read.second,X.nrow());  //preparing the element to be returned
    
    l["P-values ADF"] = pv_final;
    return l;
  }
  
  //ADF without lag order bigger than one
  adf<CaseLagOrderADF> adf_t(std::move(x),0);
  adf_t.test();                 //test
  auto pv = adf_t.p_values();   //pvalues
  auto pv_final = add_nans_vec(Eigen::Map<const KO_Traits::StoringVector>(pv.data(),pv.size()),data_read.second,X.nrow());    //preparing the element to be returned
  
  //returning element
  l["P-values ADF"] = pv_final;
  return l;
}


/*!
* @brief Function to perform pointwise ADF-test p-values for surface fts
* @param X Rcpp::NumericMatrix (matrix of double) containing the surface time series: each row (m) is the evaluation of the curve in a point of its domain, each column (n) a time instant
* @param dim_x1 number of surface's evaluations available along dimension 1
* @param dim_x2 number of surface's evaluations available along dimension 2
* @return an R list containing the pointiwise p-values
*/
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


/*!
* @brief Function to map an R list of matrices into a coherent matrix for PPCKO_2d
* @param Xt an R list of matrices such that each one represents the surface at a given instant, increasingly ordered
* @return Rcpp::NumericMatrix (matrix of double) containing the surface time series: each row (m) is the evaluation of the curve in a point of its domain, each column (n) a time instant
*/
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
  
  //number of point evaluation for the surface (has to be the same for every instant) (dummy NaNs have to be included)
  int number_point_evaluations = as<NumericMatrix>(Xt[0]).nrow()*as<NumericMatrix>(Xt[0]).ncol();
  if(number_point_evaluations==0)
  {
    std::string error_message2 = "List of empty matrices";
    throw std::invalid_argument(error_message2);
  }
  
  //matrix to be returned
  Rcpp::NumericMatrix x(number_point_evaluations,number_time_instants);
  
  for(std::size_t i = 0; i < static_cast<std::size_t>(number_time_instants); ++i)
  {
    //read each element of the list
    KO_Traits::StoringMatrix col = reader_data<T>(as<NumericMatrix>(Xt[i]),REM_NAN::NR).first;
    //map the matrix to a column vector
    KO_Traits::StoringVector first_col = from_matrix_to_col(col);
    //wrap to an R object
    Rcpp::NumericVector col_wrapped = Rcpp::wrap(first_col);
    //adding it as column
    x(_,i) = col_wrapped;
  }
   
  return x;
}


/*!
* @brief Function to map an R array into a coherent matrix for PPCKO_2d
* @param Xt an R array such that element [i,j,k] represents the surface in the evaluation (x1_i,x2_j) at instant k
* @return Rcpp::NumericMatrix (matrix of double) containing the surface time series: each row (m) is the evaluation of the curve in a point of its domain, each column (n) a time instant
*/
//
// [[Rcpp::export]]
Rcpp::NumericMatrix data_2d_wrapper_from_array(Rcpp::NumericVector Xt)
{
  using T = double;
  
  //obtaining the dimensions from the array
  IntegerVector dimensions = Xt.attr("dim");
  int dim1_size = dimensions[0];
  int dim2_size = dimensions[1];
  int number_point_evaluations = dim1_size*dim2_size;
  int number_time_instants = dimensions[2];   //this works only for 1-step time series
  
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
    
    //the mapping is forced to be done the hard way
    for (int i1 = 0; i1 < dim1_size; ++i1) {
      for (int i2 = 0; i2 < dim2_size; ++i2) {
        instant_data(i1, i2) = Xt[i1 + i2*dim1_size + i*dim1_size*dim2_size];  
      }
    }
    
    //wrapping data into R object
    KO_Traits::StoringMatrix col = reader_data<T>(instant_data,REM_NAN::NR).first;
    KO_Traits::StoringVector first_col = from_matrix_to_col(col);
    Rcpp::NumericVector col_wrapped = Rcpp::wrap(first_col);
    //add it as column
    x(_,i) = col_wrapped;
  }
  
  return x;
}