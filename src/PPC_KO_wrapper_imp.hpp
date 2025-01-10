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


#include "PPC_KO_wrapper.hpp"

#include <iostream>
#include <algorithm>

/*!
* @file PPC_KO_wrapper_imp.hpp
* @brief Implementation of the overriding of the virtual base method, depending on the child class
* @author Andrea Enrico Franzoni
*/


/*!
* @brief No cross-validation overriding: wraps the class for computations (static polymorphism) accordingly
* @details Wraps the class for computations, performs them and then update the results in the wrapper class
*/
template< SOLVER solver, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval >
void
PPC_KO_wrapper_no_cv<solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>::call_ko()
{
  
  if constexpr(k_imp == K_IMP::YES)   //k imposed
  {
    //class for computations construction
    PPC_KO_NoCV<solver,K_IMP::YES,valid_err_ret,cv_strat,cv_err_eval> KO(std::move(this->data()),m_alpha,m_k,this->number_threads());
    //solving
    KO.solve(); 
    //computing scores
    auto scores = KO.scores();
    //computing sd of scores of directions and weights
    auto sd_scores = KO.sd_scores_dir_wei();

    //if validation errors have to be stored and returned
    if constexpr( valid_err_ret == VALID_ERR_RET::YES_err){this->results() = std::make_tuple(KO.prediction(),KO.alpha(),KO.k(),scores,KO.explanatory_power(),KO.a(),KO.b(),sd_scores,KO.means(),KO.ValidErr());}
    else  {this->results() = std::make_tuple(KO.prediction(),KO.alpha(),KO.k(),scores,KO.explanatory_power(),KO.a(),KO.b(),sd_scores,KO.means());}
  }

  if constexpr(k_imp == K_IMP::NO)    //k to be found with explanatory power criterion
  {
    //class for computations construction
    PPC_KO_NoCV<solver,K_IMP::NO,valid_err_ret,cv_strat,cv_err_eval> KO(std::move(this->data()),m_alpha,m_threshold_ppc,this->number_threads());
    //solving
    KO.solve();
    //computing scores
    auto scores = KO.scores();
    //computing sd of scores of directions and weights
    auto sd_scores = KO.sd_scores_dir_wei();
    
    //if validation errors have to be stored and returned
    if constexpr( valid_err_ret == VALID_ERR_RET::YES_err){this->results() = std::make_tuple(KO.prediction(),KO.alpha(),KO.k(),scores,KO.explanatory_power(),KO.a(),KO.b(),sd_scores,KO.means(),KO.ValidErr());}
    else  {this->results() = std::make_tuple(KO.prediction(),KO.alpha(),KO.k(),scores,KO.explanatory_power(),KO.a(),KO.b(),sd_scores,KO.means());}
  }
}
  
  

/*!
* @brief Regularization parameter cross-validation overriding: wraps the class for computations (static polymorphism) accordingly
* @details Wraps the class for computations, performs them and then update the results in the wrapper class
*/
template< SOLVER solver, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval >
void
PPC_KO_wrapper_cv_alpha<solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>::call_ko()
{
  
  if constexpr(k_imp == K_IMP::YES)   //k imposed
  {
    //class for computations construction
    PPC_KO_CV_alpha<solver,K_IMP::YES,valid_err_ret,cv_strat,cv_err_eval> KO(std::move(this->data()),m_alphas,m_k,m_min_size_ts,m_max_size_ts,this->number_threads());
    //solving
    KO.solve();
    //computing scores
    auto scores = KO.scores();
    //computing sd of scores of directions and weights
    auto sd_scores = KO.sd_scores_dir_wei();
    
    //if validation errors have to be stored and returned
    if constexpr( valid_err_ret == VALID_ERR_RET::YES_err){this->results() = std::make_tuple(KO.prediction(),KO.alpha(),KO.k(),scores,KO.explanatory_power(),KO.a(),KO.b(),sd_scores,KO.means(),std::get<valid_err_cv_1_t>(KO.ValidErr()));}
    else  {this->results() = std::make_tuple(KO.prediction(),KO.alpha(),KO.k(),scores,KO.explanatory_power(),KO.a(),KO.b(),sd_scores,KO.means());}
  }

  if constexpr(k_imp == K_IMP::NO)    //k to be found with explanatory power criterion
  {
    //class for computations construction
    PPC_KO_CV_alpha<solver,K_IMP::NO,valid_err_ret,cv_strat,cv_err_eval> KO(std::move(this->data()),m_alphas,m_threshold_ppc,m_min_size_ts,m_max_size_ts,this->number_threads());
    //solving
    KO.solve();
    //computing scores
    auto scores = KO.scores();
    //computing sd of scores of directions and weights
    auto sd_scores = KO.sd_scores_dir_wei();
    
    //if validation errors have to be stored and returned
    if constexpr( valid_err_ret == VALID_ERR_RET::YES_err){this->results() = std::make_tuple(KO.prediction(),KO.alpha(),KO.k(),scores,KO.explanatory_power(),KO.a(),KO.b(),sd_scores,KO.means(),std::get<valid_err_cv_1_t>(KO.ValidErr()));}
    else  {this->results() = std::make_tuple(KO.prediction(),KO.alpha(),KO.k(),scores,KO.explanatory_power(),KO.a(),KO.b(),sd_scores,KO.means());}  
  }
}



/*!
* @brief Number of retained PPCs cross-validation overriding: wraps the class for computations (static polymorphism) accordingly
* @details Wraps the class for computations, performs them and then update the results in the wrapper class
*/
template< SOLVER solver, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval >
void
PPC_KO_wrapper_cv_k<solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>::call_ko()
{
  //class for computations construction (k_imp=K_IMP::YES by default)
  PPC_KO_CV_k<solver,K_IMP::YES,valid_err_ret,cv_strat,cv_err_eval> KO(std::move(this->data()),m_k_s,m_alpha,m_toll,m_min_size_ts,m_max_size_ts,this->number_threads());
  //solving
  KO.solve();
  //computing scores
  auto scores = KO.scores();
  //computing sd of scores of directions and weights
  auto sd_scores = KO.sd_scores_dir_wei();
  
  //if validation errors have to be stored and returned
  if constexpr( valid_err_ret == VALID_ERR_RET::YES_err){this->results() = std::make_tuple(KO.prediction(),KO.alpha(),KO.k(),scores,KO.explanatory_power(),KO.a(),KO.b(),sd_scores,KO.means(),std::get<valid_err_cv_1_t>(KO.ValidErr()));}
  else  {this->results() = std::make_tuple(KO.prediction(),KO.alpha(),KO.k(),scores,KO.explanatory_power(),KO.a(),KO.b(),sd_scores,KO.means());}
}



/*!
* @brief Regularization parameter and number of retained PPCs cross-validation overriding: wraps the class for computations (static polymorphism) accordingly
* @details Wraps the class for computations, performs them and then update the results in the wrapper class
*/
template< SOLVER solver, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval >
void
PPC_KO_wrapper_cv_alpha_k<solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>::call_ko()
{
  //class for computations construction (k_imp=K_IMP::YES by default)
  PPC_KO_CV_alpha_k<solver,K_IMP::YES,valid_err_ret,cv_strat,cv_err_eval> KO(std::move(this->data()),m_alphas,m_k_s,m_toll,m_min_size_ts,m_max_size_ts,this->number_threads());
  //solving
  KO.solve();
  //computing scores
  auto scores = KO.scores();
  //computing sd of scores of directions and weights
  auto sd_scores = KO.sd_scores_dir_wei();
  
  //if validation errors have to be stored and returned
  if constexpr( valid_err_ret == VALID_ERR_RET::YES_err){this->results() = std::make_tuple(KO.prediction(),KO.alpha(),KO.k(),scores,KO.explanatory_power(),KO.a(),KO.b(),sd_scores,KO.means(),std::get<valid_err_cv_2_t>(KO.ValidErr()));}
  else  {this->results() = std::make_tuple(KO.prediction(),KO.alpha(),KO.k(),scores,KO.explanatory_power(),KO.a(),KO.b(),sd_scores,KO.means());}
}