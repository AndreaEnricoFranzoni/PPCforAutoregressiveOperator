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

#ifndef KO_PPC_CV_ALPHA_CRTP_HPP
#define KO_PPC_CV_ALPHA_CRTP_HPP

#include "PPC_KO.hpp"
#include "PPC_KO_NoCV.hpp"


/*!
* @file PPC_KO_CV_alpha.hpp
* @brief Class for computing PPCKO algortihm with cross-validation on the regularization parameter: k can be a parameter or selected through explanatory power criterion
* @author Andrea Enrico Franzoni
*/



/*!
* @class PPC_KO_CV_alpha
* @brief Derived from 'PPC_KO_base' class for computing PPCKO algorithm with cross-validation on regularization parameter
* @tparam solver if algorithm solved inverting the regularized covariance or avoiding it through gep (not possible if retaining the number of PPCs with explanatory power criterion)
* @tparam k_imp if k is imposed or has to be found through explanatory power criterion
* @tparam valid_err_ret if validation error are stored
* @tparam cv_strat strategy for splitting training/validation sets
* @tparam err_eval how to evaluate the loss between prediction on validation set and validation set
*/
template< SOLVER solver, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval >
class PPC_KO_CV_alpha : public PPC_KO_base<PPC_KO_CV_alpha<solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>,solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>
{
private:

  /*!Input space for regularization parameter*/
  std::vector<double> m_alphas;
  /*!Fts not centered*/
  KO_Traits::StoringMatrix m_X_non_cent;
  /*!Smallest training set size (number of time instants)*/
  int m_min_size_ts;
  /*!Biggest training set size (number of time instants)*/
  int m_max_size_ts;  
  
  
public:

  /*!
  * @brief Constructor for cv version on regularization parameter if k is passed as parameter
  * @param X fts
  * @param alphas regularization parameter input space
  * @param k number of retained PPCs
  * @param min_size_ts smallest training set size (number of time instants)
  * @param max_size_ts biggest training set size (number of time instants)
  * @param number_threads number of threads for OMP
  * @details Universal constructor: move semantic used to optimazing handling big size objects
  * @note eventual usage of 'pragma' directive for OMP
  */
  template<typename STOR_OBJ>
  PPC_KO_CV_alpha(STOR_OBJ&& X, const std::vector<double> &alphas, int k, int min_size_ts, int max_size_ts, int number_threads) 
    : 
    PPC_KO_base<PPC_KO_CV_alpha,solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>(std::move(X),number_threads),
    m_alphas(alphas),
    m_X_non_cent(this->m(),this->n()),
    m_min_size_ts(min_size_ts),
    m_max_size_ts(max_size_ts)
    {
      this->k() = k; 
      
      //decentering the data, to pass them in the various cv iterations not centered
#ifdef _OPENMP
#pragma omp parallel for num_threads(this->number_threads())
#endif
      for (size_t i = 0; i < this->n(); ++i)
      {
        m_X_non_cent.col(i) = this->X().col(i).array() + this->means();
      }
    }
  
  /*!
  * @brief Constructor for cv version on regularization parameter if k has to be selected through requested explanatory power on retained PPCs
  * @param X fts 
  * @param alphas regularization parameter input space
  * @param threshold_ppc requested explanatory power for the retained PPCs
  * @param min_size_ts smallest training set size (number of time instants)
  * @param max_size_ts biggest training set size (number of time instants)
  * @param number_threads number of threads for OMP
  * @details Universal constructor: move semantic used to optimazing handling big size objects
  * @note eventual usage of 'pragma' directive for OMP
  */
  template<typename STOR_OBJ>
  PPC_KO_CV_alpha(STOR_OBJ&& X, const std::vector<double> &alphas, double threshold_ppc, int min_size_ts, int max_size_ts, int number_threads) 
    : 
    PPC_KO_base<PPC_KO_CV_alpha,solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>(std::move(X),number_threads),
    m_alphas(alphas),
    m_X_non_cent(this->m(),this->n()),
    m_min_size_ts(min_size_ts),
    m_max_size_ts(max_size_ts)
    {
      this->threshold_ppc() = threshold_ppc; 
      
      //decentering the data, to pass them in the various cv iterations not centered
#ifdef _OPENMP
#pragma omp parallel for num_threads(this->number_threads())
#endif
      for (size_t i = 0; i < this->n(); ++i)
      {
        m_X_non_cent.col(i) = this->X().col(i).array() + this->means();
      }
    }
  
  /*!
  * @brief Method to perform PPCKO if cv is performed on regularization parameter
  * @details Selects the regularization parameter through cv, and then call the .KO_algo() method of the base class  
  */
  inline
  void 
  solving()
  {
    //factory to create the cv strategy
    auto strategy_cv = Factory_cv_strat<cv_strat>::cv_strat_obj(m_min_size_ts,m_max_size_ts);

    // if k imposed
    if constexpr(k_imp == K_IMP::YES)
    {
      //lambda wrapper for the correct overload for prediction function
      auto predictor = [](KO_Traits::StoringMatrix&& data, double alpha, int k, int number_threads) { return cv_pred_func<solver,K_IMP::YES,VALID_ERR_RET::NO_err,cv_strat,cv_err_eval>(std::move(data),alpha,k,number_threads);};
      
      //cv knowing k
      CV_alpha<cv_strat,cv_err_eval,k_imp,valid_err_ret> cv(std::move(m_X_non_cent),std::move(*strategy_cv),m_alphas,this->k(),predictor,this->number_threads());
      
      //best alpha
      cv.best_param_search();
      this->alpha() = cv.param_best();      //finding the best alpha using CV
      
      //if errors to be saved
      if constexpr(valid_err_ret == VALID_ERR_RET::YES_err)
      {
        valid_err_cv_1_t m;
        m.reserve(m_alphas.size());
        for (size_t i = 0; i < m_alphas.size(); ++i)  { m.emplace_back(cv.valid_errors()[i]);}
        this->ValidErr() = m;
        m.clear();
      }
    }
    
    // if k selected through explanatory power criterion
    if constexpr(k_imp == K_IMP::NO)
    {
      //lambda wrapper for the correct overload for prediction function
      auto predictor = [](KO_Traits::StoringMatrix&& data, double alpha, double threshold_ppc, int number_threads) { return cv_pred_func<solver,K_IMP::NO,VALID_ERR_RET::NO_err,cv_strat,cv_err_eval>(std::move(data),alpha,threshold_ppc,number_threads);};
      
      //cv with k to be found with explanatory power
      CV_alpha<cv_strat,cv_err_eval,k_imp,valid_err_ret> cv(std::move(m_X_non_cent),std::move(*strategy_cv),m_alphas,this->threshold_ppc(),predictor,this->number_threads());
      
      //best alpha
      cv.best_param_search();
      this->alpha() = cv.param_best();      //finding the best alpha using CV
      
      //only if errors are saved
      if constexpr(valid_err_ret == VALID_ERR_RET::YES_err)
      {
        valid_err_cv_1_t m;
        m.reserve(m_alphas.size());
        for (size_t i = 0; i < m_alphas.size(); ++i)  { m.emplace_back(cv.valid_errors()[i]);}
        this->ValidErr() = m;
        m.clear();
      }
    }
    
    //computing regularized covariance
    this->CovReg() = this->Cov().array() + this->alpha()*this->trace_cov()*(KO_Traits::StoringMatrix::Identity(this->m(),this->m()).array());
    //PPCKO
    this->KO_algo(); 
  };
};

#endif  //KO_PPC_CV_ALPHA_CRTP_HPP