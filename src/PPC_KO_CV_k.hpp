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

#ifndef KO_PPC_CV_K_CRTP_HPP
#define KO_PPC_CV_K_CRTP_HPP

#include "PPC_KO.hpp"
#include "PPC_KO_NoCV.hpp"


/*!
* @file PPC_KO_CV_k.hpp
* @brief Class for computing PPCKO algortihm with cross-validation on the number of retained PPCs: regularization parameter is a passed parameter
* @author Andrea Enrico Franzoni
*/



/*!
* @class PPC_KO_CV_k
* @brief Derived from 'PPC_KO_base' class for computing PPCKO algorithm with cross-validation on the number of retained PPCs
* @tparam solver if algorithm solved inverting the regularized covariance or avoiding it through gep (not possible if retaining the number of PPCs with explanatory power criterion)
* @tparam k_imp if k is imposed or has to be found through explanatory power criterion
* @tparam valid_err_ret if validation error are stored
* @tparam cv_strat strategy for splitting training/validation sets
* @tparam err_eval how to evaluate the loss between prediction on validation set and validation set
*/
template< SOLVER solver, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval >
class PPC_KO_CV_k : public PPC_KO_base<PPC_KO_CV_k<solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>,solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>
{
private:

  /*!Input space for the number of retained PPCs*/
  std::vector<int> m_k_s;
  /*!Fts not centered*/
  KO_Traits::StoringMatrix m_X_non_cent;
  /*!Tolerance: the cv continues only if between two parameters, that are checked in increasing order, 
  * the absolute difference between two validation errors is bigger than tolerance*trace(covariance). 
  * If not, stops and look for k only between the tested ones 
  */  
  double m_toll;
  /*!Smallest training set size (number of time instants)*/
  int m_min_size_ts;
  /*!Biggest training set size (number of time instants)*/
  int m_max_size_ts;
  
public:
  
  /*!
  * @brief Constructor for cv version on number of retained PPCs
  * @param X fts
  * @param k_s number of retained PPCs input space
  * @param alpha regularization parameter
  * @param toll the cv on the number of retained PPCs continues only if between two parameters, that are checked in increasing order, the absolute difference between two validation errors is bigger than tolerance*trace(covariance). If not, stops and look for k only between the tested ones 
  * @param min_size_ts smallest training set size (number of time instants)
  * @param max_size_ts biggest training set size (number of time instants)
  * @param number_threads number of threads for OMP
  * @details Universal constructor: move semantic used to optimazing handling big size objects
  * @note eventual usage of 'pragma' directive for OMP
  */
  template<typename STOR_OBJ>
  PPC_KO_CV_k(STOR_OBJ&& X, std::vector<int> &k_s, double alpha, double toll, int min_size_ts, int max_size_ts, int number_threads) 
    : 
    PPC_KO_base<PPC_KO_CV_k,solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>(std::move(X),number_threads),
    m_k_s(k_s),
    m_X_non_cent(this->m(),this->n()),
    m_toll(toll),
    m_min_size_ts(min_size_ts),
    m_max_size_ts(max_size_ts)
    {
      //computing regularized covariance
      this->alpha() = alpha;
      this->CovReg() = this->Cov().array() + this->alpha()*this->trace_cov()*(KO_Traits::StoringMatrix::Identity(this->m(),this->m()).array());
      
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
  * @brief Method to perform PPCKO if cv is performed on the number of retained PPCs (k will be always imposed)
  * @details Selects the number of PPCs through cv, and then call the .KO_algo() method of the base class  
  */
  inline
  void 
  solving()
  {
    //factory to create the cv strategy
    auto strategy_cv = Factory_cv_strat<cv_strat>::cv_strat_obj(m_min_size_ts,m_max_size_ts);

    //lambda wrapper for the correct overload for prediction function
    auto predictor = [](KO_Traits::StoringMatrix&& data, double alpha, int k, int number_threads) { return cv_pred_func<solver,K_IMP::YES,VALID_ERR_RET::NO_err,cv_strat,cv_err_eval>(std::move(data),alpha,k,number_threads);};
    
    //to stop the algorithm if adding a PPC useless
    double toll_param = m_toll*this->trace_cov();
    
    //cv for k
    CV_k<cv_strat,cv_err_eval,K_IMP::YES,valid_err_ret> cv(std::move(m_X_non_cent),std::move(*strategy_cv),m_k_s,toll_param,this->alpha(),predictor,this->number_threads());
    
    //best number of PPCs
    cv.best_param_search();
    this->k() = cv.param_best();
    
    //if errors to be saved
    if constexpr (valid_err_ret == VALID_ERR_RET::YES_err)
    {
      valid_err_cv_1_t m;
      m.reserve(cv.valid_errors().size());
      for (size_t i = 0; i < cv.valid_errors().size(); ++i)  { m.emplace_back(cv.valid_errors()[i]);}
      this->ValidErr() = m;
      m.clear();
    }
    
    //PPCKO
    this->KO_algo(); 
  };
};

#endif  //KO_PPC_CV_K_CRTP_HPP