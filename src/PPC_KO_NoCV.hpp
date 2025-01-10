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

#ifndef KO_PPC_NOCV_CRTP_HPP
#define KO_PPC_NOCV_CRTP_HPP

#include "PPC_KO.hpp"


/*!
* @file PPC_KO_NoCV.hpp
* @brief Class for computing PPCKO algortihm without cross-validation: alpha as parameter, k can be a parameter or selected through explanatory power criterion
* @author Andrea Enrico Franzoni
* @note definition of the function to make prediction on the validation set during cv process is done here: they simply wrap a PPCKO without cross-validation
*/



/*!
* @class PPC_KO_NoCV
* @brief Derived from 'PPC_KO_base' class for computing PPCKO algorithm without cross-validation
* @tparam solver if algorithm solved inverting the regularized covariance or avoiding it through gep (not possible if retaining the number of PPCs with explanatory power criterion)
* @tparam k_imp if k is imposed or has to be found through explanatory power criterion
* @tparam valid_err_ret if validation error are stored
* @tparam cv_strat strategy for splitting training/validation sets
* @tparam err_eval how to evaluate the loss between prediction on validation set and validation set
*/
template< SOLVER solver, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval > 
class PPC_KO_NoCV : public PPC_KO_base<PPC_KO_NoCV<solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>,solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>
{
public:
  
  /*!
  * @brief Constructor for no cv version if k is passed as parameter
  * @param X fts
  * @param alpha regularization parameter
  * @param k number of retained PPCs
  * @param number_threads number of threads for OMP
  * @details Universal constructor: move semantic used to optimazing handling big size objects
  */
  template<typename STOR_OBJ>
  PPC_KO_NoCV(STOR_OBJ&& X, double alpha, int k, int number_threads)
    :   PPC_KO_base<PPC_KO_NoCV,solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>(std::move(X),number_threads)
    { 
      //saving parameters in the base class
      this->alpha() = alpha;
      this->k() = k;
      //computing the regularized covariance
      this->CovReg() = this->Cov().array() + this->alpha()*this->trace_cov()*(KO_Traits::StoringMatrix::Identity(this->m(),this->m()).array());
    }
  
  /*!
  * @brief Constructor for no cv version if k is selected through explanatory power criterion
  * @param X fts
  * @param alpha regularization parameter
  * @param threshold_ppc requested explanatory power of the retained PPCs
  * @param number_threads number of threads for OMP
  * @details Universal constructor: move semantic used to optimazing handling big size objects
  */
  template<typename STOR_OBJ>
  PPC_KO_NoCV(STOR_OBJ&& X, double alpha, double threshold_ppc, int number_threads)
    :   PPC_KO_base<PPC_KO_NoCV,solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>(std::move(X),number_threads)
    {
      //saving parameters in the base class
      this->alpha() = alpha;
      this->threshold_ppc() = threshold_ppc;
      //computing the regularized covariance
      this->CovReg() = this->Cov().array() + this->alpha()*this->trace_cov()*(KO_Traits::StoringMatrix::Identity(this->m(),this->m()).array());
    }
  
  /*!
  * @brief Method to perform PPCKO if no cv is performed
  * @details Wraps the .KO_algo() method of the base class since parameters are known
  */
  inline 
  void
  solving()
  {
    //PPCKO
    this->KO_algo(); 
  }
};





/*!
* @brief Function to make prediction on the validation set during cross-validation process is k is imposed (by the user or by cv process)
* @tparam solver if algorithm solved inverting the regularized covariance or avoiding it through gep (not possible if retaining the number of PPCs with explanatory power criterion)
* @tparam k_imp if k is imposed or has to be found through explanatory power criterion
* @tparam valid_err_ret if validation error are stored
* @tparam cv_strat strategy for splitting training/validation sets
* @tparam err_eval how to evaluate the loss between prediction on validation set and validation set
* @param training_data training set (not centered)
* @param alpha regularization parameter
* @param k number of retained PPCs
* @param number_threads number of threads for OMP
* @return The prediction on the validation set
* @details It creates a 'PPC_KO_NoCV' object with the given parameters, trains it with 'training_data' and makes prediction
*/
template< SOLVER solver, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval >
KO_Traits::StoringVector 
cv_pred_func(KO_Traits::StoringMatrix && training_data, double alpha, int k, int number_threads)
{  
  //k imposed: template parameter fixed
  PPC_KO_NoCV<solver,K_IMP::YES,valid_err_ret,cv_strat,cv_err_eval> iter(std::move(training_data),alpha,k,number_threads);
  iter.solving();
  
  return iter.prediction(); 
};


/*!
* @brief Function to make prediction on the validation set during cross-validation process is k is selected through explanatory power criterion
* @tparam solver if algorithm solved inverting the regularized covariance or avoiding it through gep (not possible if retaining the number of PPCs with explanatory power criterion)
* @tparam k_imp if k is imposed or has to be found through explanatory power criterion
* @tparam valid_err_ret if validation error are stored
* @tparam cv_strat strategy for splitting training/validation sets
* @tparam err_eval how to evaluate the loss between prediction on validation set and validation set
* @param training_data training set (not centered)
* @param alpha regularization parameter
* @param threshold_ppc requested explanatory power by the retained PPCs
* @param number_threads number of threads for OMP
* @return The prediction on the validation set
* @details It creates a 'PPC_KO_NoCV' object with the given parameters, trains it with 'training_data' and makes prediction
*/
template< SOLVER solver, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval >
KO_Traits::StoringVector 
cv_pred_func(KO_Traits::StoringMatrix && training_data, double alpha, double threshold_ppc, int number_threads)
{  
  //k not imposed: template parameter fixed
  PPC_KO_NoCV<solver,K_IMP::NO,valid_err_ret,cv_strat,cv_err_eval> iter(std::move(training_data),alpha,threshold_ppc,number_threads);
  iter.solving();
  
  return iter.prediction(); 
};

#endif  //KO_PPC_NOCV_CRTP_HPP