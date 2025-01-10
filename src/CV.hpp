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

#ifndef CV_CRTP_PPC_HPP
#define CV_CRTP_PPC_HPP

#include <iostream>
#include <functional>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <concepts>
#include <memory>
#include <utility>
#include <type_traits>

#include "traits_ko.hpp"
#include "strategy_cv.hpp"
#include "cv_eval_valid_err.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

/*!
* @file CV.hpp
* @brief Class for performing cross-validation
* @author Andrea Enrico Franzoni
*/


/*!
* Doing tag dispatching for the correct way of evaluating the error between prediction on validation set and validation set.
* @tparam err_eval: template parameter for the error evaluation strategy
*/
template <CV_ERR_EVAL err_eval>
using ERR_EVAL_T = std::integral_constant<CV_ERR_EVAL, err_eval>;


/*!Type for the function that predicts the validation set: k imposed: pass data, alpha and k*/
using pred_func_k_yes_t = std::function<KO_Traits::StoringVector(KO_Traits::StoringMatrix,double,int,int)>;
/*!Type for the function that predicts the validation set: k not imposed: pass data, alpha and threshold_ppc*/
using pred_func_k_no_t  = std::function<KO_Traits::StoringVector(KO_Traits::StoringMatrix,double,double,int)>;

/*!
* Type for the prediction function on validation set: depending on
* @param k_imp: enumerator K_IMP: if k is imposed or not
*/
template <K_IMP k_imp>
using pred_func_t = typename std::conditional<k_imp, pred_func_k_yes_t,pred_func_k_no_t>::type;


/*!
* @class CV_base
* @brief Template class for performing cross-validation.
* @tparam D type of the derived class (for static polymorphism thorugh CRTP):
*         - 'CV_alpha'
*         - 'CV_k'
*         - 'CV_alpha_k'
* @tparam cv_strat strategy for splitting training/validation sets
* @tparam err_eval how to evaluate the loss between prediction on validation set and validation set
* @tparam k_imp if k is imposed or has to be found through explanatory power criterion
* @tparam valid_err_ret if validation error are stored
* @details It is the base class. Polymorphism is known at compile time thanks to Curiously Recursive Template Pattern (CRTP) 
*/
template< class D, CV_STRAT cv_strat, CV_ERR_EVAL err_eval, K_IMP k_imp, VALID_ERR_RET valid_err_ret >
class CV_base
{
  
private:

  /*!Matrix storing fts*/
  KO_Traits::StoringMatrix m_Data;
  
  /*!Strategy for splitting training/validation set*/ 
  cv_strategy<cv_strat> m_strategy;

  /*!
  * @brief Evaluation of the loss between prediction on validation set and validation set
  * @param pred prediction on validation set
  * @param valid validation set
  * @param number_threads number of threads for multi-threading
  * @return the error between prediction on validation set and validation set
  */
  double err_valid_set_eval(const KO_Traits::StoringVector &pred, const KO_Traits::StoringVector &valid, int number_threads, ERR_EVAL_T<CV_ERR_EVAL::MSE>) const;
  
  /*!Number of threads for OMP*/
  int m_number_threads;
  
public:
  
  /*!
  * @brief Constructor for the class
  * @param Data fts matrix
  * @param strategy strategy for splitting training/validation sets
  * @param number_threads number of threads for OMP
  * @details Universal constructor: move semantic used to optimazing handling big size objects
  */
  template<typename STOR_OBJ,typename STRATEGY>
  CV_base(STOR_OBJ&& Data, STRATEGY && strategy, int number_threads)
    : m_Data{std::forward<STOR_OBJ>(Data)}, m_strategy{std::forward<STRATEGY>(strategy)}, m_number_threads(number_threads)  {}
  
  /*!
  * @brief Getter for the data matrix
  * @return the private m_Data
  */
  inline KO_Traits::StoringMatrix Data() const {return m_Data;}
  
  /*!
  * @brief Getter for the training/validation set splitting
  * @return the private m_strategy
  */
  inline cv_strategy<cv_strat> strategy() const {return m_strategy;}
  
  /*!
  * @brief Getter for the number of threads for OMP
  * @return the private m_number_threads
  */
  inline int number_threads() const {return m_number_threads;}
  
  /*!
  * @brief Evaluation of the loss between prediction on validation set and validation set: estimate of L2 norm loss. Tag-dispacther.
  * @param pred prediction on validation set
  * @param valid validation set
  * @param number_threads number of threads for multi-threading
  * @return the error between prediction on validation set and validation set
  */
  double err_valid_set_eval(const KO_Traits::StoringVector &pred, const KO_Traits::StoringVector &valid, int number_threads) const { return err_valid_set_eval(pred,valid,number_threads,ERR_EVAL_T<err_eval>{});};
  
};


#include "CV_err_valid_eval.hpp"

#endif  //CV_CRTP_PPC_HPP