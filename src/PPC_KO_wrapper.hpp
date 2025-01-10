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

#ifndef PPC_KO_WRAPPER_HPP
#define PPC_KO_WRAPPER_HPP

#include <vector>
#include <tuple>
#include "utility"


#include "traits_ko.hpp"
#include "PPC_KO_include.hpp"


/*!
* @file PPC_KO_wrapper.hpp
* @brief Hierarchy of classes, for different algorithm versions. Correct child class is constructed by the KO_factory, providing the solver. Wrap the real class for computations
* @author Andrea Enrico Franzoni
*/



/*!
* @class PPC_KO_wrapper
* @brief Base virtual class for wrapping class that performs PPCKO computations. Which child class is constructed is selected run-time through virtual polymorphism
* @tparam solver if algorithm solved inverting the regularized covariance or avoiding it through gep (not possible if retaining the number of PPCs with explanatory power criterion)
* @tparam k_imp if k is imposed or has to be found through explanatory power criterion
* @tparam valid_err_ret if validation error are stored
* @tparam cv_strat strategy for splitting training/validation sets
* @tparam err_eval how to evaluate the loss between prediction on validation set and validation set
* @details It is a base class. Polymorphism is known at run-time through virtual polymorphism
*/
template< SOLVER solver, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval > 
class PPC_KO_wrapper
{
private:

  /*!Data matrix storing fts*/
  KO_Traits::StoringMatrix m_data;
  /*!PPCKO algorithm results*/      
  results_t<valid_err_ret> m_results;   
  /*!Number of threads for OMP*/
  int m_number_threads;                
  

public:

  /*!
  * @brief Constructor
  * @param data matrix storing fts
  * @param number_threads number of threads for OMP
  * @details Universal constructor: move semantic used to optimazing handling big size objects
  */
  template<typename STOR_OBJ>
  PPC_KO_wrapper(STOR_OBJ&& data, int number_threads)
    :  m_data{std::forward<STOR_OBJ>(data)}, m_number_threads(number_threads)   {}
  
  /*!
  * @brief Virtual destructor
  */
  virtual ~PPC_KO_wrapper() = default;
  
  /*!
  * @brief Virtual method to call the correct PPCKO version at runtime
  */
  virtual void call_ko() = 0; 
  
  /*!
  * @brief Getter for the data matrix
  * @return the private m_data
  */
  inline KO_Traits::StoringMatrix data() const {return m_data;};
  
  /*!
  * @brief Getter for the results
  * @return the private m_results
  */
  inline results_t<valid_err_ret> results() const {return m_results;};
  
  /*!
  * @brief Getter for the number of threads for OMP
  * @return the private m_number_threads
  */
  inline int number_threads() const {return m_number_threads;};
  
  /*!
  * @brief Setter for the results
  * @return the private m_results (not-const)
  */
  inline results_t<valid_err_ret> & results() {return m_results;};
};



/*!
* @class PPC_KO_wrapper_no_cv
* @brief Derived-from-PPC_KO_wrapper class for wrapping class that performs PPCKO computations without cross-validation
* @tparam solver if algorithm solved inverting the regularized covariance or avoiding it through gep (not possible if retaining the number of PPCs with explanatory power criterion)
* @tparam k_imp if k is imposed or has to be found through explanatory power criterion
* @tparam valid_err_ret if validation error are stored
* @tparam cv_strat strategy for splitting training/validation sets
* @tparam err_eval how to evaluate the loss between prediction on validation set and validation set
* @details It is a derived class. Polymorphism is known at run-time through virtual polymorphism
*/
template< SOLVER solver, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval > 
class PPC_KO_wrapper_no_cv  : public PPC_KO_wrapper<solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>
{
private:

  /*!Regularization parameter*/
  double m_alpha;
  /*!Number of retained PPCs*/                     
  int m_k;                              
  /*!Requested explanatory power from the PPCs*/
  double m_threshold_ppc;               


public:

  /*!
  * @brief Constructor if k is imposed (k_imp = K_IMP::YES)
  * @param data matrix storing fts
  * @param alpha regularization parameter
  * @param k number of retained PPCs
  * @param number_threads number of threads for OMP
  * @details Universal constructor: move semantic used to optimazing handling big size objects
  */
  template<typename STOR_OBJ>
  PPC_KO_wrapper_no_cv(STOR_OBJ&& data, double alpha, int k, int number_threads)
    : PPC_KO_wrapper<solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>(std::move(data),number_threads), m_alpha(alpha), m_k(k) {}
  
  /*!
  * @brief Constructor if k is retained through explanatory power criterion (k_imp = K_IMP::NO)
  * @param data matrix storing fts
  * @param alpha regularization parameter
  * @param threshold_ppc requested explanatory power from the PPCs
  * @param number_threads number of threads for OMP
  * @details Universal constructor: move semantic used to optimazing handling big size objects
  */
  template<typename STOR_OBJ>
  PPC_KO_wrapper_no_cv(STOR_OBJ&& data, double alpha, double threshold_ppc, int number_threads)
    : PPC_KO_wrapper<solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>(std::move(data),number_threads), m_alpha(alpha), m_threshold_ppc(threshold_ppc)  {}
  
  /*!
  * @brief Override for calling the without cross-validation PPCKO version at runtime
  */
  void call_ko() override;
};



/*!
* @class PPC_KO_wrapper_cv_alpha
* @brief Derived-from-PPC_KO_wrapper class for wrapping class that performs PPCKO computations with cross-validation on regularization parameter
* @tparam solver if algorithm solved inverting the regularized covariance or avoiding it through gep (not possible if retaining the number of PPCs with explanatory power criterion)
* @tparam k_imp if k is imposed or has to be found through explanatory power criterion
* @tparam valid_err_ret if validation error are stored
* @tparam cv_strat strategy for splitting training/validation sets
* @tparam err_eval how to evaluate the loss between prediction on validation set and validation set
* @details It is a derived class. Polymorphism is known at run-time through virtual polymorphism
*/
template< SOLVER solver, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval > 
class PPC_KO_wrapper_cv_alpha  : public PPC_KO_wrapper<solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>
{
private:

  /*!Regularization parameter input space*/
  std::vector<double> m_alphas;   
  /*!Number of retained PPCs*/         
  int m_k;
  /*!Requested explanatory power from the PPCs*/
  double m_threshold_ppc;
  /*!Minimum size (number of time instants) of the training set*/
  int m_min_size_ts;
  /*!Maximum size (number of time instants) of the training set*/ 
  int m_max_size_ts;


public:

  /*!
  * @brief Constructor if k is imposed (k_imp = K_IMP::YES)
  * @param data matrix storing fts
  * @param alphas regularization parameter input space
  * @param k number of retained PPCs
  * @param min_size_ts minimum size (number of time instants) of the training set
  * @param max_size_ts maximum size (number of time instants) of the training set
  * @param number_threads number of threads for OMP
  * @details Universal constructor: move semantic used to optimazing handling big size objects
  */
  template<typename STOR_OBJ>
  PPC_KO_wrapper_cv_alpha(STOR_OBJ&& data, const std::vector<double> & alphas, int k, int min_size_ts, int max_size_ts, int number_threads)
    : PPC_KO_wrapper<solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>(std::move(data),number_threads), m_alphas(alphas), m_k(k), m_min_size_ts(min_size_ts), m_max_size_ts(max_size_ts) {}
  
  /*!
  * @brief Constructor if k is retained through explanatory power criterion (k_imp = K_IMP::NO)
  * @param data matrix storing fts
  * @param alphas regularization parameter input space
  * @param threshold_ppc requested explanatory power from the PPCs
  * @param min_size_ts minimum size (number of time instants) of the training set
  * @param max_size_ts maximum size (number of time instants) of the training set
  * @param number_threads number of threads for OMP
  * @details Universal constructor: move semantic used to optimazing handling big size objects
  */
  template<typename STOR_OBJ>
  PPC_KO_wrapper_cv_alpha(STOR_OBJ&& data, const std::vector<double> & alphas, double threshold_ppc, int min_size_ts, int max_size_ts, int number_threads)
    : PPC_KO_wrapper<solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>(std::move(data),number_threads), m_alphas(alphas), m_threshold_ppc(threshold_ppc), m_min_size_ts(min_size_ts), m_max_size_ts(max_size_ts) {}
  
  /*!
  * @brief Override for calling the regularization parameter cross-validation PPCKO version at runtime
  */
  void call_ko() override;
};



/*!
* @class PPC_KO_wrapper_cv_k
* @brief Derived-from-PPC_KO_wrapper class for wrapping class that performs PPCKO computations with cross-validation on the number of retained PPCs
* @tparam solver if algorithm solved inverting the regularized covariance or avoiding it through gep (not possible if retaining the number of PPCs with explanatory power criterion)
* @tparam k_imp if k is imposed or has to be found through explanatory power criterion
* @tparam valid_err_ret if validation error are stored
* @tparam cv_strat strategy for splitting training/validation sets
* @tparam err_eval how to evaluate the loss between prediction on validation set and validation set
* @details It is a derived class. Polymorphism is known at run-time through virtual polymorphism
*/
template< SOLVER solver, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval > 
class PPC_KO_wrapper_cv_k  : public PPC_KO_wrapper<solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>
{
private:

  /*!Regularization parameter*/
  double m_alpha;
  /*!Number of retained PPCs input space*/
  std::vector<int> m_k_s;
  /*!Tolerance: the cv continues only if between two parameters, that are checked in increasing order, 
  * the absolute difference between two validation errors is bigger than tolerance*trace(covariance). 
  * If not, stops and look for k only between the tested ones 
  */   
  double m_toll;
  /*!Minimum size (number of time instants) of the training set*/
  int m_min_size_ts; 
  /*!Maximum size (number of time instants) of the training set*/
  int m_max_size_ts;
  

public:

  /*!
  * @brief Constructor 
  * @param data matrix storing fts
  * @param alpha regularization parameter
  * @param k_s number of retained PPCs input space
  * @param toll the cv on the number of retained PPCs continues only if between two parameters, that are checked in increasing order, the absolute difference between two validation errors is bigger than tolerance*trace(covariance). If not, stops and look for k only between the tested ones 
  * @param min_size_ts minimum size (number of time instants) of the training set
  * @param max_size_ts maximum size (number of time instants) of the training set
  * @param number_threads number of threads for OMP
  * @details Universal constructor: move semantic used to optimazing handling big size objects
  */
  template<typename STOR_OBJ>
  PPC_KO_wrapper_cv_k(STOR_OBJ&& data, double alpha, const std::vector<int> & k_s, double toll, int min_size_ts, int max_size_ts, int number_threads)
    : PPC_KO_wrapper<solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>(std::move(data),number_threads), m_alpha(alpha), m_k_s(k_s), m_toll(toll), m_min_size_ts(min_size_ts), m_max_size_ts(max_size_ts) {}
  
  /*!
  * @brief Override for calling the number of retained PPCs cross-validation PPCKO version at runtime
  */
  void call_ko() override;
};



/*!
* @class PPC_KO_wrapper_cv_alpha_k
* @brief Derived-from-PPC_KO_wrapper class for wrapping class that performs PPCKO computations with cross-validation on both the regularization parameter and the number of retained PPCs
* @tparam solver if algorithm solved inverting the regularized covariance or avoiding it through gep (not possible if retaining the number of PPCs with explanatory power criterion)
* @tparam k_imp if k is imposed or has to be found through explanatory power criterion
* @tparam valid_err_ret if validation error are stored
* @tparam cv_strat strategy for splitting training/validation sets
* @tparam err_eval how to evaluate the loss between prediction on validation set and validation set
* @details It is a derived class. Polymorphism is known at run-time through virtual polymorphism
*/
template< SOLVER solver, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval > 
class PPC_KO_wrapper_cv_alpha_k  : public PPC_KO_wrapper<solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>
{
private:

  /*!Regularization parameter input space*/
  std::vector<double> m_alphas;
  /*!Number of retained PPCs input space*/
  std::vector<int> m_k_s;
  /*!Tolerance: the cv continues only if between two parameters, that are checked in increasing order, 
  * the absolute difference between two validation errors is bigger than tolerance*trace(covariance). 
  * If not, stops and look for k only between the tested ones 
  */  
  double m_toll;
  /*!Minimum size (number of time instants) of the training set*/
  int m_min_size_ts; 
  /*!Maximum size (number of time instants) of the training set*/
  int m_max_size_ts;
  

public:

  /*!
  * @brief Constructor 
  * @param data matrix storing fts
  * @param alphas regularization parameter input space
  * @param k_s number of retained PPCs input space
  * @param toll the cv on the number of retained PPCs continues only if between two parameters, that are checked in increasing order, the absolute difference between two validation errors is bigger than tolerance*trace(covariance). If not, stops and look for k only between the tested ones 
  * @param min_size_ts minimum size (number of time instants) of the training set
  * @param max_size_ts maximum size (number of time instants) of the training set
  * @param number_threads number of threads for OMP
  * @details Universal constructor: move semantic used to optimazing handling big size objects
  */
  template<typename STOR_OBJ>
  PPC_KO_wrapper_cv_alpha_k(STOR_OBJ&& data, const std::vector<double> &alphas, const std::vector<int> &k_s, double toll, int min_size_ts, int max_size_ts, int number_threads)
    : PPC_KO_wrapper<solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>(std::move(data),number_threads), m_alphas(alphas), m_k_s(k_s), m_toll(toll), m_min_size_ts(min_size_ts), m_max_size_ts(max_size_ts) {}
  
  /*!
  * @brief Override for calling the regularization parameter-number of retained PPCs cross-validation PPCKO version at runtime
  */
  void call_ko() override;
};


#include "PPC_KO_wrapper_imp.hpp"

#endif  //PPC_KO_WRAPPER_HPP