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

#ifndef CV_CRTP_ALPHA_K_PPC_HPP
#define CV_CRTP_ALPHA_K_PPC_HPP

#include "CV.hpp"
#include "CV_alpha_k.hpp"


/*!
* @file CV_alpha_k.hpp
* @brief Derived-from-CV_base class for performing cross-validation on both the regularization parameter and the number of retained PPCs
* @author Andrea Enrico Franzoni
*/



/*!
* @class CV_alpha_k
* @brief Template class for performing cross-validation on the regularization parameter and the number of PPCs: derived-from-CV_base thorugh CRTP ('CV_alpha_k' its template D parameter).
* @tparam cv_strat strategy for splitting training/validation sets
* @tparam err_eval how to evaluate the loss between prediction on validation set and validation set
* @tparam k_imp if k is imposed or has to be found through explanatory power criterion (necessary k_imp=K_IMP::YES since k is imposed at each cv iteration)
* @tparam valid_err_ret if validation error are stored
* @details Derived class. Polymorphism is known at compile time thanks to Curiously Recursive Template Pattern (CRTP) 
*/
template< CV_STRAT cv_strat, CV_ERR_EVAL err_eval, K_IMP k_imp, VALID_ERR_RET valid_err_ret >  
class CV_alpha_k : public CV_base< CV_alpha_k<cv_strat,err_eval,k_imp,valid_err_ret>, cv_strat, err_eval, k_imp, valid_err_ret >
{
private:

  /*!Input space for regularization parameter*/
  std::vector<double> m_alphas;
  /*!Input space for number of PPCs*/                   
  std::vector<int> m_k_s;
  /*!Optimal reguòlarization parameter*/                         
  double m_alpha_best;                            
  /*!Optimal number of retained PPCs*/ 
  int m_k_best;
  /*!Validation error for the optimal pair*/                                   
  double m_best_valid_error;                      
  /*!For each regularizatio parameter (key): its optimal number of retained PPCs (value)*/
  std::map<double,int> m_best_pairs;              
  /*!Errors for each pair regularization paramter-number of retained PPCs*/
  valid_err_cv_2_t m_valid_errors;                
  /*!Validation errors for each one of the regularization parameter with its best number of retained PPCs (regularization parameters ordered in increasing order)*/
  std::vector<double> m_valid_errors_best_pairs;  
  /*!Tolerance: the cv continues only if between two parameters, that are checked in increasing order, 
  * the absolute difference between two validation errors is bigger than tolerance*trace(covariance). 
  * If not, stops and look for k only between the tested ones 
  */ 
  double m_toll;
  /*!Function to predict validation set*/
  pred_func_t<K_IMP::YES> m_pred_f;               
  

public:

  /*!
  * @brief Constructor for derived class: constructs firstly CV_base<CV_alpha_k,...>
  * @param Data fts data matrix
  * @param strategy splitting training/validation strategy
  * @param alphas input space for the regularization parameter
  * @param k_s input space for number of retained PPCs
  * @param toll tolerance between consecutive validation errors for looking for element with bigger value in the input space
  * @param pred_f function to make validation set prediction (overloading with k imposed)
  * @param number_threads number of threads for OMP
  * @details Universal constructor: move semantic used to optimazing handling big size objects
  */
  template<typename STOR_OBJ,typename STRATEGY>
  CV_alpha_k(STOR_OBJ&& Data,
             STRATEGY && strategy,
             const std::vector<double> &alphas,
             const std::vector<int> &k_s,
             double toll,
             const pred_func_t<K_IMP::YES> & pred_f,
             int number_threads)
    : CV_base<CV_alpha_k,cv_strat,err_eval,k_imp,valid_err_ret>(std::move(Data),std::move(strategy),number_threads), 
      m_alphas(alphas),
      m_k_s(k_s),
      m_toll(toll),
      m_pred_f(pred_f)
      {}
  
  
  /*!
  * @brief Getter for the best regularization parameter
  * @return the private m_alpha_best
  */
  inline double alpha_best() const {return m_alpha_best;};
  
  /*!
  * @brief Getter for the best number of retained PPCs
  * @return the private m_k_best
  */
  inline double k_best() const {return m_k_best;};
  
  /*!
  * @brief Getter for validation errors for the best pairs
  * @return the private m_valid_errors_best_pairs
  */
  inline std::vector<double> valid_errors_best_pairs() const {return m_valid_errors_best_pairs;};
  
  /*!
  * @brief Getter for validation errors
  * @return the private m_valid_errors
  */
  inline std::vector<std::vector<double>> valid_errors() const {return m_valid_errors;};
  
  /*!
  * @brief Getter for the best validation error
  * @return the private m_best_valid_error
  */
  inline double best_valid_error() const {return m_best_valid_error;};
  
  
  /*!
  * @brief Retaining the best pair regularization parameter-number of retained PPCs
  * @details for each element of the input space for regularization parameters, a cross-validation on the number of retained PPCs is performed.
  *          Consequently, the best pair is looked for within this ones.
  * @note eventual usage of 'pragma' directive for OMP
  */
  inline 
  void 
  best_param_search()
  {
    std::size_t tot_alphas = m_alphas.size();
    
    //if OMP: going parallel
#ifdef _OPENMP
    //preparing the containers for the errors
    if constexpr(valid_err_ret == VALID_ERR_RET::YES_err)
    {
      m_valid_errors.resize(tot_alphas);  //for each alpha: valid error for each k
    }
    
    m_valid_errors_best_pairs.resize(tot_alphas); //for each alpha: valid error only for the best k

#pragma omp parallel for shared(m_alphas,m_k_s,m_toll,m_pred_f,tot_alphas) num_threads(this->number_threads()) schedule(guided)
    for(std::size_t i = 0; i < tot_alphas; ++i)
    {
      //alpha fixed: doing CV on k
      CV_k<cv_strat,err_eval,k_imp,valid_err_ret> cv(std::move(this->Data()),std::move(this->strategy()),m_k_s,m_toll,m_alphas[i],m_pred_f,this->number_threads());
      cv.best_param_search();
      
      //best k given the alpha
      m_best_pairs.insert(std::make_pair(m_alphas[i],cv.param_best()));
      //saving the validation error for the best pair
      m_valid_errors_best_pairs[i]=cv.best_valid_error();
      
      if constexpr(valid_err_ret == VALID_ERR_RET::YES_err)
      {
        //saving the validation error for each pair (further inspection)
        m_valid_errors[i]=cv.valid_errors();
      }
    }
#else
    //preparing the containers for the errors
    if constexpr(valid_err_ret == VALID_ERR_RET::YES_err)
    {
      m_valid_errors.reserve(tot_alphas);  //for each alpha: valid error for each k
    }
    
    m_valid_errors_best_pairs.reserve(tot_alphas); //for each alpha: valid error only for the best k
    
    for(const auto & alpha : m_alphas)
    {
      //alpha fixed: doing CV on k
      CV_k<cv_strat,err_eval,k_imp,valid_err_ret> cv(std::move(this->Data()),std::move(this->strategy()),m_k_s,m_toll,alpha,m_pred_f,this->number_threads());
      cv.best_param_search();
      
      //best k given the alpha
      m_best_pairs.insert(std::make_pair(alpha,cv.param_best()));
      //saving the validation error for the best pair
      m_valid_errors_best_pairs.emplace_back(cv.best_valid_error());
      
      if constexpr(valid_err_ret == VALID_ERR_RET::YES_err)
      {
        //saving the validation error for each pair (further inspection)
        m_valid_errors.emplace_back(cv.valid_errors());
      }
    }
#endif
    
    //best validation error
    auto min_err = std::min_element(m_valid_errors_best_pairs.begin(),m_valid_errors_best_pairs.end());
    m_best_valid_error = *min_err;
    
    //best alpha
    m_alpha_best = m_alphas[std::distance(m_valid_errors_best_pairs.begin(),min_err)];
    
    //best k
    m_k_best = m_best_pairs.find(m_alpha_best)->second;
  }
  
};
  

#endif  //CV_CRTP_ALPHA_K_PPC_HPP