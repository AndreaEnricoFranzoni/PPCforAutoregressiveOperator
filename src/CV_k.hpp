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

#ifndef CV_CRTP_K_PPC_HPP
#define CV_CRTP_K_PPC_HPP

#include "CV.hpp"


/*!
* @file CV_k.hpp
* @brief Derived-from-CV_base class for performing cross-validation on the number of retained PPCs
* @author Andrea Enrico Franzoni
*/



/*!
* @class CV_k
* @brief Template class for performing cross-validation on the number of retained PPCs: derived-from-CV_base thorugh CRTP ('CV_k' its template D parameter).
* @tparam cv_strat strategy for splitting training/validation sets
* @tparam err_eval how to evaluate the loss between prediction on validation set and validation set
* @tparam k_imp if k is imposed or has to be found through explanatory power criterion (necessary k_imp=K_IMP::YES since k is imposed at each cv iteration)
* @tparam valid_err_ret if validation error are stored
* @details Derived class. Polymorphism is known at compile time thanks to Curiously Recursive Template Pattern (CRTP) 
*/
template< CV_STRAT cv_strat, CV_ERR_EVAL err_eval, K_IMP k_imp, VALID_ERR_RET valid_err_ret >  
class CV_k : public CV_base< CV_k<cv_strat,err_eval,k_imp,valid_err_ret>, cv_strat, err_eval, k_imp, valid_err_ret >
{
private:

  /*!Input space for number of PPCs*/
  std::vector<int> m_params;
  /*!Validation error for each element of the number of PPCs input space*/             
  valid_err_cv_1_t m_valid_errors;
  /*!Best validation error*/     
  double m_best_valid_error;
  /*!Optimal number of retained PPCs*/            
  int m_param_best;
  /*!Tolerance: the cv continues only if between two parameters, that are checked in increasing order, 
  * the absolute difference between two validation errors is bigger than tolerance*trace(covariance). 
  * If not, stops and look for k only between the tested ones 
  */                     
  double m_toll;                        
  /*!Regularization parameter, needed for predicting validation set*/
  double m_alpha;
  /*!Function to predict validation set*/
  pred_func_t<K_IMP::YES> m_pred_f;           
  

public:
  
  /*!
  * @brief Constructor for derived class: constructs firstly CV_base<CV_k,...>
  * @param Data fts data matrix
  * @param strategy splitting training/validation strategy
  * @param params input space for number of retained PPCs
  * @param toll tolerance between consecutive validation errors for looking for element with bigger value in the input space
  * @param alpha regularization parameter
  * @param pred_f function to make validation set prediction (overloading with k imposed)
  * @param number_threads number of threads for OMP
  * @details Universal constructor: move semantic used to optimazing handling big size objects
  */
  template<typename STOR_OBJ,typename STRATEGY>
  CV_k(STOR_OBJ&& Data,
       STRATEGY && strategy,
       const std::vector<int> &params,
       double toll,
       double alpha,
       const pred_func_t<K_IMP::YES> & pred_f,
       int number_threads)
    : CV_base<CV_k,cv_strat,err_eval,k_imp,valid_err_ret>(std::move(Data),std::move(strategy),number_threads), 
      m_params(params), 
      m_toll(toll),
      m_alpha(alpha),
      m_pred_f(pred_f)
      {}
  
  
  /*!
  * @brief Getter for the best number of retained PPCs
  * @return the private m_param_best
  */
  inline int param_best() const {return m_param_best;};
  
  /*!
  * @brief Getter for validation errors
  * @return the private m_valid_errors
  */
  inline std::vector<double> valid_errors() const {return m_valid_errors;};
  
  /*!
  * @brief Getter for the best validation error
  * @return the private m_best_valid_error
  */
  inline double best_valid_error() const {return m_best_valid_error;};
  
  
   /*!
   * @brief Error for a single cross-validation iteration (parameter given), with fixed training and validation sets
   * @param param parameter that is being evaluated through cross-validation
   * @param training_set training set
   * @param validation_set validation set
   * @return the error, according to 'err_eval' class template parameter, between prediction on validation set and validation set
   */
   inline 
   double 
   error_single_cv_iter(const int &param, const KO_Traits::StoringMatrix &training_set, const KO_Traits::StoringMatrix &validation_set)
   const
   {
      //training the model, making prediction and evaluating the error on the validaiton set
      auto pred = m_pred_f(training_set,m_alpha,param,this->number_threads());
      return this->err_valid_set_eval(pred,validation_set,this->number_threads());
   }
   
   
   /*!
   * @brief Validation error for a given parameter, as the mean of the errors on the various validation sets
   * @param param element of the input space for number of retained PPCs
   * @param strat strategy for splitting training and validation set
   * @param number_cv_iter total number of different splits
   * @return the average of the errors between prediction on validation set and validation set, for each split
   * @note eventual usage of 'pragma' directive for OMP
   */
   inline 
   double 
   error_single_param(const int &param, const cv_strategy_t &strat, const std::size_t &number_cv_iter) 
   const
   {
     double err = 0.0;
     
     //if OMP: going parallel
#ifdef _OPENMP
#pragma omp parallel for shared(param,strat) num_threads(this->number_threads()) reduction(+:err)
     for (int i = 0; i < number_cv_iter; ++i)
     {
       auto train_valid_set = this->strategy().train_validation_set(this->Data(),strat[i]); 
       err += this->error_single_cv_iter(param,train_valid_set.first,train_valid_set.second);
     }
#else
     // if not OMP: STL algorithms: transform_reduce to apply at each split the evaluation of the prediction and then summing them up
     err = std::transform_reduce(strat.cbegin(),
                                 strat.cend(),
                                 0.0,
                                 std::plus{},
                                 [this,&param](auto el){ auto train_valid_set = this->strategy().train_validation_set(this->Data(),el); return this->error_single_cv_iter(param,train_valid_set.first,train_valid_set.second);});
#endif     
     
     //returning the average
     return err/(static_cast<double>(number_cv_iter));
   }
  
  
  /*!
  * @brief Selecting the best number of retained PPCs, modifying it into the class
  */
  inline 
  void 
  best_param_search() 
  { 
    //preparing the container for the validation errors
    m_valid_errors.reserve(m_params.size());
    
    //evaluating validation error for each parameter
    double previous_error(static_cast<double>(0));
    
    //if adding another PPC does not improve too much the validation error: break
    for(const auto & el : m_params)
    {
      //evaluate the error for the parameter
      double curr_err = this->error_single_param(el,this->strategy().strategy(),this->strategy().strategy().size());
      m_valid_errors.emplace_back(curr_err);
      
      if(std::abs(curr_err - previous_error) < m_toll) {break;} else {previous_error = curr_err;}
    }
    
    //Shrinking
    m_valid_errors.shrink_to_fit();

    //best validation error
    auto min_err = (std::min_element(m_valid_errors.begin(),m_valid_errors.end()));
    m_best_valid_error = *min_err;
    
    //optimal param
    m_param_best = m_params[std::distance(m_valid_errors.begin(),min_err)];
    
    //if saving errors
    if constexpr (valid_err_ret == VALID_ERR_RET::NO_err)
    {
      m_valid_errors.clear();
    }
  }
};


#endif //CV_CRTP_K_PPC_HPP