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

#ifndef CV_CRTP_ALPHA_PPC_HPP
#define CV_CRTP_ALPHA_PPC_HPP

#include "CV.hpp"



/*!
* @file CV_alpha.hpp
* @brief Derived-from-CV_base class for performing cross-validation on the regularization parameter alpha
* @author Andrea Enrico Franzoni
*/


/*!
* @class CV_alpha
* @brief Template class for performing cross-validation on the regularization parameter: derived-from-CV_base thorugh CRTP ('CV_alpha' its template D parameter).
* @tparam cv_strat strategy for splitting training/validation sets
* @tparam err_eval how to evaluate the loss between prediction on validation set and validation set
* @tparam k_imp if k is imposed or has to be found through explanatory power criterion
* @tparam valid_err_ret if validation error are stored
* @details Derived class. Polymorphism is known at compile time thanks to Curiously Recursive Template Pattern (CRTP) 
*/
template< CV_STRAT cv_strat, CV_ERR_EVAL err_eval, K_IMP k_imp, VALID_ERR_RET valid_err_ret >  
class CV_alpha : public CV_base< CV_alpha<cv_strat,err_eval,k_imp,valid_err_ret>, cv_strat, err_eval, k_imp, valid_err_ret >
{
private:

  /*!Input space for regularization parameter*/
  std::vector<double> m_params; 
  /*!Validation error for each element of the regularization parameter input space*/        
  valid_err_cv_1_t m_valid_errors;      
  /*!Best validation error*/
  double m_best_valid_error;  
  /*!Optimal regularization parameter*/          
  double m_param_best;                  
  /*!Number of PPCs, needed for predicting validation set*/
  int m_k = 0;
  /*!Threshold for explanatory power of PPCs, needed for predicting validation set*/
  double m_threshold_ppc = 0.0;
  /*!Function to predict validation set*/
  pred_func_t<k_imp> m_pred_f;             
  

public:
  
  /*!
  * @brief Constructor if number of PPCs k is known (k_imp=K_IMP::YES), for derived class: constructs firstly CV_base<CV_alpha,...>
  * @param Data fts data matrix
  * @param strategy splitting training/validation strategy
  * @param params input space for regularization parameter
  * @param k number of retained PPCs
  * @param pred_f function to make validation set prediction (overloading with k imposed)
  * @param number_threads number of threads for OMP
  * @details Universal constructor: move semantic used to optimazing handling big size objects
  */
  template<typename STOR_OBJ,typename STRATEGY>
  CV_alpha(STOR_OBJ&& Data,
           STRATEGY && strategy,
           const std::vector<double> &params,
           int k,
           const pred_func_t<k_imp> & pred_f,
           int number_threads)
    : CV_base<CV_alpha,cv_strat,err_eval,k_imp,valid_err_ret>(std::move(Data),std::move(strategy),number_threads), 
      m_params(params), 
      m_k(k),
      m_pred_f(pred_f)
      {}
  
  
  /*!
  * @brief Constructor if number of PPCs k is not known (k_imp=K_IMP::NO), but has to be retained using explanatory power criterion, for derived class: constructs firstly CV_base<CV_alpha,...>
  * @param Data fts data matrix
  * @param strategy splitting training/validation strategy
  * @param params input space for regularization parameter
  * @param threshold_ppc requested explanatory power for PPCs
  * @param pred_f function to make validation set prediction (overloading with k not imposed)
  * @param number_threads number of threads for OMP
  * @details Universal constructor: move semantic used to optimazing handling big size objects
  */
  template<typename STOR_OBJ,typename STRATEGY>
  CV_alpha(STOR_OBJ&& Data,
           STRATEGY && strategy,
           const std::vector<double> &params,
           double threshold_ppc,
           const pred_func_t<k_imp> & pred_f,
           int number_threads)
    : CV_base<CV_alpha,cv_strat,err_eval,k_imp,valid_err_ret>(std::move(Data),std::move(strategy),number_threads), 
      m_params(params), 
      m_threshold_ppc(threshold_ppc),
      m_pred_f(pred_f)
      {}
  
  
  /*!
  * @brief Getter for the best regularization parameter
  * @return the private m_param_best
  */
  inline double param_best() const {return m_param_best;};
  
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
   error_single_cv_iter(const double &param, const KO_Traits::StoringMatrix &training_set, const KO_Traits::StoringMatrix &validation_set)
   const
   {  
      //traiing the model, making the prediction and then evaluating it
      if constexpr( k_imp == YES)     //k is imposed
      {
        auto pred = m_pred_f(training_set,param,m_k,this->number_threads());
        return this->err_valid_set_eval(pred,validation_set,this->number_threads());
      }
      else                            //explanatory power for retained PPCs
      {
        auto pred = m_pred_f(training_set,param,m_threshold_ppc,this->number_threads());
        return this->err_valid_set_eval(pred,validation_set,this->number_threads());
      } 
   }
   
   
   /*!
   * @brief Validation error for a given parameter, as the mean of the errors on the various validation sets
   * @param param element of the input space for regularization parameter
   * @param strat strategy for splitting training and validation set
   * @param number_cv_iter total number of different splits
   * @return the average of the errors between prediction on validation set and validation set, for each split
   * @note eventual usage of 'pragma' directive for OMP
   */
   inline 
   double                                 
   error_single_param(const double &param, const cv_strategy_t &strat, const std::size_t &number_cv_iter) 
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
  * @brief Selecting the best regularization parameter, modifying it into the class
  * @note eventual usage of 'pragma' directive for OMP
  */
  inline 
  void 
  best_param_search() 
  { 
    int tot_params = m_params.size();
    
    //preparing the container for the errors: resize to use transform
    m_valid_errors.resize(tot_params);

    //if OMP: going parallel
#ifdef _OPENMP
#pragma omp parallel for shared(m_params,tot_params) num_threads(this->number_threads())
    for(std::size_t i = 0; i < tot_params; ++i)
    {
      m_valid_errors[i] = this->error_single_param(m_params[i],this->strategy().strategy(),this->strategy().strategy().size());
    }
#else
    // if not OMP: STL algorithms: transform to apply at each element the input space the function to retain their validation error
    std::transform(m_params.cbegin(),
                   m_params.cend(),
                   m_valid_errors.begin(),
                   [this](double const &param_i){return this->error_single_param(param_i,this->strategy().strategy(),this->strategy().strategy().size());});
#endif
    
    //best validation error
    auto min_err = (std::min_element(m_valid_errors.begin(),m_valid_errors.end()));
    m_best_valid_error = *min_err;
    
    //optimal param
    m_param_best = m_params[std::distance(m_valid_errors.begin(),min_err)];
    
    //free memory
    if constexpr(valid_err_ret == VALID_ERR_RET::NO_err)
    {
      m_valid_errors.clear();
    }
  }
  
};
                           

#endif //CV_CRTP_ALPHA_PPC_HPP