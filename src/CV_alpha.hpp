#ifndef CV_CRTP_ALPHA_PPC_HPP
#define CV_CRTP_ALPHA_PPC_HPP

#include "CV.hpp"

#include <iostream>



template< CV_STRAT cv_strat, CV_ERR_EVAL err_eval, K_IMP k_imp, VALID_ERR_RET valid_err_ret >  
class CV_alpha : public CV_base< CV_alpha<cv_strat,err_eval,k_imp,valid_err_ret>, cv_strat, err_eval, k_imp, valid_err_ret >
{
private:
  std::vector<double> m_params;         //alphas
  valid_err_cv_1_t m_valid_errors;      //validation error for each alpha
  double m_best_valid_error;            //best validation error
  double m_param_best;                  //best alpha -> optimal parameter
  //parameters needed for KO
  int m_k = 0;
  double m_threshold_ppc = 0.0;
  pred_func_t<k_imp> m_pred_f;             //function to make predictions
  
public:
  
  //constructor if k is known
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
  
  //costructor if k is not known
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
  
  
  //Getters
  /*!
   * Getter for m_param_best
   */
  inline double param_best() const {return m_param_best;};
  
  /*!
   * Getter for m_valid_errors
   */
  inline std::vector<double> valid_errors() const {return m_valid_errors;};
  
  /*!
   * Getter for m_best_valid_error
   */
  inline double best_valid_error() const {return m_best_valid_error;};
  
  
   //error in a single CV iteration (specific training and validation set)
   inline 
   double 
   error_single_cv_iter(const double &param, const KO_Traits::StoringMatrix &training_set, const KO_Traits::StoringMatrix &validation_set)
   const
   {  
      //train the model to make prediction
      //auto pred = m_pred_f(training_set,m_threshold_ppc,param,m_k);
      if constexpr( k_imp == YES)
      {
        auto pred = m_pred_f(training_set,param,m_k,this->number_threads());
        return this->err_valid_set_eval(pred,validation_set,this->number_threads());
      }
      else
      {
        auto pred = m_pred_f(training_set,param,m_threshold_ppc,this->number_threads());
        return this->err_valid_set_eval(pred,validation_set,this->number_threads());
      }
      
   }
   
   
   //validation error for a single parameter
   inline 
   double                                 //vector di pair di vector
   error_single_param(const double &param, const cv_strategy_t &strat, const std::size_t &number_cv_iter) 
   const
   { 
     double err = 0.0;
     
#ifdef _OPENMP
#pragma omp parallel for shared(param,strat) num_threads(this->number_threads()) reduction(+:err)
     for (int i = 0; i < number_cv_iter; ++i)
     {
       auto train_valid_set = this->strategy().train_validation_set(this->Data(),strat[i]); 
       err += this->error_single_cv_iter(param,train_valid_set.first,train_valid_set.second);
     }
#else     
     err = std::transform_reduce(strat.cbegin(),
                                 strat.cend(),
                                 0.0,
                                 std::plus{},
                                 [this,&param](auto el){ auto train_valid_set = this->strategy().train_validation_set(this->Data(),el); return this->error_single_cv_iter(param,train_valid_set.first,train_valid_set.second);});
#endif
  
     return err/(static_cast<double>(number_cv_iter));
   }
   
  
  //how to select the best alpha
  inline 
  void 
  best_param_search() 
  { 
    int tot_params = m_params.size();
    
    //preparing the container for the errors: resize to use transform
    m_valid_errors.resize(tot_params);
    
#ifdef _OPENMP
#pragma omp parallel for shared(m_params,tot_params) num_threads(this->number_threads())
    for(std::size_t i = 0; i < tot_params; ++i)
    {
      m_valid_errors[i] = this->error_single_param(m_params[i],this->strategy().strategy(),this->strategy().strategy().size());
    }
#else
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