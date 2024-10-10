#ifndef CV_CRTP_ALPHA_PPC_HPP
#define CV_CRTP_ALPHA_PPC_HPP

#include "CV_CRTP.hpp"


namespace CV_PPC
{

template<DEF_PARAMS_PPC::cv_strat_type cv_strat, DEF_PARAMS_PPC::cv_err_eval_type err_eval>  
class CV_KO_PPC_alpha_CRTP : public CV_KO_PPC_CRTP<CV_KO_PPC_alpha_CRTP<cv_strat,err_eval>, cv_strat, err_eval>
{
private:
  std::vector<double> m_params;     //alphas
  std::vector<double> m_valid_errors;     //error for each alpha
  double m_best_valid_error;
  double m_param_best;              //best alpha
  //parameters needed for KO
  double m_threshold_ppc;
  int m_k;
  pred_func_t m_pred_f;             //function to make predictions
  
public:
  
  //constructor
  CV_KO_PPC_alpha_CRTP(KO_Traits::StoringMatrix&& Data,
                       const std::vector<double> &params,
                       double threshold_ppc,
                       int k,
                       const pred_func_t & pred_f)
    : CV_KO_PPC_CRTP<CV_KO_PPC_alpha_CRTP, cv_strat, err_eval>(std::move(Data)), 
      m_params(params), 
      m_threshold_ppc(threshold_ppc),
      m_k(k),
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
   error_single_cv_iter(const double &param, const KO_Traits::StoringMatrix &train_set, const KO_Traits::StoringMatrix &validation_set)
   const
   { 
   //train the model to make prediction
   auto pred = m_pred_f(train_set,m_threshold_ppc,param,m_k);
   return this->err_valid_set_eval(pred,validation_set);
   }
   
   
   //validation error for a single parameter
   inline 
   double 
   error_single_param(const double &param, const cv_strategy_t &strat, const std::size_t &number_cv_iter) 
   const
   { 
   double err = std::transform_reduce(strat.cbegin(),
   strat.cend(),
   0.0,
   std::plus{},
   [this,&param](auto el){ auto train_valid_set = this->train_validation_set(el); return this->error_single_cv_iter(param,train_valid_set.first,train_valid_set.second);});
   
   return err/(static_cast<double>(number_cv_iter));
   }
   
  
  //how to select the best alpha
  inline 
  void 
  best_param_search() 
  { 
    //strategy for doing CV
    cv_strategy_t strategy_cv = this->train_validation_set_strategy(this->Data());
    std::size_t number_cv_iter = strategy_cv.size();   //how many cv iterations are done
    
    //preparing the container for the errors: resize to use transform
    m_valid_errors.resize(m_params.size());
    
    //evaluating validation error for each parameter
    std::transform(m_params.cbegin(),
                   m_params.cend(),
                   m_valid_errors.begin(),
                   [this,&strategy_cv,&number_cv_iter](double const &param_i){return this->error_single_param(param_i,strategy_cv,number_cv_iter);});
    
    //best validation error
    auto min_err = (std::min_element(m_valid_errors.begin(),m_valid_errors.end()));
    m_best_valid_error = *min_err;
    
    //optimal param
    m_param_best = m_params[std::distance(m_valid_errors.begin(),min_err)];
    
    //free memory
    strategy_cv.clear();
  }
  
};
                           
}   //end namespace CV_PPC

#endif //CV_CRTP_ALPHA_PPC_HPP