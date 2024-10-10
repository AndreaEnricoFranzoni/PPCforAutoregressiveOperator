#ifndef CV_CRTP_K_PPC_HPP
#define CV_CRTP_K_PPC_HPP

#include "CV_CRTP.hpp"


namespace CV_PPC
{

template<DEF_PARAMS_PPC::cv_strat_type cv_strat, DEF_PARAMS_PPC::cv_err_eval_type err_eval>  
class CV_KO_PPC_k_CRTP : public CV_KO_PPC_CRTP<CV_KO_PPC_k_CRTP<cv_strat,err_eval>, cv_strat, err_eval>
{
private:
  std::vector<int> m_params;        //k_s
  std::vector<double> m_valid_errors;     //error for each k
  double m_best_valid_error;        //best validation error
  int m_param_best;                 //best k
  double m_toll;                    //tolerance for stopping if extra components are meaningless
  //parameters needed for KO
  double m_threshold_ppc;
  double m_alpha;
  pred_func_t m_pred_f;             //function to make predictions
  
public:
  
  //constructor
  CV_KO_PPC_k_CRTP(KO_Traits::StoringMatrix&& Data,
                   const std::vector<int> &params,
                   double toll,
                   double threshold_ppc,
                   double alpha,
                   const pred_func_t & pred_f)
    : CV_KO_PPC_CRTP<CV_KO_PPC_k_CRTP, cv_strat, err_eval>(std::move(Data)), 
      m_params(params), 
      m_toll(toll),
      m_threshold_ppc(threshold_ppc),
      m_alpha(alpha),
      m_pred_f(pred_f)
      {}
  
  
  //Getters
  /*!
   * Getter for m_param_best
   */
  inline int param_best() const {return m_param_best;};
  
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
   error_single_cv_iter(const int &param, const KO_Traits::StoringMatrix &train_set, const KO_Traits::StoringMatrix &validation_set)
   const
   {
   //train the model to make prediction
   auto pred = m_pred_f(train_set,m_threshold_ppc,m_alpha,param);
   return this->err_valid_set_eval(pred,validation_set);
   }
   
   
   //validation error for a single parameter
   inline 
   double 
   error_single_param(const int &param, const cv_strategy_t &strat, const std::size_t &number_cv_iter) 
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
    auto strategy_cv = this->train_validation_set_strategy(this->Data());
    std::size_t number_cv_iter = strategy_cv.size();   //how many cv iterations are done
    
    //preparing the container for the validation errors
    m_valid_errors.reserve(m_params.size());
    
    //evaluating validation error for each parameter
    double previous_error(static_cast<double>(0));
    
    //if adding another PPC does not improve too much the validation error: break
    for(const auto & el : m_params)
    {
      //evaluate the error for the parameter
      double curr_err = this->error_single_param(el,strategy_cv,number_cv_iter);
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
    
    //free memory
    strategy_cv.clear();
  }
  
};

}   //end namespace CV_PPC

#endif //CV_CRTP_K_PPC_HPP



/*
 * //strategy for doing CV
 auto strategy_cv = this->train_validation_set_strategy(this->Data());
 std::size_t number_cv_iter = strategy_cv.size();   //how many cv iterations are done
 
 //preparing the container for the validation errors
 m_valid_errors.resize(m_params.size());
 
 //evaluating validation error for each parameter
 std::size_t counter_k = 0;
 double previous_error = static_cast<double>(0);
 
 //if adding another PPC does not improve too much the validation error: break
 for(const auto & el : m_params)
 {
 m_valid_errors[counter_k]=(this->error_single_param(el,strategy_cv,number_cv_iter));
 if(std::abs(m_valid_errors[counter_k] - previous_error) < m_toll)     //already found a point in which you do not improve anymore
 {
 //m_errors.resize(el);
 //break;
 for(std::size_t i = counter_k; i < m_params.size(); ++i){m_valid_errors[i]=(m_valid_errors[counter_k]);}
 }
 else
 {
 previous_error = m_valid_errors[counter_k];
 counter_k++;
 }
 }
 
 //best validation error
 auto min_err = (std::min_element(m_valid_errors.begin(),m_valid_errors.end()));
 m_best_valid_error = *min_err;
 
 //optimal param
 m_param_best = m_params[std::distance(m_valid_errors.begin(),min_err)];
 
 //free memory
 strategy_cv.clear();
 */