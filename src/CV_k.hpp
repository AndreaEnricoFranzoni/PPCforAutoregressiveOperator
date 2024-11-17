#ifndef CV_CRTP_K_PPC_HPP
#define CV_CRTP_K_PPC_HPP

#include "CV.hpp"


template< CV_STRAT cv_strat, CV_ERR_EVAL err_eval, K_IMP k_imp, VALID_ERR_RET valid_err_ret >  
class CV_k : public CV_base< CV_k<cv_strat,err_eval,k_imp,valid_err_ret>, cv_strat, err_eval, k_imp, valid_err_ret >
{
private:
  std::vector<int> m_params;            //k_s
  valid_err_cv_1_t m_valid_errors;      //error for each k
  double m_best_valid_error;            //best validation error
  int m_param_best;                     //best k
  double m_toll;                        //tolerance for stopping if extra components are meaningless
  //parameters needed for KO
  double m_alpha;
  pred_func_t<K_IMP::YES> m_pred_f;             //function to make predictions
  
public:
  
  //constructor
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
   error_single_cv_iter(const int &param, const KO_Traits::StoringMatrix &training_set, const KO_Traits::StoringMatrix &validation_set)
   const
   {
      //train the model to make prediction
      //auto pred = m_pred_f(training_set,0.95,m_alpha,param);
      auto pred = m_pred_f(training_set,m_alpha,param,this->number_threads());
      return this->err_valid_set_eval(pred,validation_set,this->number_threads());
   }
   
   
   //validation error for a single parameter
   inline 
   double 
   error_single_param(const int &param, const cv_strategy_t &strat, const std::size_t &number_cv_iter) 
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