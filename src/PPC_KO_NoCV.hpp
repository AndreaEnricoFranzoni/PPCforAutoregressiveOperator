#ifndef KO_PPC_NOCV_CRTP_HPP
#define KO_PPC_NOCV_CRTP_HPP

#include "PPC_KO.hpp"


//class for doing KO given paramters
template< SOLVER solver, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval > 
class PPC_KO_NoCV : public PPC_KO_base<PPC_KO_NoCV<solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>,solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>
{
public:
  
  //k already known
  template<typename STOR_OBJ>
  PPC_KO_NoCV(STOR_OBJ&& X, double alpha, int k, int number_threads)
    :   PPC_KO_base<PPC_KO_NoCV,solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>(std::move(X),number_threads)
    { 
      this->alpha() = alpha;
      this->k() = k;
      this->CovReg() = this->Cov().array() + this->alpha()*this->trace_cov()*(KO_Traits::StoringMatrix::Identity(this->m(),this->m()).array());
    }
  
  //k to be found with explanatory power
  template<typename STOR_OBJ>
  PPC_KO_NoCV(STOR_OBJ&& X, double alpha, double threshold_ppc, int number_threads)
    :   PPC_KO_base<PPC_KO_NoCV,solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>(std::move(X),number_threads)
    {
      this->alpha() = alpha;
      this->threshold_ppc() = threshold_ppc;
      this->CovReg() = this->Cov().array() + this->alpha()*this->trace_cov()*(KO_Traits::StoringMatrix::Identity(this->m(),this->m()).array());
    }
  
  //run KO
  inline 
  void
  solving()
  {
    this->KO_algo(); 
  }
};






//function to be used in CV to make predictions with the training set
//if k is known
template< SOLVER solver, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval >
KO_Traits::StoringVector 
cv_pred_func(KO_Traits::StoringMatrix && training_data, double alpha, int k, int number_threads)
{  
  PPC_KO_NoCV<solver,K_IMP::YES,valid_err_ret,cv_strat,cv_err_eval> iter(std::move(training_data),alpha,k,number_threads);
  iter.solving();
  
  return iter.prediction(); 
};


//overloading if k has to be found
template< SOLVER solver, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval >
KO_Traits::StoringVector 
cv_pred_func(KO_Traits::StoringMatrix && training_data, double alpha, double threshold_ppc, int number_threads)
{  
  std::cout << "Here" << std::endl;
  PPC_KO_NoCV<solver,K_IMP::NO,valid_err_ret,cv_strat,cv_err_eval> iter(std::move(training_data),alpha,threshold_ppc,number_threads);
  iter.solving();
  
  return iter.prediction(); 
};

#endif  //KO_PPC_NOCV_CRTP_HPP