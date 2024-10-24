#ifndef KO_PPC_NOCV_CRTP_HPP
#define KO_PPC_NOCV_CRTP_HPP

#include "PPC_KO.hpp"


//class for doing KO given paramters
template< DOM_DIM dom_dim, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval > 
class PPC_KO_NoCV : public PPC_KO_base<PPC_KO_NoCV<dom_dim,k_imp,valid_err_ret,cv_strat,cv_err_eval>,dom_dim,k_imp,valid_err_ret,cv_strat,cv_err_eval>
{
public:
  
  //k already known
  template<typename STOR_OBJ>
  PPC_KO_NoCV(STOR_OBJ&& X, double alpha, int k)
    :   PPC_KO_base<PPC_KO_NoCV,dom_dim,k_imp,valid_err_ret,cv_strat,cv_err_eval>(std::move(X))
    { 
      this->alpha() = alpha;
      this->k() = k;
      this->CovReg() = this->Cov().array() + this->alpha()*this->trace_cov()*(KO_Traits::StoringMatrix::Identity(this->m(),this->m()).array());
      
      
      
       /*
        * std::cout << "PPC KO NoCV BEGIN: dim dom: " << dom_dim << ", k imp: " << k_imp << ", ver: " << valid_err_ret << std::endl;
        std::cout << "My data has " << this->X().rows() << " rows and " << this->X().cols() << " cols" << std::endl;
        std::cout << this->X() << std::endl;
        std::cout << "alpha: " << this->alpha() << ", k:" << this->k() <<",thresh: " << this->threshold_ppc() <<  std::endl;
        std::cout << "PPC KO NoCV END: dim dom: " << std::endl;
        
        */
      
    }
  
  //k to be found with explanatory power
  template<typename STOR_OBJ>
  PPC_KO_NoCV(STOR_OBJ&& X, double alpha, double threshold_ppc)
    :   PPC_KO_base<PPC_KO_NoCV,dom_dim,k_imp,valid_err_ret,cv_strat,cv_err_eval>(std::move(X))
    {
      this->alpha() = alpha;
      this->threshold_ppc() = threshold_ppc;
      this->CovReg() = this->Cov().array() + this->alpha()*this->trace_cov()*(KO_Traits::StoringMatrix::Identity(this->m(),this->m()).array());
      
      
       /*
        * std::cout << "PPC KO NoCV BEGIN: dim dom: " << dom_dim << ", k imp: " << k_imp << ", ver: " << valid_err_ret << std::endl;
        std::cout << "My data has " << this->X().rows() << " rows and " << this->X().cols() << " cols" << std::endl;
        std::cout << this->X() << std::endl;
        std::cout << "alpha: " << this->alpha() << ", k:" << this->k() <<",thresh: " << this->threshold_ppc() <<  std::endl;
        std::cout << "PPC KO NoCV END: dim dom: " << std::endl;
        */
       
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
template< DOM_DIM dom_dim, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval >
KO_Traits::StoringVector 
cv_pred_func(KO_Traits::StoringMatrix && training_data, double alpha, int k)
{
  std::cout << "k imp: " << k << std::endl;
  
  PPC_KO_NoCV<dom_dim,k_imp,valid_err_ret,cv_strat,cv_err_eval> iter(std::move(training_data),alpha,k);
  iter.solving();
  
  return iter.prediction(); 
};


//overloading if k has to be found
template< DOM_DIM dom_dim, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval >
KO_Traits::StoringVector 
cv_pred_func(KO_Traits::StoringMatrix && training_data, double alpha, double threshold_ppc)
{
  std::cout << "thresh: " << threshold_ppc << std::endl;
  
  PPC_KO_NoCV<dom_dim,k_imp,valid_err_ret,cv_strat,cv_err_eval> iter(std::move(training_data),alpha,threshold_ppc);
  iter.solving();
  
  return iter.prediction(); 
};

#endif  //KO_PPC_NOCV_CRTP_HPP