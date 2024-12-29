#ifndef KO_PPC_CV_K_CRTP_HPP
#define KO_PPC_CV_K_CRTP_HPP

#include "PPC_KO.hpp"
#include "PPC_KO_NoCV.hpp"


//class for doing CV for k
template< SOLVER solver, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval >
class PPC_KO_CV_k : public PPC_KO_base<PPC_KO_CV_k<solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>,solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>
{
private:
  std::vector<int> m_k_s;
  KO_Traits::StoringMatrix m_X_non_cent; 
  double m_toll;
  int m_min_size_ts;
  int m_max_size_ts;
  
public:
  
  template<typename STOR_OBJ>
  PPC_KO_CV_k(STOR_OBJ&& X, std::vector<int> &k_s, double alpha, double toll, int min_size_ts, int max_size_ts, int number_threads) 
    : 
    PPC_KO_base<PPC_KO_CV_k,solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>(std::move(X),number_threads),
    m_k_s(k_s),
    m_X_non_cent(this->m(),this->n()),
    m_toll(toll),
    m_min_size_ts(min_size_ts),
    m_max_size_ts(max_size_ts)
    {
      this->alpha() = alpha;
      this->CovReg() = this->Cov().array() + this->alpha()*this->trace_cov()*(KO_Traits::StoringMatrix::Identity(this->m(),this->m()).array());
      
#ifdef _OPENMP
#pragma omp parallel for num_threads(this->number_threads())
#endif
      for (size_t i = 0; i < this->n(); ++i)
      {
        m_X_non_cent.col(i) = this->X().col(i).array() + this->means();
      }
    }
  
  
  
  inline
  void 
  solving()
  {
  
    //scrivere una l'oggetto CV strategy
    auto strategy_cv = Factory_cv_strat<cv_strat>::cv_strat_obj(m_min_size_ts,m_max_size_ts);

    //lambda wrapper for the correct overload
    auto predictor = [](KO_Traits::StoringMatrix&& data, double alpha, int k, int number_threads) { return cv_pred_func<solver,K_IMP::YES,VALID_ERR_RET::NO_err,cv_strat,cv_err_eval>(std::move(data),alpha,k,number_threads);};
    
    //to stop the algorithm if not too much difference in adding a PPC
    double toll_param = m_toll*this->trace_cov();
    
    //CV for k
    CV_k<cv_strat,cv_err_eval,K_IMP::YES,valid_err_ret> cv(std::move(m_X_non_cent),std::move(*strategy_cv),m_k_s,toll_param,this->alpha(),predictor,this->number_threads());
    
    //best number of PPCs
    cv.best_param_search();
    this->k() = cv.param_best();
    
    //if errors to be saved
    if constexpr (valid_err_ret == VALID_ERR_RET::YES_err)
    {
      valid_err_cv_1_t m;
      m.reserve(cv.valid_errors().size());
      for (size_t i = 0; i < cv.valid_errors().size(); ++i)  { m.emplace_back(cv.valid_errors()[i]);}
      this->ValidErr() = m;
      m.clear();
    }
    
    this->KO_algo(); 
  };
};


#endif  //KO_PPC_CV_K_CRTP_HPP