#ifndef KO_PPC_CV_K_CRTP_HPP
#define KO_PPC_CV_K_CRTP_HPP

#include "PPC_KO.hpp"
#include "PPC_KO_NoCV.hpp"


//class for doing CV for k
template< DOM_DIM dom_dim, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval >
class PPC_KO_CV_k : public PPC_KO_base<PPC_KO_CV_k<dom_dim,k_imp,valid_err_ret,cv_strat,cv_err_eval>,dom_dim,k_imp,valid_err_ret,cv_strat,cv_err_eval>
{
private:
  std::vector<int> m_k_s;
  KO_Traits::StoringMatrix m_X_non_cent; 
  double m_toll;
  int min_size_ts;
  int max_size_ts;
  
public:
  
  template<typename STOR_OBJ>
  PPC_KO_CV_k(STOR_OBJ&& X, std::vector<int> &k_s, double alpha, double toll, int min_size_ts, int max_size_ts) 
    : 
    PPC_KO_base<PPC_KO_CV_k,dom_dim,k_imp,valid_err_ret,cv_strat,cv_err_eval>(std::move(X)),
    m_k_s(k_s),
    m_X_non_cent(this->m(),this->n()),
    m_toll(toll),
    m_min_size_ts(min_size_ts),
    m_max_size_ts(max_size_ts)
    {
      this->alpha() = alpha;
      this->CovReg() = this->Cov().array() + this->alpha()*this->trace_cov()*(KO_Traits::StoringMatrix::Identity(this->m(),this->m()).array());
      
      for (size_t i = 0; i < this->n(); ++i)
      {
        m_X_non_cent.col(i) = this->X().col(i).array() + this->means();
      }
      
      
       /*
        * std::cout << "PPC KO CV k CRTP: dim dom: " << dom_dim << ", k imp: " << k_imp << ", ver: " << valid_err_ret << std::endl;
        std::cout << "My data centered have " << this->X().rows() << " rows and " << this->X().cols() << " cols" << std::endl;
        std::cout << this->X() << std::endl;
        std::cout << "My data non centered have " << m_X_non_cent.rows() << " rows and " << m_X_non_cent.cols() << " cols" << std::endl;
        std::cout << m_X_non_cent << std::endl;
        std::cout << "alpha:" << this->alpha() <<  std::endl;
        std::cout << "ks:" << std::endl;
        for(std::size_t i = 0; i < m_k_s.size(); ++i){std::cout<<m_k_s[i]<<std::endl;}
        std::cout << "PPC KO CV k END" << std::endl;
        */
       
    }
  
  
  
  inline
  void 
  solving()
  {
  
    //scrivere una l'oggetto CV strategy
    auto strategy_cv = Factory_cv_strat<cv_strat>::cv_strat_obj(m_min_size_ts,m_max_size_ts);

    for (size_t i = 0; i < (*strategy_cv).strategy().size(); ++i)
    {
      std::cout << "Train: " << (*strategy_cv).strategy()[i].first.front() << std::endl;
      std::cout << "Valid: " << (*strategy_cv).strategy()[i].second.front() << std::endl;
    }
    
    
    //lambda wrapper for the correct overload
    auto predictor = [](KO_Traits::StoringMatrix&& data, double alpha, int k) { return cv_pred_func<DOM_DIM::uni_dim,K_IMP::YES,VALID_ERR_RET::NO_err,CV_STRAT::AUGMENTING_WINDOW,CV_ERR_EVAL::MSE>(std::move(data),alpha,k);};
    
    double toll_param = m_toll*this->trace_cov();
    //CV_PPC::CV_KO_PPC_k_CRTP<cv_strat,cv_err_eval,valid_err_ret> cv(std::move(m_X_non_cent),std::move(*strategy_cv),m_k_s,m_toll,this->alpha(),PPC::ko_single_cv_CRTP<dom_dim,K_IMP::YES,valid_err_ret,cv_strat,cv_err_eval>);
    CV_k<cv_strat,cv_err_eval,K_IMP::YES,valid_err_ret> cv(std::move(m_X_non_cent),std::move(*strategy_cv),m_k_s,toll_param,this->alpha(),predictor);
    
    cv.best_param_search();
    this->k() = cv.param_best();
    
    //da qui
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
