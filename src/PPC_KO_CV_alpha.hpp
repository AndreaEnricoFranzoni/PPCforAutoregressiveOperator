#ifndef KO_PPC_CV_ALPHA_CRTP_HPP
#define KO_PPC_CV_ALPHA_CRTP_HPP

#include "PPC_KO.hpp"
#include "PPC_KO_NoCV.hpp"


//class for doing for alpha
template< DOM_DIM dom_dim, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval >
class PPC_KO_CV_alpha : public PPC_KO_base<PPC_KO_CV_alpha<dom_dim,k_imp,valid_err_ret,cv_strat,cv_err_eval>,dom_dim,k_imp,valid_err_ret,cv_strat,cv_err_eval>
{
private:
  std::vector<double> m_alphas;
  KO_Traits::StoringMatrix m_X_non_cent;
  int m_min_size_ts;
  int m_max_size_ts;  
  
  
public:
  //k already known
  template<typename STOR_OBJ>
  PPC_KO_CV_alpha(STOR_OBJ&& X, const std::vector<double> &alphas, int k, int min_size_ts, int max_size_ts) 
    : 
    PPC_KO_base<PPC_KO_CV_alpha,dom_dim,k_imp,valid_err_ret,cv_strat,cv_err_eval>(std::move(X)),
    m_alphas(alphas),
    m_X_non_cent(this->m(),this->n()),
    m_min_size_ts(min_size_ts),
    m_max_size_ts(max_size_ts)
    {
      this->k() = k; 
      
      for (size_t i = 0; i < this->n(); ++i)
      {
        m_X_non_cent.col(i) = this->X().col(i).array() + this->means();
      }
         
      
      
       /*
        * std::cout << "PPC KO CV ALPHA CRTP: dim dom: " << dom_dim << ", k imp: " << k_imp << ", ver: " << valid_err_ret << std::endl;
        std::cout << "My data centered have " << this->X().rows() << " rows and " << this->X().cols() << " cols" << std::endl;
        std::cout << this->X() << std::endl;
        std::cout << "My data non centered have " << this->X().rows() << " rows and " << this->X().cols() << " cols" << std::endl;
        std::cout << m_X_non_cent << std::endl;
        std::cout << "k:" << this->k() <<",thresh: " << this->threshold_ppc() <<  std::endl;
        std::cout << "Alphas:" << std::endl;
        for(std::size_t i = 0; i < m_alphas.size(); ++i){std::cout<<m_alphas[i]<<std::endl;}
        std::cout << "PPC KO CV ALPHA END" << std::endl;
        */
       
      
      
      
      
      
    }
  
  //k to be found
  template<typename STOR_OBJ>
  PPC_KO_CV_alpha(STOR_OBJ&& X, const std::vector<double> &alphas, double threshold_ppc, int min_size_ts, int max_size_ts) 
    : 
    PPC_KO_base<PPC_KO_CV_alpha,dom_dim,k_imp,valid_err_ret,cv_strat,cv_err_eval>(std::move(X)),
    m_alphas(alphas),
    m_X_non_cent(this->m(),this->n()),
    m_min_size_ts(min_size_ts),
    m_max_size_ts(max_size_ts)
    {
      this->threshold_ppc() = threshold_ppc; 
      
      for (size_t i = 0; i < this->n(); ++i)
      {
        m_X_non_cent.col(i) = this->X().col(i).array() + this->means();
      }
      
    
    
      /*
       * std::cout << "PPC KO CV ALPHA CRTP: dim dom: " << dom_dim << ", k imp: " << k_imp << ", ver: " << valid_err_ret << std::endl;
       std::cout << "My data centered have " << this->X().rows() << " rows and " << this->X().cols() << " cols" << std::endl;
       std::cout << this->X() << std::endl;
       std::cout << "My data non centered have " << m_X_non_cent.rows() << " rows and " << m_X_non_cent.cols() << " cols" << std::endl;
       std::cout << m_X_non_cent << std::endl;
       std::cout << "k:" << this->k() <<",thresh: " << this->threshold_ppc() <<  std::endl;
       std::cout << "Alphas:" << std::endl;
       for(std::size_t i = 0; i < m_alphas.size(); ++i){std::cout<<m_alphas[i]<<std::endl;}
       std::cout << "PPC KO CV ALPHA END" << std::endl;
       */
      
      
      
    }
  

  inline
  void 
  solving()
  {
    //factory to create the cv strategy
    auto strategy_cv = Factory_cv_strat<cv_strat>::cv_strat_obj(m_min_size_ts,m_max_size_ts);

    for (size_t i = 0; i < (*strategy_cv).strategy().size(); ++i)
    {
      std::cout << "Train: " << (*strategy_cv).strategy()[i].first.front() << std::endl;
      std::cout << "Valid: " << (*strategy_cv).strategy()[i].second.front() << std::endl;
    }
    
    if constexpr(k_imp == K_IMP::YES)
    {
      //lambda wrapper for the correct overload for predictor
      auto predictor = [](KO_Traits::StoringMatrix&& data, double alpha, int k) { return cv_pred_func<DOM_DIM::uni_dim,K_IMP::YES,VALID_ERR_RET::NO_err,CV_STRAT::AUGMENTING_WINDOW,CV_ERR_EVAL::MSE>(std::move(data),alpha,k);};
      
      //CV already knowing k
      CV_alpha<cv_strat,cv_err_eval,k_imp,valid_err_ret> cv(std::move(m_X_non_cent),std::move(*strategy_cv),m_alphas,this->k(),predictor);
      
      //best alpha
      cv.best_param_search();
      this->alpha() = cv.param_best();      //finding the best alpha using CV
      
      //if I want to save the errors
      if constexpr(valid_err_ret == VALID_ERR_RET::YES_err)
      {
        valid_err_cv_1_t m;
        m.reserve(m_alphas.size());
        for (size_t i = 0; i < m_alphas.size(); ++i)  { m.emplace_back(cv.valid_errors()[i]);}
        this->ValidErr() = m;
        m.clear();
      }
    }
    
    if constexpr(k_imp == K_IMP::NO)
    {
      //lambda wrapper to resolve the overload
      auto predictor = [](KO_Traits::StoringMatrix&& data, double alpha, double threshold_ppc) { return cv_pred_func<DOM_DIM::uni_dim,K_IMP::NO,VALID_ERR_RET::NO_err,CV_STRAT::AUGMENTING_WINDOW,CV_ERR_EVAL::MSE>(std::move(data),alpha,threshold_ppc);};
      
      //CV already with k to be found with explanatory power
      CV_alpha<cv_strat,cv_err_eval,k_imp,valid_err_ret> cv(std::move(m_X_non_cent),std::move(*strategy_cv),m_alphas,this->threshold_ppc(),predictor);
      
      //best alpha
      cv.best_param_search();
      this->alpha() = cv.param_best();      //finding the best alpha using CV
      
      //only if errors are saved
      if constexpr(valid_err_ret == VALID_ERR_RET::YES_err)
      {
        valid_err_cv_1_t m;
        m.reserve(m_alphas.size());
        for (size_t i = 0; i < m_alphas.size(); ++i)  { m.emplace_back(cv.valid_errors()[i]);}
        this->ValidErr() = m;
        m.clear();
      }
    }
    
    //only the evaluation of the regularized covariance was missing since there was not any regularization parameter before
    this->CovReg() = this->Cov().array() + this->alpha()*this->trace_cov()*(KO_Traits::StoringMatrix::Identity(this->m(),this->m()).array());
    
    this->KO_algo(); 
  };
};

#endif  //KO_PPC_CV_ALPHA_CRTP_HPP