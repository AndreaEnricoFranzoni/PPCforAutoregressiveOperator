#include "PPC_KO_wrapper.hpp"

#include <iostream>
#include <algorithm>

//no cv
template< SOLVER solver, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval >
void
PPC_KO_wrapper_no_cv<solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>::call_ko()
{
  
  if constexpr(k_imp == K_IMP::YES)   //k already known
  {
    PPC_KO_NoCV<solver,K_IMP::YES,valid_err_ret,cv_strat,cv_err_eval> KO(std::move(this->data()),m_alpha,m_k,this->number_threads());
    
    KO.solve(); 
    auto scores = KO.scores();
    auto sd_scores = KO.sd_scores_dir_wei();
    
    if constexpr( valid_err_ret == VALID_ERR_RET::YES_err){this->results() = std::make_tuple(KO.prediction(),KO.alpha(),KO.k(),scores,KO.explanatory_power(),KO.a(),KO.b(),sd_scores,KO.means(),KO.ValidErr());}
    else  {this->results() = std::make_tuple(KO.prediction(),KO.alpha(),KO.k(),scores,KO.explanatory_power(),KO.a(),KO.b(),sd_scores,KO.means());}
  }
  if constexpr(k_imp == K_IMP::NO)    //k to be found with explanatory power
  {
    PPC_KO_NoCV<solver,K_IMP::NO,valid_err_ret,cv_strat,cv_err_eval> KO(std::move(this->data()),m_alpha,m_threshold_ppc,this->number_threads());
    
    KO.solve();
    auto scores = KO.scores();
    auto sd_scores = KO.sd_scores_dir_wei();
    
    if constexpr( valid_err_ret == VALID_ERR_RET::YES_err){this->results() = std::make_tuple(KO.prediction(),KO.alpha(),KO.k(),scores,KO.explanatory_power(),KO.a(),KO.b(),sd_scores,KO.means(),KO.ValidErr());}
    else  {this->results() = std::make_tuple(KO.prediction(),KO.alpha(),KO.k(),scores,KO.explanatory_power(),KO.a(),KO.b(),sd_scores,KO.means());}
  }
}
  
  

//cv alpha
template< SOLVER solver, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval >
void
PPC_KO_wrapper_cv_alpha<solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>::call_ko()
{
  
  if constexpr(k_imp == K_IMP::YES)   //k known
  {
    PPC_KO_CV_alpha<solver,K_IMP::YES,valid_err_ret,cv_strat,cv_err_eval> KO(std::move(this->data()),m_alphas,m_k,m_min_size_ts,m_max_size_ts,this->number_threads());
    
    KO.solve();
    auto scores = KO.scores();
    auto sd_scores = KO.sd_scores_dir_wei();
    
    
    if constexpr( valid_err_ret == VALID_ERR_RET::YES_err){this->results() = std::make_tuple(KO.prediction(),KO.alpha(),KO.k(),scores,KO.explanatory_power(),KO.a(),KO.b(),sd_scores,KO.means(),std::get<valid_err_cv_1_t>(KO.ValidErr()));}
    else  {this->results() = std::make_tuple(KO.prediction(),KO.alpha(),KO.k(),scores,KO.explanatory_power(),KO.a(),KO.b(),sd_scores,KO.means());}
  }
  if constexpr(k_imp == K_IMP::NO)    //k to be found with explanatory power
  {
    std::cout << "Wrapper: " << m_threshold_ppc << std::endl;
    PPC_KO_CV_alpha<solver,K_IMP::NO,valid_err_ret,cv_strat,cv_err_eval> KO(std::move(this->data()),m_alphas,m_threshold_ppc,m_min_size_ts,m_max_size_ts,this->number_threads());
    
    KO.solve();
    auto scores = KO.scores();
    auto sd_scores = KO.sd_scores_dir_wei();
    
    if constexpr( valid_err_ret == VALID_ERR_RET::YES_err){this->results() = std::make_tuple(KO.prediction(),KO.alpha(),KO.k(),scores,KO.explanatory_power(),KO.a(),KO.b(),sd_scores,KO.means(),std::get<valid_err_cv_1_t>(KO.ValidErr()));}
    else  {this->results() = std::make_tuple(KO.prediction(),KO.alpha(),KO.k(),scores,KO.explanatory_power(),KO.a(),KO.b(),sd_scores,KO.means());}  
  }
}



//cv k
template< SOLVER solver, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval >
void
PPC_KO_wrapper_cv_k<solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>::call_ko()
{
  PPC_KO_CV_k<solver,K_IMP::YES,valid_err_ret,cv_strat,cv_err_eval> KO(std::move(this->data()),m_k_s,m_alpha,m_toll,m_min_size_ts,m_max_size_ts,this->number_threads());
  
  KO.solve();
  auto scores = KO.scores();
  auto sd_scores = KO.sd_scores_dir_wei();
  
  if constexpr( valid_err_ret == VALID_ERR_RET::YES_err){this->results() = std::make_tuple(KO.prediction(),KO.alpha(),KO.k(),scores,KO.explanatory_power(),KO.a(),KO.b(),sd_scores,KO.means(),std::get<valid_err_cv_1_t>(KO.ValidErr()));}
  else  {this->results() = std::make_tuple(KO.prediction(),KO.alpha(),KO.k(),scores,KO.explanatory_power(),KO.a(),KO.b(),sd_scores,KO.means());}
}



//cv alpha-k
template< SOLVER solver, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval >
void
PPC_KO_wrapper_cv_alpha_k<solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>::call_ko()
{
  PPC_KO_CV_alpha_k<solver,K_IMP::YES,valid_err_ret,cv_strat,cv_err_eval> KO(std::move(this->data()),m_alphas,m_k_s,m_toll,m_min_size_ts,m_max_size_ts,this->number_threads());
  
  KO.solve();
  auto scores = KO.scores();
  auto sd_scores = KO.sd_scores_dir_wei();
  
  
  if constexpr( valid_err_ret == VALID_ERR_RET::YES_err){this->results() = std::make_tuple(KO.prediction(),KO.alpha(),KO.k(),scores,KO.explanatory_power(),KO.a(),KO.b(),sd_scores,KO.means(),std::get<valid_err_cv_2_t>(KO.ValidErr()));}
  else  {this->results() = std::make_tuple(KO.prediction(),KO.alpha(),KO.k(),scores,KO.explanatory_power(),KO.a(),KO.b(),sd_scores,KO.means());}
}
