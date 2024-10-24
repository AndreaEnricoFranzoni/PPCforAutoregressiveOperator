#include "PPC_KO_wrapper.hpp"

#include "mesh.hpp"
#include <iostream>
#include <algorithm>

//no cv
template< DOM_DIM dom_dim, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval >
void
PPC_KO_wrapper_no_cv<dom_dim,k_imp,valid_err_ret,cv_strat,cv_err_eval>::call_ko()
{
  //k already known
  if constexpr(k_imp == K_IMP::YES)
  {
    
    
    
    
     /*
      * std::cout << "WRAPPER BEGIN NoCV" << std::endl;
      
      std::cout << "Data" << std::endl;
      std::cout << this->data() << std::endl;
      std::cout << "Alpha: " << m_alpha << std::endl;
      std::cout << "k: " << m_k << ", thresh:" << m_threshold_ppc << std::endl;
      std::cout << "WRAPPER END NoCV" << std::endl;
      */
     
    
    
    PPC_KO_NoCV<dom_dim,k_imp,valid_err_ret,cv_strat,cv_err_eval> KO(std::move(this->data()),m_alpha,m_k);
    
    KO.solve(); 
    
    
    Geometry::Domain1D domain_func_data(0,1);
    Geometry::Mesh1D   grid_func_data(domain_func_data,this->data().rows()-static_cast<int>(1));
    std::cout << "Grid of nodes" << this->grid_func_data().numNodes() << std::endl;
    //std::for_each(this->grid_func_data().begin(),this->grid_func_data().end(),[](auto el){std::cout << el << std::endl;});
    auto scores = KO.scores(grid_func_data);
    
    if constexpr( valid_err_ret == VALID_ERR_RET::YES_err){this->results() = std::make_tuple(KO.prediction(),KO.alpha(),KO.k(),scores,KO.rho(),KO.ValidErr());}
    else  {this->results() = std::make_tuple(KO.prediction(),KO.alpha(),KO.k(),scores,KO.rho());}
    
  }
  //k to be found with explanatory power
  if constexpr(k_imp == K_IMP::NO)
  {
    
    
    
     /*
      * std::cout << "WRAPPER BEGIN NoCV" << std::endl;
      
      std::cout << "Data" << std::endl;
      std::cout << this->data() << std::endl;
      std::cout << "Alpha: " << m_alpha << std::endl;
      std::cout << "k: " << m_k << ", thresh:" << m_threshold_ppc << std::endl;
      std::cout << "WRAPPER END NoCV" << std::endl;
      */
    
    
    
    
    
    
    
    
    PPC_KO_NoCV<dom_dim,k_imp,valid_err_ret,cv_strat,cv_err_eval> KO(std::move(this->data()),m_alpha,m_threshold_ppc);
    
    KO.solve();
    
    
    
    Geometry::Domain1D domain_func_data(0,1);
    Geometry::Mesh1D   grid_func_data(domain_func_data,this->data().rows()-static_cast<int>(1));
    //std::cout << "Grid of nodes" << this->grid_func_data().numNodes() << std::endl;
    //std::for_each(this->grid_func_data().begin(),this->grid_func_data().end(),[](auto el){std::cout << el << std::endl;});
    auto scores = KO.scores(grid_func_data);
    
    if constexpr( valid_err_ret == VALID_ERR_RET::YES_err){this->results() = std::make_tuple(KO.prediction(),KO.alpha(),KO.k(),scores,KO.rho(),KO.ValidErr());}
    else  {this->results() = std::make_tuple(KO.prediction(),KO.alpha(),KO.k(),scores,KO.rho());}
  }
}
  
  

//cv alpha
template< DOM_DIM dom_dim, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval >
void
PPC_KO_wrapper_cv_alpha<dom_dim,k_imp,valid_err_ret,cv_strat,cv_err_eval>::call_ko()
{
  //k known
  if constexpr(k_imp == K_IMP::YES)
  {
    
     /*
      * std::cout << "WRAPPER BEGIN CV ALPHA" << std::endl;
      std::cout << "Data" << std::endl;
      std::cout << this->data() << std::endl;
      std::cout << "Alphas: " << std::endl;
      for(std::size_t i = 0; i < m_alphas.size(); ++i){std::cout<<m_alphas[i]<<std::endl;}
      std::cout << "k " << m_k << std::endl;
      std::cout << "WRAPPER END CV ALPHA" << std::endl;
      */
     
    
    PPC_KO_CV_alpha<dom_dim,k_imp,valid_err_ret,cv_strat,cv_err_eval> KO(std::move(this->data()),m_alphas,m_k);
    
    KO.solve();
    
    
    Geometry::Domain1D domain_func_data(0,1);
    Geometry::Mesh1D   grid_func_data(domain_func_data,this->data().rows()-static_cast<int>(1));
    //std::cout << "Grid of nodes" << this->grid_func_data().numNodes() << std::endl;
    //std::for_each(this->grid_func_data().begin(),this->grid_func_data().end(),[](auto el){std::cout << el << std::endl;});
    auto scores = KO.scores(grid_func_data);
    
    
    if constexpr( valid_err_ret == VALID_ERR_RET::YES_err){this->results() = std::make_tuple(KO.prediction(),KO.alpha(),KO.k(),scores,KO.rho(),std::get<valid_err_cv_1_t>(KO.ValidErr()));}
    else  {this->results() = std::make_tuple(KO.prediction(),KO.alpha(),KO.k(),scores,KO.rho());}
  }
  //k to be found with explanatory power
  if constexpr(k_imp == K_IMP::NO)
  {
    
    
    
    
     /*
      * std::cout << "WRAPPER BEGIN CV ALPHA" << std::endl;
      std::cout << "Data" << std::endl;
      std::cout << this->data() << std::endl;
      std::cout << "Alphas: " << std::endl;
      for(std::size_t i = 0; i < m_alphas.size(); ++i){std::cout<<m_alphas[i]<<std::endl;}
      std::cout << "thresh" << m_threshold_ppc << std::endl;
      std::cout << "WRAPPER END CV ALPHA" << std::endl;
      */
     
    
    PPC_KO_CV_alpha<dom_dim,k_imp,valid_err_ret,cv_strat,cv_err_eval> KO(std::move(this->data()),m_alphas,m_threshold_ppc);
    
    KO.solve();
    
    
    Geometry::Domain1D domain_func_data(0,1);
    Geometry::Mesh1D   grid_func_data(domain_func_data,this->data().rows()-static_cast<int>(1));
    //std::cout << "Grid of nodes" << this->grid_func_data().numNodes() << std::endl;
    //std::for_each(this->grid_func_data().begin(),this->grid_func_data().end(),[](auto el){std::cout << el << std::endl;});
    auto scores = KO.scores(grid_func_data);
    
    
    
    if constexpr( valid_err_ret == VALID_ERR_RET::YES_err){this->results() = std::make_tuple(KO.prediction(),KO.alpha(),KO.k(),scores,KO.rho(),std::get<valid_err_cv_1_t>(KO.ValidErr()));}
    else  {this->results() = std::make_tuple(KO.prediction(),KO.alpha(),KO.k(),scores,KO.rho());}  }
}



//cv k
template< DOM_DIM dom_dim, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval >
void
PPC_KO_wrapper_cv_k<dom_dim,k_imp,valid_err_ret,cv_strat,cv_err_eval>::call_ko()
{

  
   /*
    * std::cout << "WRAPPER BEGIN CV K" << std::endl;
    
    std::cout << "Data" << std::endl;
    std::cout << this->data() << std::endl;
    std::cout << "k_s: " << std::endl;
    for(std::size_t i = 0; i < m_k_s.size(); ++i){std::cout<<m_k_s[i]<<std::endl;}
    std::cout << "alpha: " << m_alpha << std::endl;
    std::cout << "WRAPPER END CV K" << std::endl;
    */
   
  
  
  
  
  PPC_KO_CV_k<dom_dim,K_IMP::YES,valid_err_ret,cv_strat,cv_err_eval> KO(std::move(this->data()),m_k_s,m_alpha,m_toll);
  
  KO.solve();
  
  
  Geometry::Domain1D domain_func_data(0,1);
  Geometry::Mesh1D   grid_func_data(domain_func_data,this->data().rows()-static_cast<int>(1));
  //std::cout << "Grid of nodes" << this->grid_func_data().numNodes() << std::endl;
  //std::for_each(this->grid_func_data().begin(),this->grid_func_data().end(),[](auto el){std::cout << el << std::endl;});
  auto scores = KO.scores(grid_func_data);
  
  
  
  if constexpr( valid_err_ret == VALID_ERR_RET::YES_err){this->results() = std::make_tuple(KO.prediction(),KO.alpha(),KO.k(),scores,KO.rho(),std::get<valid_err_cv_1_t>(KO.ValidErr()));}
  else  {this->results() = std::make_tuple(KO.prediction(),KO.alpha(),KO.k(),scores,KO.rho());}
  }



//cv alpha-k
template< DOM_DIM dom_dim, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval >
void
PPC_KO_wrapper_cv_alpha_k<dom_dim,k_imp,valid_err_ret,cv_strat,cv_err_eval>::call_ko()
{
  
  
  
   /*
    * std::cout << "WRAPPER BEGIN CV ALPHA-K" << std::endl;
    std::cout << "Data" << std::endl;
    std::cout << this->data() << std::endl;
    std::cout << "alphas: " << std::endl;
    for(std::size_t i = 0; i < m_alphas.size(); ++i){std::cout<<m_alphas[i]<<std::endl;}
    std::cout << "k_s: " << std::endl;
    for(std::size_t i = 0; i < m_k_s.size(); ++i){std::cout<<m_k_s[i]<<std::endl;}
    std::cout << "WRAPPER END CV ALPHA-K" << std::endl;
    */
   
  
  PPC_KO_CV_alpha_k<dom_dim,K_IMP::YES,valid_err_ret,cv_strat,cv_err_eval> KO(std::move(this->data()),m_alphas,m_k_s,m_toll);
  
  KO.solve();
  
  
  
  Geometry::Domain1D domain_func_data(0,1);
  Geometry::Mesh1D   grid_func_data(domain_func_data,this->data().rows()-static_cast<int>(1));
  //std::cout << "Grid of nodes" << this->grid_func_data().numNodes() << std::endl;
  //std::for_each(this->grid_func_data().begin(),this->grid_func_data().end(),[](auto el){std::cout << el << std::endl;});
  auto scores = KO.scores(grid_func_data);
  
  
  if constexpr( valid_err_ret == VALID_ERR_RET::YES_err){this->results() = std::make_tuple(KO.prediction(),KO.alpha(),KO.k(),scores,KO.rho(),std::get<valid_err_cv_2_t>(KO.ValidErr()));}
  else  {this->results() = std::make_tuple(KO.prediction(),KO.alpha(),KO.k(),scores,KO.rho());}
}