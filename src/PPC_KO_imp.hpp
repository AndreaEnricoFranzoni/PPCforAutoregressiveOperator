#include "PPC_KO.hpp"

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include "spectra/include/Spectra/MatOp/DenseSymMatProd.h"
#include "spectra/include/Spectra/MatOp/DenseCholesky.h"
#include "spectra/include/Spectra/SymEigsSolver.h"
#include "spectra/include/Spectra/SymGEigsSolver.h"





//how many PPCs have to be retained
template< class D, DOM_DIM dom_dim, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval >
std::tuple<int,KO_Traits::StoringVector,KO_Traits::StoringMatrix>
PPC_KO_base<D, dom_dim, k_imp, valid_err_ret, cv_strat, cv_err_eval>::PPC_retained()
{
  
  //exact method
  if constexpr(dom_dim == DOM_DIM::uni_dim)
  {
    
    //Square root inverse of reg covariance: self-adjoint:exploiting it
    Eigen::SelfAdjointEigenSolver<KO_Traits::StoringMatrix> eigensolver_cov_reg(m_CovReg);
    KO_Traits::StoringMatrix m_CovRegRoot = eigensolver_cov_reg.operatorInverseSqrt();
    
    //Phi hat: self-adjoint:exploiting it
    KO_Traits::StoringMatrix phi_hat = m_CovRegRoot*m_GammaSquared*m_CovRegRoot.transpose();
    m_tot_exp_pow = phi_hat.trace();
    
    std::cout << "Usando metodo diretto" << std::endl;
    Spectra::DenseSymMatProd<double> op(phi_hat);
    
    if constexpr( k_imp == K_IMP::NO )    //number of PPCs to be detected
    {
      for(std::size_t i = 0; i < m_m; ++i)
      {
        int n_ppcs = i+1;
        Spectra::SymEigsSolver<Spectra::DenseSymMatProd<double>> eigsolver_phi(op, n_ppcs, 2*n_ppcs);
        eigsolver_phi.init();
        int nconv = eigsolver_phi.compute(Spectra::SortRule::LargestAlge);
        
        if(eigsolver_phi.eigenvalues().sum()/m_tot_exp_pow >= m_threshold_ppc)
        {
          return std::make_tuple(n_ppcs,eigsolver_phi.eigenvalues(),eigsolver_phi.eigenvectors());
        }
      }
    } 
    else              // number of PPCs already known
    {
      Spectra::SymEigsSolver<Spectra::DenseSymMatProd<double>> eigsolver_phi(op, m_k, 2*m_k);
      eigsolver_phi.init();
      int nconv = eigsolver_phi.compute(Spectra::SortRule::LargestAlge);
      
      return std::make_tuple(m_k,eigsolver_phi.eigenvalues(),eigsolver_phi.eigenvectors());
    }
  }
  //generalized method
  else if constexpr(dom_dim == DOM_DIM::bi_dim)
  {

    std::cout << "Usando metodo generalizzato" << std::endl;
    
    // Construct matrix operation objects using the wrapper classes
    Spectra::DenseSymMatProd<double> op(m_GammaSquared);
    Spectra::DenseCholesky<double>  Bop(m_CovReg);
    m_tot_exp_pow = 100.0;
    
    if constexpr(k_imp == K_IMP::NO)
    {
      for(std::size_t i = 0; i < m_m; ++i)
      {
        int n_ppcs = i+1;
        Spectra::SymGEigsSolver<Spectra::DenseSymMatProd<double>, Spectra::DenseCholesky<double>, Spectra::GEigsMode::Cholesky> eigsolver_ppc(op, Bop, n_ppcs, 2*n_ppcs);
        
        eigsolver_ppc.init();
        int nconv = eigsolver_ppc.compute(Spectra::SortRule::LargestAlge);
        
        if(eigsolver_ppc.eigenvalues().sum()/m_tot_exp_pow >= m_threshold_ppc)
        {
          return std::make_tuple(n_ppcs,eigsolver_ppc.eigenvalues(),eigsolver_ppc.eigenvectors());
        }
      }
    }
    else
    {
      Spectra::SymGEigsSolver<Spectra::DenseSymMatProd<double>, Spectra::DenseCholesky<double>, Spectra::GEigsMode::Cholesky> eigsolver_ppc(op, Bop, m_k, 2*m_k);
      
      eigsolver_ppc.init();
      int nconv = eigsolver_ppc.compute(Spectra::SortRule::LargestAlge);
      
      return std::make_tuple(m_k,eigsolver_ppc.eigenvalues(),eigsolver_ppc.eigenvectors());
    }
  }
}





//KO algorithm (all parameters been set)
template< class D, DOM_DIM dom_dim, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval >
void
PPC_KO_base<D, dom_dim, k_imp, valid_err_ret, cv_strat, cv_err_eval>::KO_algo()
{ 
  std::cout << "SOLVING" << std::endl;
  //finding the PPCs
  auto ppcs_ret = this->PPC_retained();
  
  if constexpr( k_imp == K_IMP::NO)
  {
    m_k = std::get<0>(ppcs_ret);
  }
  
  //explanatory power
  m_explanatory_power.resize(m_k);
  std::partial_sum(std::get<1>(ppcs_ret).begin(),std::get<1>(ppcs_ret).end(),m_explanatory_power.begin());        
  std::for_each(m_explanatory_power.begin(),m_explanatory_power.end(),[this](auto &el){el=el/m_tot_exp_pow;});
  
  //Weights (b_i)
  if constexpr(dom_dim == DOM_DIM::uni_dim){m_b = m_CovRegRoot*std::get<2>(ppcs_ret);}  else{m_b = std::get<2>(ppcs_ret);}
  
  //Directions (a_i)
  m_a = m_CrossCov*m_b;
  //predictor estimate
  m_rho = m_a*(m_b.transpose());
  
  
  
  
  
  
  /*
   *   //Square root inverse of reg covariance: self-adjoint:exploiting it
   Eigen::SelfAdjointEigenSolver<KO_Traits::StoringMatrix> eigensolver_cov_reg(m_CovReg);
   KO_Traits::StoringMatrix m_CovRegRoot = eigensolver_cov_reg.operatorInverseSqrt();
   
   //Phi hat: self-adjoint:exploiting it
   KO_Traits::StoringMatrix phi_hat = m_CovRegRoot*m_GammaSquared*m_CovRegRoot.transpose();
   double tot_exp_pow = phi_hat.trace();
   
   //PPCs retained from phi
   auto ppcs_ret = this->PPC_retained(phi_hat,tot_exp_pow);
   
   if constexpr( k_imp == K_IMP::NO)
   {
   m_k = std::get<0>(ppcs_ret);
   }
   
   //explanatory power
   m_explanatory_power.resize(m_k);
   std::partial_sum(std::get<1>(ppcs_ret).begin(),std::get<1>(ppcs_ret).end(),m_explanatory_power.begin());        
   std::for_each(m_explanatory_power.begin(),m_explanatory_power.end(),[&tot_exp_pow](auto &el){el=el/tot_exp_pow;});
   
   //Weights (b_i)
   m_b = m_CovRegRoot*std::get<2>(ppcs_ret);
   //Directions (a_i)
   m_a = m_CrossCov*m_b;
   //predictor estimate
   m_rho = m_a*(m_b.transpose());
   */
}



//one step ahead prediction
template< class D, DOM_DIM dom_dim, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval >
KO_Traits::StoringArray
PPC_KO_base<D, dom_dim, k_imp, valid_err_ret, cv_strat, cv_err_eval>::prediction()
const 
{
  return (m_rho*m_X.col(m_n-1)).array() + m_means;
}



//evaluating the scores on the PPCs
template< class D, DOM_DIM dom_dim, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval >
std::vector<double>
PPC_KO_base<D, dom_dim, k_imp, valid_err_ret, cv_strat, cv_err_eval>::scores()
const
{ 
  std::vector<double> scores;
  scores.reserve(m_k);
  
  for(std::size_t i = 0; i < m_k; ++i)
  {
    scores.emplace_back((m_X.col(m_n-1)).dot((m_a.col(i))));
  }
  
  return scores;
}


//evaluating the sd of the scores on directions and weights
template< class D, DOM_DIM dom_dim, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval >
std::vector<std::array<double,2>>
PPC_KO_base<D, dom_dim, k_imp, valid_err_ret, cv_strat, cv_err_eval>::sd_scores_dir_wei()
const
{
  std::vector<std::array<double,2>> standard_dev;
  standard_dev.reserve(m_k);
  
  //the number of computed scalar prods is all the instants - 1
  std::size_t n = m_X.cols() - 1;
  
  // for each one of the PPC
  for(std::size_t comp = 0; comp < m_k; ++comp){
    
    std::vector<double> scores_dir;
    scores_dir.reserve(n);
    std::vector<double> scores_wei;
    scores_wei.reserve(n);
    
    //computing the scores of directions (it has to be for the next value) and weights (current value)
    for(std::size_t i = 0; i < n; ++i){
      scores_dir.emplace_back(m_X.col(i+1).dot(m_a.col(comp)));
      scores_wei.emplace_back(m_X.col(i).dot(m_b.col(comp)));
    }
    
    //compute the standard deviation
    // mean
    double mean_dir = std::accumulate(scores_dir.begin(), scores_dir.end(), 0.0)/n;
    double mean_wei = std::accumulate(scores_wei.begin(), scores_wei.end(), 0.0)/n;
    
    // variance 
    double variance_dir = std::transform_reduce(scores_dir.begin(), 
                                                scores_dir.end(), 
                                                0.0,
                                                std::plus{},
                                                [mean_dir](auto el) {
                                                return std::pow(el-mean_dir,2);})/n;
    
    double variance_wei = std::transform_reduce(scores_wei.begin(), 
                                                scores_wei.end(), 
                                                0.0,
                                                std::plus{},
                                                [mean_wei](auto el) {
                                                return std::pow(el-mean_wei,2);})/n;
    
    // standard deviation
    standard_dev.emplace_back(std::array<double, 2>{std::sqrt(variance_dir),std::sqrt(variance_wei)});
    
    scores_dir.clear();
    scores_wei.clear();
  }
  
  return standard_dev;
}