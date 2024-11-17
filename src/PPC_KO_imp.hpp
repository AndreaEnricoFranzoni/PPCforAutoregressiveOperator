#include "PPC_KO.hpp"

#include <Eigen/Core>
#include "spectra/include/Spectra/MatOp/DenseSymMatProd.h"
#include "spectra/include/Spectra/SymEigsSolver.h"


//how many PPCs have to be retained
template< class D, DOM_DIM dom_dim, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval >
std::tuple<int,KO_Traits::StoringVector,KO_Traits::StoringMatrix>
PPC_KO_base<D, dom_dim, k_imp, valid_err_ret, cv_strat, cv_err_eval>::PPC_retained(const KO_Traits::StoringMatrix& phi, double tot_exp_pow)
const
{
  Spectra::DenseSymMatProd<double> op(phi);
  
  if constexpr( k_imp == K_IMP::NO )
  {
    for(std::size_t i = 0; i < m_m; ++i)
    {
      int n_ppcs = i+1;
      Spectra::SymEigsSolver<Spectra::DenseSymMatProd<double>> eigsolver_phi(op, n_ppcs, 2*n_ppcs);
      eigsolver_phi.init();
      int nconv = eigsolver_phi.compute(Spectra::SortRule::LargestAlge);
      
      if(eigsolver_phi.eigenvalues().sum()/tot_exp_pow >= m_threshold_ppc)
      {
        return std::make_tuple(n_ppcs,eigsolver_phi.eigenvalues(),eigsolver_phi.eigenvectors());
      }
    }
  } 
  else
  {
    Spectra::SymEigsSolver<Spectra::DenseSymMatProd<double>> eigsolver_phi(op, m_k, 2*m_k);
    eigsolver_phi.init();
    int nconv = eigsolver_phi.compute(Spectra::SortRule::LargestAlge);
    
    return std::make_tuple(m_k,eigsolver_phi.eigenvalues(),eigsolver_phi.eigenvectors());
  }
}



//KO algorithm (all parameters been set)
template< class D, DOM_DIM dom_dim, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval >
void
PPC_KO_base<D, dom_dim, k_imp, valid_err_ret, cv_strat, cv_err_eval>::KO_algo()
{ 
  //Square root inverse of reg covariance: self-adjoint:exploiting it
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