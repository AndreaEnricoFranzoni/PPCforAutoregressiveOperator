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
    std::cout << "FIndining them" << std::endl;
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
    std::cout << "Already them" << std::endl;
    Spectra::SymEigsSolver<Spectra::DenseSymMatProd<double>> eigsolver_phi(op, m_k, 2*m_k);
    eigsolver_phi.init();
    int nconv = eigsolver_phi.compute(Spectra::SortRule::LargestAlge);
    
    return std::make_tuple(m_k,eigsolver_phi.eigenvalues(),eigsolver_phi.eigenvectors());
  }
}














//doing the inverse square root of a spd matrix
//template< class D, DOM_DIM dom_dim, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval >
//KO_Traits::StoringMatrix
//PPC_KO_base<D, dom_dim, k_imp, valid_err_ret, cv_strat, cv_err_eval>::matrix_inverse_square_root(const KO_Traits::StoringMatrix& mat)
//const
//{   
  //covariance matrix: is symmetric. Since only real values: self-adjoint
  //expoliting Eigen library to do spectral decomposition efficiently
//  Eigen::SelfAdjointEigenSolver<KO_Traits::StoringMatrix> eigensolver(mat);
  
  //if (eigensolver.info() != Eigen::Success) abort();  //TODO: THROW EXCEPTION
  
  //eigenvalues (descending order)
 // KO_Traits::StoringVector eigvals = eigensolver.eigenvalues().reverse();
  
  //eigenvectors (ordererd wrt descending order of their eigenvalues)
  //KO_Traits::StoringMatrix eigvcts = eigensolver.eigenvectors().rowwise().reverse();
  
  //spectral theorem
  //std::transform(eigvals.begin(),eigvals.end(),
    //             eigvals.begin(),
      //           [](double s){return static_cast<double>(1)/std::sqrt(s);});
  
  //spectral theorem for inverse square root    
  //return eigvcts*(eigvals.asDiagonal())*(eigvcts.transpose());
//}








//KO algorithm (all parameters been set)
template< class D, DOM_DIM dom_dim, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval >
void
PPC_KO_base<D, dom_dim, k_imp, valid_err_ret, cv_strat, cv_err_eval>::KO_algo()
{ 
  //Square root inverse of reg covariance
  //m_CovRegRoot = this->matrix_inverse_square_root(m_CovReg);
  Eigen::SelfAdjointEigenSolver<KO_Traits::StoringMatrix> eigensolver_cov_reg(m_CovReg);
  KO_Traits::StoringMatrix m_CovRegRoot = eigensolver_cov_reg.operatorInverseSqrt();
  
  
  
  //Phi hat
  KO_Traits::StoringMatrix phi_hat = m_CovRegRoot*m_GammaSquared*m_CovRegRoot.transpose();
  double tot_exp_pow = phi_hat.trace();
  
  //Spectral decomposition of phi_hat: self-adjoint, so exploiting it 
  //Eigen::SelfAdjointEigenSolver<KO_Traits::StoringMatrix> eigensolver_phi(phi_hat);
  //KO_Traits::StoringVector eigvals_phi = eigensolver_phi.eigenvalues().reverse();
  /*
   * //retaining only the first k components of eigenvalues and eigenvectors
   const KO_Traits::StoringVector D_hat = eigvals_phi.head(m_k);
   const KO_Traits::StoringMatrix V_hat = eigensolver_phi.eigenvectors().rowwise().reverse().leftCols(m_k);
   
   
   std::cout << "Eigen:\n" << D_hat << std::endl;
   std::cout << "Eigen:\n" << V_hat << std::endl;
   
   
   //explanatory power
   m_explanatory_power.resize(m_k);
   std::partial_sum(D_hat.begin(),D_hat.end(),m_explanatory_power.begin());        
   std::for_each(m_explanatory_power.begin(),m_explanatory_power.end(),[&tot_exp_pow](auto &el){el=el/tot_exp_pow;});
   
   */
  
  
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
  
  //Predictive factors (b_i)
  m_b = m_CovRegRoot*std::get<2>(ppcs_ret);

  
  //Predictive factors (b_i)
  //m_b = m_CovRegRoot*V_hat;
  //Predictive loadings (a_i)
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
  
  //se bisogna fare con il prodotto scalare con l'integrale
  //PPC::scores<dom_dim> scores_(m_X.col(m_n-1),m_a,grid_eval);
  //scores_.evaluating_scores();
  //return scores_.scores_evaluations();
}