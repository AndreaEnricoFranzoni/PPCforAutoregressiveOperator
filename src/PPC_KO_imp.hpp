#include "PPC_KO.hpp"


//how many PPCs have to be retained
template< class D, DOM_DIM dom_dim, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval >
int
PPC_KO_base<D, dom_dim, k_imp, valid_err_ret, cv_strat, cv_err_eval>::PPC_retained(const KO_Traits::StoringArray& eigvals)
const
{
  KO_Traits::StoringArray explained_power(m_m);         //as Eigen::array to perform well division by scalar
  std::partial_sum(eigvals.begin(),eigvals.end(),explained_power.begin());        //cumulative prediction power
  
  explained_power = explained_power/explained_power(this->m()-1); //normalizing
  
  //retaining the number of components that explained a fixed threshold
  return std::distance(explained_power.begin(),std::find_if(explained_power.begin(),explained_power.end(),[this](double s){return s > this->threshold_ppc();})) + static_cast<int>(1);
}


//doing the inverse square root of a spd matrix
template< class D, DOM_DIM dom_dim, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval >
KO_Traits::StoringMatrix
PPC_KO_base<D, dom_dim, k_imp, valid_err_ret, cv_strat, cv_err_eval>::matrix_inverse_square_root(const KO_Traits::StoringMatrix& mat)
const
{   
  //covariance matrix: is symmetric. Since only real values: self-adjoint
  //expoliting Eigen library to do spectral decomposition efficiently
  Eigen::SelfAdjointEigenSolver<KO_Traits::StoringMatrix> eigensolver(mat);
  
  if (eigensolver.info() != Eigen::Success) abort();  //TODO: THROW EXCEPTION
  
  //eigenvalues (descending order)
  KO_Traits::StoringVector eigvals = eigensolver.eigenvalues().reverse();
  
  //eigenvectors (ordererd wrt descending order of their eigenvalues)
  KO_Traits::StoringMatrix eigvcts = eigensolver.eigenvectors().rowwise().reverse();
  
  //spectral theorem
  std::transform(eigvals.begin(),eigvals.end(),
                 eigvals.begin(),
                 [](double s){return static_cast<double>(1)/std::sqrt(s);});
  
  //spectral theorem for inverse square root    
  return eigvcts*(eigvals.asDiagonal())*(eigvcts.transpose());
}


//KO algorithm (all parameters been set)
template< class D, DOM_DIM dom_dim, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval >
void
PPC_KO_base<D, dom_dim, k_imp, valid_err_ret, cv_strat, cv_err_eval>::KO_algo()
{ 
  //Square root inverse of reg covariance
  m_CovRegRoot = this->matrix_inverse_square_root(m_CovReg);
  
  //Phi hat
  KO_Traits::StoringMatrix phi_hat = m_CovRegRoot*m_GammaSquared*m_CovRegRoot.transpose();
  
  //Spectral decomposition of phi_hat: self-adjoint, so exploiting it 
  Eigen::SelfAdjointEigenSolver<KO_Traits::StoringMatrix> eigensolver_phi(phi_hat);
  KO_Traits::StoringVector eigvals_phi = eigensolver_phi.eigenvalues().reverse();
  
  if constexpr(k_imp == K_IMP::NO)  //if k not imposed form user: to be retrieved from phi
  {
    //PPCs retained from phi
    std::cout << "K WITH EXPL PW"<< std::endl;
    m_k = this->PPC_retained(eigvals_phi);
  }
  else{std::cout << "K IMP"<< std::endl;}
  
  
  //retaining only the first k components of eigenvalues and eigenvectors
  const KO_Traits::StoringVector D_hat = eigvals_phi.head(m_k);
  const KO_Traits::StoringMatrix V_hat = eigensolver_phi.eigenvectors().rowwise().reverse().leftCols(m_k);
  
  
  //Predictive factors (b_i)
  m_b = m_CovRegRoot*V_hat;
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
PPC_KO_base<D, dom_dim, k_imp, valid_err_ret, cv_strat, cv_err_eval>::scores(const Geometry::Mesh1D &grid_eval)
const
{ 
  //constexpr unsigned long int DIM_DOM=1ul;
  //devo calcolarli su dati centrati? Al momento, sto facendo così
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