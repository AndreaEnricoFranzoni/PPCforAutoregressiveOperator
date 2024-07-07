#include "PPC_KO.hpp"


/**********************************
 ********    BASE    ***************
 ***********************************/

int
PPC::PPC_KO_base::PPC_retained(const KO_Traits::StoringArray& cov_reg_eigvals)
{
  KO_Traits::StoringArray explained_power(this->m());         //as Eigen::array to perform well division by scalar
  std::partial_sum(cov_reg_eigvals.begin(),cov_reg_eigvals.end(),explained_power.begin());        //cumulative prediction power
  
  explained_power = explained_power/explained_power(this->m()-1); //normalizing
  
  //retaining the number of components that explained a fixed threshold
  int p = std::distance(explained_power.begin(),std::find_if(explained_power.begin(),explained_power.end(),[this](double s){return s > this->p_threshold();})) + static_cast<int>(1);
  
  //setting the number of components
  this->k() = p;
  
  return p;
}


KO_Traits::StoringMatrix
PPC::PPC_KO_base::matrix_inverse_root(const KO_Traits::StoringMatrix& gamma_alpha) 
{   
  //covariance matrix: is symmetric. Since only real values: self-adjoint
  //expoliting Eigen library to do spectral decomposition efficiently
  Eigen::SelfAdjointEigenSolver<KO_Traits::StoringMatrix> eigensolver(gamma_alpha);
  
  if (eigensolver.info() != Eigen::Success) abort();  //TODO: THROW EXCEPTION
  
  //eigenvalues (descending order)
  KO_Traits::StoringVector eigvals = eigensolver.eigenvalues().reverse();
  
  //eigenvectors (ordererd wrt ascending order of their eigenvalues)
  KO_Traits::StoringMatrix eigvcts = eigensolver.eigenvectors().rowwise().reverse();
  
  //retained components: first p
  int p = this->PPC_retained(eigvals);
  
  //taking inverse square root of the eigenvalues as diagonal matrix (only of the p biggest)
  //NB: PAR, DA SISTEMARE PER R (std::execution::par non funziona)
  std::transform(eigvals.head(p).begin(),eigvals.head(p).end(),
                 eigvals.head(p).begin(),
                 [](double s){return static_cast<double>(1)/std::sqrt(s);});
  
  //spectral theorem for inverse square root    
  return (eigvcts.leftCols(p))*(eigvals.head(p).asDiagonal())*(eigvcts.leftCols(p).transpose());
}


KO_Traits::StoringMatrix
const
PPC::PPC_KO_base::phi_estimate(const KO_Traits::StoringMatrix& gamma_alpha_rooted)
{   
  //gamma_alpha rooted inverse is self-adjoint
  return gamma_alpha_rooted*this->GammaSquared()*gamma_alpha_rooted.transpose();
}


void
PPC::PPC_KO_base::KO_algo()
{
  //Square root inverse of reg covariance
  this->CovRegRoot() = this->matrix_inverse_root(this->CovReg());
  
  //Phi hat
  KO_Traits::StoringMatrix phi_hat = phi_estimate(this->CovRegRoot());
  
  //Spectral decomposition of phi_hat: self-adjoint, so exploiting it 
  Eigen::SelfAdjointEigenSolver<KO_Traits::StoringMatrix> eigensolver_phi(phi_hat);
  
  //reatining only the first k components of eigenvalues and eigenvectors
  const KO_Traits::StoringVector D_hat = eigensolver_phi.eigenvalues().reverse().head(this->k());
  const KO_Traits::StoringMatrix V_hat = eigensolver_phi.eigenvectors().rowwise().reverse().leftCols(this->k());
  
  //Predictive factors (b_i)
  this->b() = this->CovRegRoot()*V_hat;
  //Predictive loadings (a_i)
  this->a() = this->CrossCov()*this->b();
  
  //predictor estimate
  this->rho() = this->a()*(this->b().transpose());
}


KO_Traits::StoringArray
PPC::PPC_KO_base::prediction() const
{
  return (this->rho()*this->X().rightCols(1)).array() + this->means();
}

/**********************************
 ********    NO CV   ***************
 ***********************************/

void
PPC::KO_NO_CV::solve()
{
  std::cout << "Kargin-Onatski using regularization parameter equal to " << this->alpha()  << std::endl;
  
  this->KO_algo();
}





/**********************************
 ********       CV   ***************
 ***********************************/
void
PPC::KO_CV::CV()
{
  std::cout << "Exclusive method" << std::endl;
}


void
PPC::KO_CV::solve()
{
  std::cout << "Kargin-Onatski using CV for regularization parameter" << std::endl;
  
  this->CV();
  
  
}

