#include "PPC_KO.hpp"
#include <execution>

/**********************************
 ********    BASE    ***************
 ***********************************/

int
PPC::PPC_KO_base::PPC_retained(const KO_Traits::StoringArray& cov_reg_eigvals)
const
{
  KO_Traits::StoringArray explained_power(this->m());         //as Eigen::array to perform well division by scalar
  std::partial_sum(cov_reg_eigvals.begin(),cov_reg_eigvals.end(),explained_power.begin());        //cumulative prediction power
  
  explained_power = explained_power/explained_power(this->m()-1); //normalizing
  
  //retaining the number of components that explained a fixed threshold
  int p = std::distance(explained_power.begin(),std::find_if(explained_power.begin(),explained_power.end(),[this](double s){return s > this->p_threshold();})) + static_cast<int>(1);
  
  return p;
}


KO_Traits::StoringMatrix
PPC::PPC_KO_base::matrix_inverse_root(const KO_Traits::StoringMatrix& gamma_alpha)
const
{   
  //covariance matrix: is symmetric. Since only real values: self-adjoint
  //expoliting Eigen library to do spectral decomposition efficiently
  Eigen::SelfAdjointEigenSolver<KO_Traits::StoringMatrix> eigensolver(gamma_alpha);
  
  if (eigensolver.info() != Eigen::Success) abort();  //TODO: THROW EXCEPTION
  
  //eigenvalues (descending order)
  KO_Traits::StoringVector eigvals = eigensolver.eigenvalues().reverse();
  
  //eigenvectors (ordererd wrt ascending order of their eigenvalues)
  KO_Traits::StoringMatrix eigvcts = eigensolver.eigenvectors().rowwise().reverse();
  
  
  //taking inverse square root of the eigenvalues as diagonal matrix (only of the p biggest)
  //NB: PAR, DA SISTEMARE PER R (std::execution::par non funziona)
  std::transform(eigvals.begin(),eigvals.end(),
                 eigvals.begin(),
                 [](double s){return static_cast<double>(1)/std::sqrt(s);});
  
  
  //spectral theorem for inverse square root    
  return eigvcts*(eigvals.asDiagonal())*(eigvcts.transpose());
  
  /*
   //VECCHIA VERSIONE, PER TIRARE FUORI LE PPC DA QUI
   //if PPCs retained are chosen in the square root inverse of cov reg
   //int p = this->PPC_retained(eigvals);
   //if PPCs retained are chosen from phi
   int p = this->m();
   
   //taking inverse square root of the eigenvalues as diagonal matrix (only of the p biggest)
   //NB: PAR, DA SISTEMARE PER R (std::execution::par non funziona)
   std::transform(eigvals.head(p).begin(),eigvals.head(p).end(),
   eigvals.head(p).begin(),
   [](double s){return static_cast<double>(1)/std::sqrt(s);});
   
   //spectral theorem for inverse square root    
   return (eigvcts.leftCols(p))*(eigvals.head(p).asDiagonal())*(eigvcts.leftCols(p).transpose());
   */
}


KO_Traits::StoringMatrix
PPC::PPC_KO_base::phi_estimate()
const
{   
  //gamma_alpha rooted inverse is self-adjoint
  return this->CovRegRoot()*this->GammaSquared()*this->CovRegRoot().transpose();
}


void
PPC::PPC_KO_base::KO_algo()
{
  //Square root inverse of reg covariance
  this->CovRegRoot() = this->matrix_inverse_root(this->CovReg());
  
  //Phi hat
  KO_Traits::StoringMatrix phi_hat = phi_estimate();
  
  //Spectral decomposition of phi_hat: self-adjoint, so exploiting it 
  Eigen::SelfAdjointEigenSolver<KO_Traits::StoringMatrix> eigensolver_phi(phi_hat);
  KO_Traits::StoringVector eigvals_phi = eigensolver_phi.eigenvalues().reverse();
  
  //PPCs retained from phi
  int p = this->PPC_retained(eigvals_phi);
  //setting number of retained PPCs
  this->k() = p;
  
  //reatining only the first k components of eigenvalues and eigenvectors
  const KO_Traits::StoringVector D_hat = eigvals_phi.head(this->k());
  const KO_Traits::StoringMatrix V_hat = eigensolver_phi.eigenvectors().rowwise().reverse().leftCols(this->k());
  
  //Predictive factors (b_i)
  this->b() = this->CovRegRoot()*V_hat;
  //Predictive loadings (a_i)
  this->a() = this->CrossCov()*this->b();
  
  //predictor estimate
  this->rho() = this->a()*(this->b().transpose());
}


KO_Traits::StoringArray
PPC::PPC_KO_base::prediction()
  const 
{
  return (this->rho()*this->X().rightCols(1)).array() + this->means();
}



/**********************************
 ********    NO CV   ***************
 ***********************************/

void
PPC::KO_NO_CV::solve()
{
  this->KO_algo();
}



/**********************************
 ********       CV   ***************
 ***********************************/



double
  PPC::KO_CV::alpha_best_CV() 
  {
    //construction of the object to make CV
    //passing data that are non normalized
    CV_PPC::CV_KO ko_cv(std::move(this->X_non_norm()),this->n_disc());
    //doing Cv to retrieve the best alpha
    ko_cv.best_param();
    
    
    
    //da qui
    this->ValidErr().reserve(m_n_disc);
    for (size_t i = 0; i < ko_cv.errors().size(); ++i)
    {
      this->ValidErr().push_back(ko_cv.errors()[i]);
    }
    ko_cv.errors().clear();
    //a qui serve solo se si vogliono salvare gli errori di validazione anche nella parte finale:
    //sennò, togliere questa parte, rimettere questa funzione const, e togliere il membro privato dalla classe padre astratta
    
    
    
    //returning the best alpha
    return ko_cv.alpha_best();
  }



void
PPC::KO_CV::solve()
{
  //finding the best alpha using CV
  double alpha_best = this->alpha_best_CV();
  
  //setting alpha as alpha best
  this->alpha() = alpha_best;
  
  //only the evaluation of the regularized covariance was missing since there was not any regularization parameter before
  this->CovReg() = this->Cov().array() + this->alpha()*(KO_Traits::StoringMatrix::Identity(this->m(),this->m()).array());
  
  this->KO_algo(); 
}