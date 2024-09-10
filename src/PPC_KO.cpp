#include "PPC_KO.hpp"
#include <execution>

/**********************************
 ********    BASE    ***************
 ***********************************/

int
PPC::PPC_KO_base::PPC_retained(const KO_Traits::StoringArray& eigvals)
const
{
  KO_Traits::StoringArray explained_power(this->m());         //as Eigen::array to perform well division by scalar
  std::partial_sum(eigvals.begin(),eigvals.end(),explained_power.begin());        //cumulative prediction power
  
  explained_power = explained_power/explained_power(this->m()-1); //normalizing
  
  //retaining the number of components that explained a fixed threshold
  int p = std::distance(explained_power.begin(),std::find_if(explained_power.begin(),explained_power.end(),[this](double s){return s > this->p_threshold();})) + static_cast<int>(1);
  
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
  
  
  
  if(this->p_imposed())       //I impose the number of components to invert cov reg
  { 
    
    m_p = this->k();    //TODO: CHECK THAT IF P_IMPOSED, K >= 1
    std::transform(eigvals.head(m_p).begin(),eigvals.head(m_p).end(),
                   eigvals.head(m_p).begin(),
                   [](double s){return static_cast<double>(1)/std::sqrt(s);});
    
    //spectral theorem for inverse square root    
    return (eigvcts.leftCols(m_p))*(eigvals.head(m_p).asDiagonal())*(eigvcts.leftCols(m_p).transpose());
  }
  else                        //I have to find the numbe rof components to invert cov reg
  {

    if(m_p_as_k==true)
    {
      //if PPCs retained are chosen in the square root inverse of cov reg
      int p = this->PPC_retained(eigvals);
      m_p = p;
      //if PPCs retained are chosen from phi
      //int p = this->m();
      
      //taking inverse square root of the eigenvalues as diagonal matrix (only of the p biggest)
      //NB: PAR, DA SISTEMARE PER R (std::execution::par non funziona)
      std::transform(eigvals.head(p).begin(),eigvals.head(p).end(),
                     eigvals.head(p).begin(),
                     [](double s){return static_cast<double>(1)/std::sqrt(s);});
      
      //spectral theorem for inverse square root    
      return (eigvcts.leftCols(p))*(eigvals.head(p).asDiagonal())*(eigvcts.leftCols(p).transpose());
      
    }
    else if(m_p_as_k==false)
    { 
      m_p = this->m();
      //taking inverse square root of the eigenvalues as diagonal matrix (only of the p biggest)
      //NB: PAR, DA SISTEMARE PER R (std::execution::par non funziona)
      std::transform(eigvals.begin(),eigvals.end(),
                     eigvals.begin(),
                     [](double s){return static_cast<double>(1)/std::sqrt(s);});
      
      
      //spectral theorem for inverse square root    
      return eigvcts*(eigvals.asDiagonal())*(eigvcts.transpose());
    }
  }
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
  
  if(!m_k_imposed)
  {
    if(m_p_as_k==true)
    {
      this->k() = this->p();
    }
    else if(m_p_as_k==false)
    {
      //PPCs retained from phi
      int p = this->PPC_retained(eigvals_phi);
      //setting number of retained PPCs
      this->k() = p;
    }
  }



  /*
   * //PPCs retained from phi
   int p = this->PPC_retained(eigvals_phi);
   //setting number of retained PPCs
   this->k() = p;
   */
  
  
  
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
 ********  CV alpha   ***************
 ***********************************/



double
PPC::KO_CV_alpha::alpha_best_CV() 
{ 
  //passing data that are non normalized to make CV
  CV_PPC::CV_KO<double> ko_cv(std::move(this->X_non_norm()),this->n_disc(), m_alpha_min, m_alpha_max, this->p_threshold(), this->p_as_k(), this->p_imposed(), 0.75, this->k(), this->cv_iter_f(), this->ef());
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
  return ko_cv.param_best();
}



void
PPC::KO_CV_alpha::solve()
{
  
  this->alpha() = this->alpha_best_CV();      //finding the best alpha using CV
  
  //only the evaluation of the regularized covariance was missing since there was not any regularization parameter before
  this->CovReg() = this->Cov().array() + this->alpha()*(KO_Traits::StoringMatrix::Identity(this->m(),this->m()).array());
  
  this->KO_algo(); 
}




/**********************************
 ********  CV k   ***************
 ***********************************/

int
PPC::KO_CV_k::k_best_CV()
{ 
  //passing data that are non normalized to make CV
  CV_PPC::CV_KO<int> ko_cv(std::move(this->X_non_norm()),this->X().rows(), static_cast<int>(1), static_cast<int>(this->m()), this->p_threshold(), this->p_as_k(), this->p_imposed(), this->alpha(), static_cast<int>(1), this->cv_iter_f(), this->ef());
  //doing Cv to retrieve the best alpha
  ko_cv.best_param();
  
  
  
  //da qui
  this->ValidErr().reserve(this->X().rows());
  for (size_t i = 0; i < ko_cv.errors().size(); ++i)
  {
    this->ValidErr().push_back(ko_cv.errors()[i]);
  }
  ko_cv.errors().clear();
  //a qui serve solo se si vogliono salvare gli errori di validazione anche nella parte finale:
  //sennò, togliere questa parte, rimettere questa funzione const, e togliere il membro privato dalla classe padre astratta
  
  
  
  //returning the best alpha
  return ko_cv.param_best();
}


void
PPC::KO_CV_k::solve()
{
  this->k() = this->k_best_CV();
  this->k_imposed() = true;
  
  this->KO_algo();
}