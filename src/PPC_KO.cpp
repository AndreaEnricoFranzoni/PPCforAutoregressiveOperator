#include "PPC_KO.hpp"


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
  return std::distance(explained_power.begin(),std::find_if(explained_power.begin(),explained_power.end(),[this](double s){return s > this->p_threshold();})) + static_cast<int>(1);
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
  
  //eigenvectors (ordererd wrt descending order of their eigenvalues)
  KO_Traits::StoringMatrix eigvcts = eigensolver.eigenvectors().rowwise().reverse();
  
  //spectral theorem
  std::transform(eigvals.begin(),eigvals.end(),
                 eigvals.begin(),
                 [](double s){return static_cast<double>(1)/std::sqrt(s);});
  
  //spectral theorem for inverse square root    
  return eigvcts*(eigvals.asDiagonal())*(eigvcts.transpose());
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
  
  if(!m_k_imposed)  //if k not imposed form user: to be retrieved from phi
  {
    //PPCs retained from phi
    int number_PPCs = this->PPC_retained(eigvals_phi);
    //setting number of retained PPCs
    this->k() = number_PPCs;
  }

  //retaining only the first k components of eigenvalues and eigenvectors
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



std::vector<double>
PPC::PPC_KO_base::scores()
const
{ 
  //constexpr int N = 250;
  //devo calcolarli su dati centrati?
  PPC::scores scores_(this->X().col(m_n-1),this->a(), DEF_PARAMS_PPC::a_interval, DEF_PARAMS_PPC::b_interval, m_m - static_cast<int>(1));
  scores_.evaluating_scores();
  
  return scores_.scores_evaluations();
  
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
  CV_PPC::CV_KO<double> ko_cv(std::move(this->X_non_norm()),DEF_PARAMS_PPC::dim_grid_alpha, this->p_threshold(), 0.75, this->k(), this->cv_iter_f(), this->ef());
  //doing Cv to retrieve the best alpha
  ko_cv.best_param();
    
    
    
    //da qui
    this->ValidErr().reserve(DEF_PARAMS_PPC::dim_grid_alpha);
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
  this->CovReg() = this->Cov().array() + this->alpha()*this->trace_cov()*(KO_Traits::StoringMatrix::Identity(this->m(),this->m()).array());
  
  this->KO_algo(); 
}





/**********************************
 ********  CV k   ***************
 ***********************************/

int
PPC::KO_CV_k::k_best_CV()
{ 
  //passing data that are non normalized to make CV
  CV_PPC::CV_KO<int> ko_cv(std::move(this->X_non_norm()), this->m(), this->p_threshold(), this->alpha(), static_cast<int>(1), this->cv_iter_f(), this->ef());
  //doing Cv to retrieve the best k
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
  
  
  
  //returning the best k
  return ko_cv.param_best();
}



void
PPC::KO_CV_k::solve()
{
  this->k() = this->k_best_CV();
  this->k_imposed() = true;
  
  this->KO_algo();
}





/**********************************
 ********  CV    *****************
 ***********************************/

std::pair<double,int>
PPC::KO_CV::params_best_CV()
{
  
  CV_PPC::CV_KO_2 ko_cv(std::move(this->X_non_norm()), DEF_PARAMS_PPC::dim_grid_alpha, this->p_threshold(), this->cv_iter_on_k_f());
  
  ko_cv.best_params();
  
  return ko_cv.params_best();
}



void 
PPC::KO_CV::solve()
{
  //find the best parameters alpha and k
  auto best_parameters_KO = this->params_best_CV();
  
  this->alpha() = best_parameters_KO.first;
  this->CovReg() = this->Cov().array() + this->alpha()*this->trace_cov()*(KO_Traits::StoringMatrix::Identity(this->m(),this->m()).array());
  
  this->k() = best_parameters_KO.second;
  this->k_imposed() = true;
  
  this->KO_algo();
}













/**********************************
 ********  CV CRTP ALPHA    *****************
 ***********************************/


void
PPC::KO_CV_CRTP_ALPHA::solve()
{
  std::vector<double> alphas;
  alphas.resize(21);
  
  std::iota(alphas.begin(),alphas.end(),static_cast<double>(DEF_PARAMS_PPC::min_exp_alphas));
  std::transform(alphas.begin(),alphas.end(),alphas.begin(),[](double el){return(pow(static_cast<double>(10),el));});
  
  
  
  
  CV_PPC::CV_KO_PPC_alpha_CRTP<DEF_PARAMS_PPC::cv_strat_type::AUGMENTING_WINDOW,DEF_PARAMS_PPC::cv_err_eval_type::MSE> cv(std::move(this->X_non_cent()), alphas, this->p_threshold(), this->k(), PPC::ko_single_cv);
  cv.best_param_search();
  
  
  
  //da qui
  this->ValidErr().resize(21);
  for (size_t i = 0; i < 21; ++i)
  {
    this->ValidErr()[i]=cv.valid_errors()[i];
  }
  
  
  
  
  
  this->alpha() = cv.param_best();      //finding the best alpha using CV
  
  //only the evaluation of the regularized covariance was missing since there was not any regularization parameter before
  this->CovReg() = this->Cov().array() + this->alpha()*this->trace_cov()*(KO_Traits::StoringMatrix::Identity(this->m(),this->m()).array());
  
  this->KO_algo(); 
}






/**********************************
 ********  CV CRTP k    *****************
 ***********************************/
void
PPC::KO_CV_CRTP_K::solve()
{
  std::vector<int> k_s;
  k_s.resize(this->m());
  std::iota(k_s.begin(),k_s.end(),static_cast<int>(1));
  
  
  CV_PPC::CV_KO_PPC_k_CRTP<DEF_PARAMS_PPC::cv_strat_type::AUGMENTING_WINDOW,DEF_PARAMS_PPC::cv_err_eval_type::MSE> cv(std::move(this->X_non_cent()), k_s, DEF_PARAMS_PPC::toll_cv_k, this->p_threshold(), this->alpha(), PPC::ko_single_cv);
  
  cv.best_param_search();
  
  
  
  
  //da qui
  this->ValidErr().resize(cv.valid_errors().size());
  for (size_t i = 0; i < cv.valid_errors().size(); ++i)
  {
    this->ValidErr()[i]=(cv.valid_errors()[i]);
  }
  
  
  
  
  
  this->k() = cv.param_best();
  this->k_imposed() = true;
  
  this->KO_algo();
}


















/**********************************
 ********  CV CRTP ALPHA k    *****************
 ***********************************/
void 
PPC::KO_CV_CRTP_ALPHA_K::solve()
{
  
  std::vector<double> alphas;
  alphas.resize(21);
  
  std::iota(alphas.begin(),alphas.end(),static_cast<double>(DEF_PARAMS_PPC::min_exp_alphas));
  std::transform(alphas.begin(),alphas.end(),alphas.begin(),[](double el){return(pow(static_cast<double>(10),el));});
  
  std::vector<int> k_s;
  k_s.resize(this->m());
  std::iota(k_s.begin(),k_s.end(),static_cast<int>(1));
  
  
  
  
  
  
  CV_PPC::CV_KO_PPC_alpha_k_CRTP<DEF_PARAMS_PPC::cv_strat_type::AUGMENTING_WINDOW,DEF_PARAMS_PPC::cv_err_eval_type::MSE> cv(std::move(this->X_non_cent()), alphas, k_s, DEF_PARAMS_PPC::toll_cv_k, this->p_threshold(), PPC::ko_single_cv);
  cv.best_param_search();
  
  
  this->alpha() = cv.alpha_best();
  this->CovReg() = this->Cov().array() + this->alpha()*this->trace_cov()*(KO_Traits::StoringMatrix::Identity(this->m(),this->m()).array());
  
  this->k() = cv.k_best();
  this->k_imposed() = true;
  
  this->KO_algo();
}