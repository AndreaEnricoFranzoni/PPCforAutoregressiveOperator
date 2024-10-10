#ifndef KO_PPC_HPP
#define KO_PPC_HPP

#include <iostream>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <execution>
#include <vector>
#include <functional>
#include <utility>

#include "KO_Traits.hpp"
#include "default_parameters.hpp"
#include "error_function.hpp"
#include "CV_KO.hpp"
#include "CV_KO_2.hpp"
#include "scores.hpp"


#include "CV_CRTP_include.hpp"



namespace PPC       //PrincipalPredictiveComponents
{

/*!
 *Class for doing Principal Predictive Components according to Kargin-Onatski method (PPC)
 *The dataframe passed as input contains temporal series stored column-wise:
 *each row is a time series (m is their total)
 *each column is a time instant (n is their total)
 */

class PPC_KO_base
{
private:
  
  std::size_t m_m;                            //number of evaluation of the functional object (m)
  std::size_t m_n;                            //number of time instants (n)
  KO_Traits::StoringMatrix m_X;               //matrix storing time series (m x n)
  KO_Traits::StoringArray m_means;            //vector storing the time series means (m x 1)  
  KO_Traits::StoringMatrix m_Cov;             //matrix estimating the covariance operator (m x m)
  double m_trace_cov;                         //trace of the covariance estimator
  KO_Traits::StoringMatrix m_CrossCov;        //matrix estimating the cross-covariance operator (m x m)
  KO_Traits::StoringMatrix m_CovReg;          //matrix containing the covariance + alpha*trace(cov)*I (m x m)
  KO_Traits::StoringMatrix m_GammaSquared;    //matrix containing the squared of the cross-covariance matrix (m x m)
  KO_Traits::StoringMatrix m_CovRegRoot;      //matrix containing the inverse squared root of regularized covariance (m x m)
  KO_Traits::StoringMatrix m_a;               //matrix containing predictive loadings (each col)  (m x k)
  KO_Traits::StoringMatrix m_b;               //matrix containing predictive factors (each col)   (m x k)        
  KO_Traits::StoringMatrix m_rho;             //matrix containing the estimate of the operator for doing 1-step ahead prediction (m x m)
  double m_p_threshold;                       //threshold according to how much predictive power has to be retained by the PPCs
  double m_alpha;                             //regularization parameter
  int m_k;                                    //number of PPCs retained
  bool m_k_imposed;                           //if true: it means that k has been passed as parameter, it is imposed from outside
                                              //if false: k has to be found (or from p or from eigvalus of phi) DA TENERE PERCHE' MI SERVE PER FARE CV

  std::vector<double> m_valid_err;        //just for debugging
  
public:
  
  PPC_KO_base(KO_Traits::StoringMatrix&& X,
              double threshold_ppc)
    :   
    m_X{std::forward<KO_Traits::StoringMatrix>(X)},
    m_m(X.rows()),
    m_n(X.cols()),
    m_p_threshold(threshold_ppc)
  
    {  
      
      //evaluating row mean and saving it in the m_means
      m_means = (m_X.rowwise().sum())/m_n;
      
      //normalizing
      for (size_t i = 0; i < m_n; ++i)
      {
        m_X.col(i) = m_X.col(i).array() - m_means;
      }
      
      // X * X' 
      m_Cov =  ((m_X*m_X.transpose()).array())/static_cast<double>(m_n);
      
      // trace of covariance
      m_trace_cov = m_Cov.trace();
      
      // (X[,2:n]*(X[,1:(n-1)])')/(n-1)
      m_CrossCov =  ((m_X.rightCols(m_n-1)*m_X.leftCols(m_n-1).transpose()).array())/(static_cast<double>(m_n-1));
      
      //squared of cross covariance
      m_GammaSquared = m_CrossCov.transpose()*m_CrossCov;
    }
  
  /*!
   * Virtual destructor (needed for inheritance)
   */
  virtual ~PPC_KO_base() = default;
  
  /*!
   * Virtual method to do prediction: will be defined in the children classes
   */
  virtual void solve() = 0; 
  
  
  
  /*!
   * Getter for m_valid_err
   */
  inline auto ValidErr() const {return m_valid_err;};
  
  /*!
   * Setter for m_CovReg (needed for CV)
   */
  inline auto & ValidErr() {return m_valid_err;};
  
  
  
  
  
  /*!
   * Getter for m_m
   */
  inline std::size_t m() const {return m_m;};
  
  /*!
   * Getter for m_n
   */
  inline std::size_t n() const {return m_n;};
  
  /*!
   * Getter for m_X
   */
  inline KO_Traits::StoringMatrix X() const {return m_X;};
  
  /*!
   * Getter for m_means
   */
  inline KO_Traits::StoringArray means() const {return m_means;};
  
  /*!
   * Getter for m_Cov
   */
  inline KO_Traits::StoringMatrix Cov() const {return m_Cov;};
  
  /*!
   * Getter for m_trace_cov
   */
  inline double trace_cov() const {return m_trace_cov;};
  
  /*!
   * Getter for m_CrossCov
   */
  inline KO_Traits::StoringMatrix CrossCov() const {return m_CrossCov;};
  
  /*!
   * Getter for m_CovReg
   */
  inline KO_Traits::StoringMatrix CovReg() const {return m_CovReg;};
  
  /*!
   * Setter for m_CovReg (needed for CV)
   */
  inline KO_Traits::StoringMatrix & CovReg() {return m_CovReg;};
  
  /*!
   * Getter for m_GammaSquared
   */
  inline KO_Traits::StoringMatrix GammaSquared() const {return m_GammaSquared;};
  
  /*!
   * Setter for m_GammaSquared (needed for CV)
   */
  inline KO_Traits::StoringMatrix & GammaSquared() {return m_GammaSquared;};
  
  /*!
   * Getter for m_CovRegRoot
   */
  inline KO_Traits::StoringMatrix CovRegRoot() const {return m_CovRegRoot;};
  
  /*!
   * Setter for m_CovRegRoot
   */
  inline KO_Traits::StoringMatrix & CovRegRoot() {return m_CovRegRoot;};
  
  /*!
   * Getter for m_a
   */
  inline KO_Traits::StoringMatrix a() const {return m_a;};
  
  /*!
   * Setter for m_a
   */
  inline KO_Traits::StoringMatrix & a() {return m_a;};
  
  /*!
   * Getter for m_b
   */
  inline KO_Traits::StoringMatrix b() const {return m_b;};
  
  /*!
   * Setter for m_b
   */
  inline KO_Traits::StoringMatrix & b() {return m_b;};
  
  /*!
   * Getter for m_rho
   */
  inline KO_Traits::StoringMatrix rho() const {return m_rho;};
  
  /*!
   * Setter for m_rho
   */
  inline KO_Traits::StoringMatrix &rho() {return m_rho;};
  
  /*!
   * Getter for m_p_threshold
   */
  inline double p_threshold() const {return m_p_threshold;};
  
  /*!
   * Getter for m_alpha
   */
  inline double alpha() const {return m_alpha;};
  
  /*!
   * Setter for m_alpha
   */
  inline double & alpha() {return m_alpha;};
  
  /*!
   * Getter for m_k
   */
  inline int k() const {return m_k;};
  
  /*!
   * Setter for m_k
   */
  inline int & k() {return m_k;};
  
  /*!
   * Getter for m_k_imposed
   */
  inline bool k_imposed() const {return m_k_imposed;};
  
  /*!
   * Setter for m_k_imposed
   */
  inline bool & k_imposed() {return m_k_imposed;};
  

  //number of PPCs retained
  int PPC_retained(const KO_Traits::StoringArray & eigvals) const;
  
  //Inverse square for regularized covariance (k is chosen in this function)
  KO_Traits::StoringMatrix matrix_inverse_root(const KO_Traits::StoringMatrix& gamma_alpha);
  
  //operator Phi estimate
  KO_Traits::StoringMatrix phi_estimate() const;
  
  //KO algorithm: called by the method solve(); if no CV: it is simply called once. If there is CV, it is called for each value of the parameters
  void KO_algo();
  
  //doing one-step ahead prediction
  KO_Traits::StoringArray prediction() const;
  
  //scores
  std::vector<double> scores() const;
};



//version without CV
class KO_NO_CV final : public PPC_KO_base
{
  
public:
  
  KO_NO_CV(KO_Traits::StoringMatrix&& X,
           double threshold_ppc,
           double alpha,
           int k)
    :   
    PPC::PPC_KO_base(std::move(X), threshold_ppc)
  
    { 
      this->alpha() = alpha;
      this->k() = k;
      this->k()==0 ? this->k_imposed() = false : this->k_imposed() = true;
      //for the case without CV, covariance regularized will be evaluated only once since regularization parameter has been passed as input parameter
      this->CovReg() = this->Cov().array() + this->alpha()*this->trace_cov()*(KO_Traits::StoringMatrix::Identity(this->m(),this->m()).array());
    }
  
  //virtual ~KO_NO_CV(){};
  
  void solve() override;              //overrided method to solve without CV
};





//function neede do a single run of CV for one parameter, fixing the other one
inline
KO_Traits::StoringVector 
ko_single_cv(KO_Traits::StoringMatrix && training_data,
             double threshold_ppc, 
             double alpha,
             int k)
{
  PPC::KO_NO_CV iter(std::move(training_data),threshold_ppc,alpha,k);
  iter.solve();
 
  return iter.prediction();
};
 



//version with CV for alpha
class KO_CV_alpha final : public PPC_KO_base
{
  
private:
  KO_Traits::StoringMatrix m_X_non_norm;          //data non normalized, necessary since is necessary to pass them at every cv iteration
  std::function<KO_Traits::StoringVector(KO_Traits::StoringMatrix,double,double,int)> m_cv_iter_f = ko_single_cv;
  std::function<double(KO_Traits::StoringVector)> m_ef = EF_PPC::mse<double>;
  
public:
  
  KO_CV_alpha(KO_Traits::StoringMatrix&& X,
              double threshold_ppc,
              int k)
    :   
    PPC::PPC_KO_base(std::move(X), threshold_ppc),
    m_X_non_norm(this->m(),this->n())
  
    {
      for (size_t i = 0; i < this->n(); ++i)
      {
        m_X_non_norm.col(i) = this->X().col(i).array() + this->means();
      }
      
      this->k() = k;            //k può essere passato, ma se non dovesse essere passato, significa che per l'alpha, guardo la explanatory power, non che ci faccia CV
      this->k()==0 ? this->k_imposed() = false : this->k_imposed() = true;
    }
  
  //virtual ~KO_CV_alpha(){};
  
  /*!
   * Getter for m_X_non_norm
   */
  inline KO_Traits::StoringMatrix X_non_norm() const {return m_X_non_norm;};
  
  /*!
   * Getter for m_cv_iter_f
   */
  inline std::function<KO_Traits::StoringVector(KO_Traits::StoringMatrix,double,double,int)> cv_iter_f() const {return m_cv_iter_f;};
  
  /*!
   * Getter for m_ef
   */
  inline std::function<double(KO_Traits::StoringVector)> ef() const {return m_ef;};
  
  
  double alpha_best_CV();             //to obtain the best alpha parameter for regularization
  void solve() override;              //overrided method to solve with CV
  
};



//version with CV for k
class KO_CV_k final : public PPC_KO_base
{
  
private:
  KO_Traits::StoringMatrix m_X_non_norm;          //data non normalized, necessary since is necessary to pass them at every cv iteration
  std::function<KO_Traits::StoringVector(KO_Traits::StoringMatrix,double,double,int)> m_cv_iter_f = ko_single_cv;
  std::function<double(KO_Traits::StoringVector)> m_ef = EF_PPC::mse<double>;
  
public:
  
  KO_CV_k(KO_Traits::StoringMatrix&& X,
          double threshold_ppc,
          double alpha)
    :   
    PPC::PPC_KO_base(std::move(X), threshold_ppc),
    m_X_non_norm(this->m(),this->n())
  
    {
      this->alpha() = alpha;
      this->CovReg() = this->Cov().array() + this->alpha()*this->trace_cov()*(KO_Traits::StoringMatrix::Identity(this->m(),this->m()).array());
      
      for (size_t i = 0; i < this->n(); ++i)
      {
        m_X_non_norm.col(i) = this->X().col(i).array() + this->means();
      }
    }
  
  //virtual ~KO_CV_k(){};
  
  /*!
   * Getter for m_X_non_norm
   */
  inline KO_Traits::StoringMatrix X_non_norm() const {return m_X_non_norm;};
  
  /*!
   * Getter for m_cv_iter_f
   */
  inline std::function<KO_Traits::StoringVector(KO_Traits::StoringMatrix,double,double,int)> cv_iter_f() const {return m_cv_iter_f;};
  
  /*!
   * Getter for m_ef
   */
  inline std::function<double(KO_Traits::StoringVector)> ef() const {return m_ef;};
  

  int k_best_CV();                    //to obtain the best number k of PPCs
  void solve() override;              //overrided method to solve with CV
  
};







//function needed in the cv for pair alpha-k: for every alpha, the CV on k begins: gives back the best pair
inline
std::pair<int,double>
ko_single_cv_on_k(KO_Traits::StoringMatrix && training_data,
                  double threshold_ppc, 
                  double alpha)
{
  //given an alpha, do the CV on k
  PPC::KO_CV_k iter(std::move(training_data),threshold_ppc,alpha);
  iter.solve();
  
  //find the best k associated to that alpha and its validation error
  int k_best = iter.k();
  double err_best = *(std::min_element(iter.ValidErr().begin(),iter.ValidErr().end()));
 
  return std::make_pair(k_best,err_best);
};






//versione with CV for both parameters
class KO_CV final : public PPC_KO_base
{
  KO_Traits::StoringMatrix m_X_non_norm;          //data non normalized, necessary since is necessary to pass them at every cv iteration
  std::function<std::pair<int,double>(KO_Traits::StoringMatrix,double,double)> m_cv_iter_on_k_f = ko_single_cv_on_k;

public:
  
  KO_CV(KO_Traits::StoringMatrix&& X,
        double threshold_ppc)
    :   
    PPC::PPC_KO_base(std::move(X), threshold_ppc),
    m_X_non_norm(this->m(),this->n())
  
    {
      for (size_t i = 0; i < this->n(); ++i)
      {
        m_X_non_norm.col(i) = this->X().col(i).array() + this->means();
      }
    }
  
  /*!
   * Getter for m_X_non_norm
   */
  inline KO_Traits::StoringMatrix X_non_norm() const {return m_X_non_norm;};
  

  /*!
   * Getter for m_cv_iter_on_k_f
   */
  inline std::function<std::pair<int,double>(KO_Traits::StoringMatrix,double,double)> cv_iter_on_k_f() const {return m_cv_iter_on_k_f;};
  

  std::pair<double,int> params_best_CV();             //to obtain the best parameters (alpha for regularization and k for number of PPCs)
  void solve() override;                              //overrided method to solve with CV
};






//test for CV CRTP alpha
class KO_CV_CRTP_ALPHA : public  PPC_KO_base
{
private:
  KO_Traits::StoringMatrix m_X_non_cent; 
  
public:
  KO_CV_CRTP_ALPHA(KO_Traits::StoringMatrix&& X, double threshold_ppc) 
    : 
    PPC::PPC_KO_base(std::move(X), threshold_ppc),
    m_X_non_cent(this->m(),this->n())
    {
      for (size_t i = 0; i < this->n(); ++i)
      {
        m_X_non_cent.col(i) = this->X().col(i).array() + this->means();
      }
    }
  
  
  /*!
   * Getter for m_X_non_cent
   */
  inline KO_Traits::StoringMatrix X_non_cent() const {return m_X_non_cent;};
  
  void solve() override;
};










//test for CV CRTP k
class KO_CV_CRTP_K : public  PPC_KO_base
{
private:
  KO_Traits::StoringMatrix m_X_non_cent; 
  
public:
  KO_CV_CRTP_K(KO_Traits::StoringMatrix&& X, double threshold_ppc, double alpha) 
    : 
    PPC::PPC_KO_base(std::move(X), threshold_ppc),
    m_X_non_cent(this->m(),this->n())
    {
      this->alpha() = alpha;
      this->CovReg() = this->Cov().array() + this->alpha()*this->trace_cov()*(KO_Traits::StoringMatrix::Identity(this->m(),this->m()).array());
      
      for (size_t i = 0; i < this->n(); ++i)
      {
        m_X_non_cent.col(i) = this->X().col(i).array() + this->means();
      }
    }
  
  
  /*!
   * Getter for m_X_non_cent
   */
  inline KO_Traits::StoringMatrix X_non_cent() const {return m_X_non_cent;};
  
  void solve() override;
};









//test for CV CRTP alpha-k
class KO_CV_CRTP_ALPHA_K : public  PPC_KO_base
{
private:
  KO_Traits::StoringMatrix m_X_non_cent; 
  
public:
  KO_CV_CRTP_ALPHA_K(KO_Traits::StoringMatrix&& X, double threshold_ppc) 
    : 
    PPC::PPC_KO_base(std::move(X), threshold_ppc),
    m_X_non_cent(this->m(),this->n())
    {
      for (size_t i = 0; i < this->n(); ++i)
      {
        m_X_non_cent.col(i) = this->X().col(i).array() + this->means();
      }
    }
  
  
  /*!
   * Getter for m_X_non_cent
   */
  inline KO_Traits::StoringMatrix X_non_cent() const {return m_X_non_cent;};
  
  void solve() override;
};


}//end namespace PPC

#endif /*KO_PPC_HPP*/