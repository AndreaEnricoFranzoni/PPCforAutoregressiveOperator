#ifndef KO_PPC_HPP
#define KO_PPC_HPP

#include "KO_Traits.hpp"
#include <iostream>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <execution>

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
  
  std::size_t m_m;                            //number of temporal series (m)
  std::size_t m_n;                            //number of time instants (n)
  KO_Traits::StoringMatrix m_X;               //matrix storing time series (m x n)
  KO_Traits::StoringArray m_means;            //vector storing the time series means (m x 1)  
  KO_Traits::StoringMatrix m_Cov;             //matrix estimating the covariance operator (m x m)
  KO_Traits::StoringMatrix m_CrossCov;        //matrix estimating the cross-covariance operator (m x m)
  KO_Traits::StoringMatrix m_CovReg;          //matrix containing the covariance + alpha*I (m x m) (used only for the non-CV version)
  KO_Traits::StoringMatrix m_GammaSquared;    //matrix containing the squared of the cross-covariance matrix (m x m) (used only for the non-CV version)
  KO_Traits::StoringMatrix m_CovRegRoot;      //matrix containing the inverse squared root of regularized covariance (m x m)
  KO_Traits::StoringMatrix m_a;               //matrix containing predictive loadings (each col)  (m x k)
  KO_Traits::StoringMatrix m_b;               //matrix containing predictive factors (each col)   (m x k)        
  KO_Traits::StoringMatrix m_rho;             //matrix containing the estimate of the operator for doing 1-step ahead prediction (m x m)
  double m_p_threshold = 0.95;                //threshold according to how much predictive power has to be reatined by the PPCs
  double m_alpha = 0.75;                      //regularization parameter  
  int m_k;                                    //number of PPCs retained                  
  
public:
  
  PPC_KO_base(KO_Traits::StoringMatrix&& X)
    :   
    m_X{std::forward<KO_Traits::StoringMatrix>(X)},
    m_m(X.rows()),
    m_n(X.cols())  
  
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
      
      // (X[,2:n]*(X[,1:(n-1)])')/(n-1)
      m_CrossCov =  ((m_X.rightCols(m_n-1)*m_X.leftCols(m_n-1).transpose()).array())/(static_cast<double>(m_n-1));      
      
    }
  
  /*!
   * Virtual destructor
   */
  virtual ~PPC_KO_base() = default;
  
  /*!
   * Virtual method to do prediction: will be defined in the children classes
   */
  virtual void solve() = 0; 
  
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
  
  //number of PPCs retained
  int PPC_retained(const KO_Traits::StoringArray & cov_reg_eigvals);
  
  //Inverse square for reg covariance
  KO_Traits::StoringMatrix matrix_inverse_root(const KO_Traits::StoringMatrix& gamma_alpha);
  
  //operator Phi estimate
  KO_Traits::StoringMatrix const phi_estimate(const KO_Traits::StoringMatrix& gamma_alpha_rooted);
  
  //KO algorithm: called by the method solve(); if no CV: it is simply called once. If there is CV, it is called for each value of alpha
  void KO_algo();
  
  //doing one-step ahead prediction
  KO_Traits::StoringArray prediction() const;
};



//version without CV
class KO_NO_CV final : public PPC_KO_base
{
  
public:
  
  KO_NO_CV(KO_Traits::StoringMatrix&& X)
    :   
    PPC::PPC_KO_base(std::move(X))
    {   
      this->CovReg() = this->Cov().array() + this->alpha()*(KO_Traits::StoringMatrix::Identity(this->m(),this->m()).array()); 
      this->GammaSquared() = (this->CrossCov().transpose())*this->CrossCov();
    }
  //virtual ~KO_NO_CV(){};
  
  
  void solve() override;
  
  
};



//version with CV
class KO_CV final : public PPC_KO_base
{
private:
  static constexpr double alpha_min = 1.7e-308;
  static constexpr double alpha_max = 1.5;
  std::size_t m_n_disc = 100;
  
public:
  
  KO_CV(KO_Traits::StoringMatrix&& X)
    :   
    PPC::PPC_KO_base(std::move(X))
    {}
  
  void CV();
  void solve() override;
  
  /*!
   * Getter for m_n_disc
   */
  inline auto n_disc() const {return m_n_disc;};
  
  /*!
   * Setter for m_n_disc
   */
  inline auto & n_disc() {return m_n_disc;};
};

}           //end namespace PPC


#endif /*KO_PPC_HPP*/