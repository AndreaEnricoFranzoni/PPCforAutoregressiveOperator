#ifndef KO_PPC_CRTP_HPP
#define KO_PPC_CRTP_HPP

#include <iostream>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <execution>
#include <vector>
#include <functional>
#include <utility>
#include <tuple>
#include <cmath>
#include <array>

#include "traits_ko.hpp"
#include "CV_include.hpp"
#include "Factory_cv_strategy.hpp"
#include "strategy_cv.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif


//templates params
//derived class (CRTP), domain dimension, if k is passed as param, if validation errors have to be stored and returned
template< class D, DOM_DIM dom_dim, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval > 
class PPC_KO_base
{
  
private:
  
  std::size_t m_m;                            //number of evaluation of the functional object (m)
  std::size_t m_n;                            //number of time instants (n)
  KO_Traits::StoringMatrix m_X;               //matrix storing time series (m x n), data centered
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
  std::vector<double> m_explanatory_power;    //vector containing the cumulative explanatory power (will have size k)
  double m_alpha;                             //regularization parameter
  int m_k;                                    //number of PPCs retained
  double m_threshold_ppc;                     //threshold according to how much predictive power has to be retained by the PPCs
  valid_err_variant m_valid_err;              //validation errors
  int m_number_threads;                       //number of threads for OMP
  
  
  
public:
  
  template<typename STOR_OBJ>
  PPC_KO_base(STOR_OBJ&& X,int number_threads)
    :   
    m_X{std::forward<STOR_OBJ>(X)},
    m_m(X.rows()),
    m_n(X.cols()),
    m_number_threads(number_threads)
    {  
      
      //evaluating row mean and saving it in the m_means
      m_means = (m_X.rowwise().sum())/m_n;
      
      //normalizing
#ifdef _OPENMP
#pragma omp parallel for num_threads(m_number_threads)
#endif
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
  
  
  
  //Getters and setters
  
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
   * Setter for m_CovReg (needed for CV)
   */
  inline KO_Traits::StoringMatrix & CovReg() {return m_CovReg;};
  
  /*!
   * Getter for m_rho
   */
  inline KO_Traits::StoringMatrix rho() const {return m_rho;};
  
  /*!
   * Getter for m_a
   */
  inline KO_Traits::StoringMatrix a() const {return m_a;};
  
  /*!
   * Getter for m_b
   */
  inline KO_Traits::StoringMatrix b() const {return m_b;};
  
  /*!
   * Getter for m_explanatory_power
   */
  inline std::vector<double> explanatory_power() const {return m_explanatory_power;};
  
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
   * Getter for m_threshold_ppc
   */
  inline double threshold_ppc() const {return m_threshold_ppc;};
  
  /*!
   * Setter for m_threshold_ppc
   */
  inline double & threshold_ppc() {return m_threshold_ppc;};
  
  /*!
   * Getter for m_valid_err
   */
  inline valid_err_variant ValidErr() const {return m_valid_err;};
  
  /*!
   * Setter for m_valid_err
   */
  inline valid_err_variant & ValidErr() {return m_valid_err;};
  
  /*!
   * Getter for m_number_threads
   */
  inline int number_threads() const {return m_number_threads;};
  

  //methods common to all child classes
  //number of PPCs retained, eigvalues and eigvcts di phi
  std::tuple<int,KO_Traits::StoringVector,KO_Traits::StoringMatrix> PPC_retained(const KO_Traits::StoringMatrix & phi, double tot_exp_pow) const;

  //KO algorithm, once parameters have been set
  void KO_algo();
  
  //one-step ahead prediction
  KO_Traits::StoringArray prediction() const;
  
  //scores evaluation
  std::vector<double> scores() const;
  
  //sd scores dir and wei
  std::vector<std::array<double,2>> sd_scores_dir_wei() const;
  
  //solving
  inline
  void 
  solve() 
  {
    static_cast<D*>(this)->solving();   //solving depends on child class: downcasting of CRTP
  }
};
  

#include "PPC_KO_imp.hpp"

#endif  //KO_PPC_CRTP_HPP