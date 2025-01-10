// Copyright (c) 2024 Andrea Enrico Franzoni (andreaenrico.franzoni@gmail.com)
//
// This file is part of PPCKO
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of PPCKO and associated documentation files (the PPCKO software), to deal
// PPCKO without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of PPCKO, and to permit persons to whom PPCKO is
// furnished to do so, subject to the following conditions:
//
// PPCKO IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH PPCKO OR THE USE OR OTHER DEALINGS IN
// PPCKO.

#ifndef KO_PPC_CRTP_HPP
#define KO_PPC_CRTP_HPP


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


/*!
* @file PPC_KO.hpp
* @brief Class for computing PPCKO algortihm. Hierarchy of classes, for different algorithm versions. Static polymorphism
* @author Andrea Enrico Franzoni
*/



/*!
* @class PPC_KO_base
* @brief Base class for computing PPCKO algorithm computations. Polymorphism is known at compile-time (CRTP)
* @tparam D type of the derived class (for static polymorphism thorugh CRTP):
*         - 'PPC_KO_NoCV': alpha has to be passed as paramter. k can be passed as paramter or selected through explanatory power criterion
*         -
*         -
*         - 
* @tparam solver if algorithm solved inverting the regularized covariance or avoiding it through gep (not possible if retaining the number of PPCs with explanatory power criterion)
* @tparam k_imp if k is imposed or has to be found through explanatory power criterion
* @tparam valid_err_ret if validation error are stored
* @tparam cv_strat strategy for splitting training/validation sets
* @tparam err_eval how to evaluate the loss between prediction on validation set and validation set
* @details It is the base class. Polymorphism is known at compile time thanks to Curiously Recursive Template Pattern (CRTP) 
*/
template< class D, SOLVER solver, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval > 
class PPC_KO_base
{
  
private:
   
  /*!Number of evaluations, for each instant, of the curve/surface (number of rows of data matrix) (m)*/ 
  std::size_t m_m;
  /*!Number of time instants of the fts (number of columns of data matrix) (n)*/                            
  std::size_t m_n;                            
  /*!Fts: data will be centered as soon as object construction (matrix: m x n)*/
  KO_Traits::StoringMatrix m_X;               
  /*!Fts mean function (array: m x 1)*/
  KO_Traits::StoringArray m_means;            
  /*!Covariance operator estimate (matrix: m x m)*/
  KO_Traits::StoringMatrix m_Cov;             
  /*!Trace of the covariance operator estimate*/
  double m_trace_cov;                         
  /*!Cross-covariance operator estimate (matrix: m x m)*/
  KO_Traits::StoringMatrix m_CrossCov;        
  /*!Regularized sample covariance (sample covariance + alpha*trace(cov)*I) (matrix: m x m)*/
  KO_Traits::StoringMatrix m_CovReg;          
  /*!Square of the cross-covariance operator estimate (matrix: m x m)*/
  KO_Traits::StoringMatrix m_GammaSquared;    
  /*!Inverse square root of the regularized covariance (matrix: m x m)*/
  KO_Traits::StoringMatrix m_CovRegRoot;      
  /*!Predictive loading (PPCs directions) (matrix: m x k)*/
  KO_Traits::StoringMatrix m_a;           
  /*!Predictive factors factor (PPCs weights) (matrix: m x k)*/
  KO_Traits::StoringMatrix m_b;               
  /*!Estimate of the autoregressive operator for doing one-step ahead prediction (matrix: m x m)*/      
  KO_Traits::StoringMatrix m_rho;             
  /*!Cumulative explanatory power of PPCs (vector of size k)*/
  std::vector<double> m_explanatory_power;    
  /*!Total explanatory power intrinsic in the fts*/
  double m_tot_exp_pow;               
  /*!Regularization parameter*/
  double m_alpha;                             
  /*!Number of retained PPCs*/
  int m_k;                                   
  /*!Requested explanatory power from the retained PPCs*/
  double m_threshold_ppc;                     
  /*!Validation errors*/
  valid_err_variant m_valid_err;             
  /*!Number of threads for OMP*/
  int m_number_threads;                      
  
  
public:
  
  /*!
  * @brief Constructor: centers data, evaluate mean function, sample covariance, sample cross-covariance and its square
  * @param X fts
  * @param number_threads number of threads for OMP
  * @details Universal constructor: move semantic used to optimazing handling big size objects
  * @note eventual usage of 'pragma' directive for OMP
  */
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
      
      //centering
#ifdef _OPENMP
#pragma omp parallel for num_threads(m_number_threads)
#endif
      for (size_t i = 0; i < m_n; ++i)
      {
        m_X.col(i) = m_X.col(i).array() - m_means;
      }
      
      // covariance operator estimate: (X * X')/n
      m_Cov =  ((m_X*m_X.transpose()).array())/static_cast<double>(m_n);
      
      // trace of covariance
      m_trace_cov = m_Cov.trace();
      
      // cross-covariance operator estimate: (X[,2:n]*(X[,1:(n-1)])')/(n-1)
      m_CrossCov =  ((m_X.rightCols(m_n-1)*m_X.leftCols(m_n-1).transpose()).array())/(static_cast<double>(m_n-1));
      
      // square of cross covariance estimate
      m_GammaSquared = m_CrossCov.transpose()*m_CrossCov;
    }
  
  
  /*!
  * @brief Getter for the number of evaluation of the curve/surface
  * @return the private m_m
  */
  inline std::size_t m() const {return m_m;};
  
  /*!
  * @brief Getter for the number of time instants
  * @return the private m_n
  */
  inline std::size_t n() const {return m_n;};
  
  /*!
  * @brief Getter for fts data matrix (centered)
  * @return the private m_X
  */
  inline KO_Traits::StoringMatrix X() const {return m_X;};
  
  /*!
  * @brief Getter for the mean function
  * @return the private m_means
  */
  inline KO_Traits::StoringArray means() const {return m_means;};
  
  /*!
  * @brief Getter for the covariance operator estimate
  * @return the private m_Cov
  */
  inline KO_Traits::StoringMatrix Cov() const {return m_Cov;};
  
  /*!
  * @brief Getter for the covariance operator estimate's trace
  * @return the private m_trace_cov
  */
  inline double trace_cov() const {return m_trace_cov;};
  
  /*!
  * @brief Setter for the regularized sample covariance operator
  * @return the private m_CovReg (non-const)
  */
  inline KO_Traits::StoringMatrix & CovReg() {return m_CovReg;};
  
  /*!
  * @brief Getter for the autoregressive operator estimate
  * @return the private m_rho
  */
  inline KO_Traits::StoringMatrix rho() const {return m_rho;};
  
  /*!
  * @brief Getter for the predictive loadings (PPCs directions)
  * @return the private m_a
  */
  inline KO_Traits::StoringMatrix a() const {return m_a;};
  
  /*!
  * @brief Getter for the predictive factors factor (PPCs weights)
  * @return the private m_b
  */
  inline KO_Traits::StoringMatrix b() const {return m_b;};
  
  /*!
  * @brief Getter for the cumulative explanatory power of the PPCs
  * @return the private m_explanatory_power
  */
  inline std::vector<double> explanatory_power() const {return m_explanatory_power;};
  
  /*!
  * @brief Getter for the regularization parameter
  * @return the private m_alpha
  */
  inline double alpha() const {return m_alpha;};
  
  /*!
  * @brief Setter for the regularization parameter
  * @return the private m_alpha (non-const)
  */
  inline double & alpha() {return m_alpha;};
  
  /*!
  * @brief Getter for the number of retained PPCs
  * @return the private m_k
  */
  inline int k() const {return m_k;};
  
  /*!
  * @brief Setter for the number of retained PPCs
  * @return the private m_k (non-const)
  */
  inline int & k() {return m_k;};
  
  /*!
  * @brief Getter for the requested explanatory power of the retained PPCs
  * @return the private m_threshold_ppc
  */
  inline double threshold_ppc() const {return m_threshold_ppc;};
  
  /*!
  * @brief Setter for the requested explanatory power of the retained PPCs
  * @return the private m_threshold_ppc (non-const)
  */
  inline double & threshold_ppc() {return m_threshold_ppc;};
  
  /*!
  * @brief Getter for the validation errors
  * @return the private m_valid_err
  */
  inline valid_err_variant ValidErr() const {return m_valid_err;};
  
  /*!
  * @brief Setter for the validation errors
  * @return the private m_valid_err (non-const)
  */
  inline valid_err_variant & ValidErr() {return m_valid_err;};
  
  /*!
  * @brief Getter for the number of threads for OMP
  * @return the private m_number_threads
  */
  inline int number_threads() const {return m_number_threads;};
  

  /*!
  * @brief Retaining the the PPCs: pairs eigenvalue-eigenvector and their number
  * @return a tuple containing: the number of retained PPCs, the eigenvalues of phi/of GEP, the eigenvectors of phi/of GEP
  * @details The PPCs are computed according to the solver strategy using 'Spectra'. Only the first k pairs eigenvalue/eigenvactor are evaluated,
  *          corresponding to the k laregest eigenvalues, if k imposed. If instead (only for 'SOLVER::ex_solver') are computed
  *          using the explanatory power criterion, the k pairs eigenvalues-eigenvectors are computed increasing the number of computed
  *          ones until the requested explanatory power is reached
  */
  std::tuple<int,KO_Traits::StoringVector,KO_Traits::StoringMatrix> PPC_retained();
  
  /*!
  * @brief Performing PPCKO algorithm once regularization parameter is selected and k or it is fixed or to be retained through explanatory power.
  *        Computes PPCs, direction and weight, their number and their cumulative explanatory power, and the estimate of the autoregressive operator
  * @details Modifying the private members of the class corresponding to the computed quantities
  */
  void KO_algo();
  
  /*!
  * @brief Performs one-step ahead prediction of the fts. The mean function is added
  * @return the array of the prediction
  */
  KO_Traits::StoringArray prediction() const;
  
  /*!
  * @brief Computes the scores of the PPCs, defined as scalar product between the direction and the fts at the last instant
  * @return a vector containing the score of each PPC
  */
  std::vector<double> scores() const;
  
  /*!
  * @brief Computes the standard deviation of the scores of directions and weights
  * @return a vector (of size equal to the number of PPCs) containing arrays with two elements (standard deviation of direction and weight score of the PPC)
  * @details - Scores of directions are computed as the scalar product of the direction and the fts in the instants between 2 and n. 
  *          - Scores of weights are computed as the scalar product of the weight and the fts in the instants between 1 and n-1     
  */
  std::vector<std::array<double,2>> sd_scores_dir_wei() const;
  
  /*!
  * @brief Method to solve PPCKO according to which cross-validation is performed, if any
  * @details Entails downcasting of base class with a static cast of pointer to the derived-known-at-compile-time class, CRTP fashion
  */
  inline
  void 
  solve() 
  {
    static_cast<D*>(this)->solving();   //solving depends on child class: downcasting with CRTP of base to derived
  }
};
  

#include "PPC_KO_imp.hpp"

#endif  //KO_PPC_CRTP_HPP