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

#include "PPC_KO.hpp"

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include "spectra/include/Spectra/MatOp/DenseSymMatProd.h"
#include "spectra/include/Spectra/MatOp/DenseCholesky.h"
#include "spectra/include/Spectra/SymEigsSolver.h"
#include "spectra/include/Spectra/SymGEigsSolver.h"


/*!
* @file PPC_KO_imp.hpp
* @brief Definition of methods of the base class for computing PPCKO
* @author Andrea Enrico Franzoni
*/



/*!
* @brief Retaining the the PPCs: pairs eigenvalue-eigenvector and their number
* @return a tuple containing: the number of retained PPCs, the eigenvalues of phi/of GEP, the eigenvectors of phi/of GEP
* @details The PPCs are computed according to the solver strategy using 'Spectra'. Only the first k pairs eigenvalue/eigenvactor are evaluated,
*          corresponding to the k laregest eigenvalues, if k imposed. If instead (only for 'SOLVER::ex_solver') are computed
*          using the explanatory power criterion, the k pairs eigenvalues-eigenvectors are computed increasing the number of computed
*          ones until the requested explanatory power is reached
*/
template< class D, SOLVER solver, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval >
std::tuple<int,KO_Traits::StoringVector,KO_Traits::StoringMatrix>
PPC_KO_base<D, solver, k_imp, valid_err_ret, cv_strat, cv_err_eval>::PPC_retained()
{

  //exact solver: can be used for k not imp (selected through explanatory power criterion) and k imp (by the user of by cv process)
  if constexpr(solver == SOLVER::ex_solver)
  {
    //Square root inverse of reg covariance: self-adjoint:exploiting it
    Eigen::SelfAdjointEigenSolver<KO_Traits::StoringMatrix> eigensolver_cov_reg(m_CovReg);
    m_CovRegRoot = eigensolver_cov_reg.operatorInverseSqrt();

    //Phi estimate: self-adjoint:exploiting it
    KO_Traits::StoringMatrix phi_hat = m_CovRegRoot*m_GammaSquared*m_CovRegRoot.transpose();
    //sum of phi eigenvalues: its trace
    m_tot_exp_pow = phi_hat.trace();

    //PPCS are found through Spectra, for efficiency
    Spectra::DenseSymMatProd<double> op(phi_hat);
    
    if constexpr( k_imp == K_IMP::NO )    //number of PPCs to be selected through explanatory power
    {
      //compute i pairs, with i staring from 1, increasing i until the requested explnatory power is reached
      for(std::size_t i = 0; i < m_m; ++i)
      {
        int n_ppcs = i+1;
        //Spectra framework
        Spectra::SymEigsSolver<Spectra::DenseSymMatProd<double>> eigsolver_phi(op, n_ppcs, 2*n_ppcs);
        eigsolver_phi.init();
        int nconv = eigsolver_phi.compute(Spectra::SortRule::LargestAlge);

        //if explanatory power reached: return
        if(eigsolver_phi.eigenvalues().sum()/m_tot_exp_pow >= m_threshold_ppc)
        {
          return std::make_tuple(n_ppcs,eigsolver_phi.eigenvalues(),eigsolver_phi.eigenvectors());
        }
      }
    } 
    else              // number of PPCs already known (imposed by the user of by cv process)
    {
      //Spectra framework
      Spectra::SymEigsSolver<Spectra::DenseSymMatProd<double>> eigsolver_phi(op, m_k, 2*m_k);
      eigsolver_phi.init();
      int nconv = eigsolver_phi.compute(Spectra::SortRule::LargestAlge);
      
      return std::make_tuple(m_k,eigsolver_phi.eigenvalues(),eigsolver_phi.eigenvectors());
    }
  }

  //gep solver: quicker, but only if you impose k by the user of by cv process
  else if constexpr(solver == SOLVER::gep_solver)
  {
    if constexpr(k_imp == K_IMP::YES)    
    {
      //preparing GEP: m_GammaSquared*v = lambda*m_CovReg*v, v geigvct, lambda geigval
      Spectra::DenseSymMatProd<double> op(m_GammaSquared);
      Spectra::DenseCholesky<double>  Bop(m_CovReg);    //since it is a covariance: sdp: Cholesky dec for efficiency

      //Spectra framework
      Spectra::SymGEigsSolver<Spectra::DenseSymMatProd<double>, Spectra::DenseCholesky<double>, Spectra::GEigsMode::Cholesky> eigsolver_ppc(op, Bop, m_k, 2*m_k);
      eigsolver_ppc.init();
      int nconv = eigsolver_ppc.compute(Spectra::SortRule::LargestAlge);
      //since the total sum of the eigenvalues is not for free: at least we can compare magnitude between the retained ones
      m_tot_exp_pow = eigsolver_ppc.eigenvalues().sum();
      
      return std::make_tuple(m_k,eigsolver_ppc.eigenvalues(),eigsolver_ppc.eigenvectors());
    }
  }
}





/*!
* @brief Performing PPCKO algorithm once regularization parameter is selected and k or it is fixed or to be retained through explanatory power.
*        Computes PPCs, direction and weight, their number and their cumulative explanatory power, and the estimate of the autoregressive operator
* @details Modifying the private members of the class corresponding to the computed quantities.
*/
template< class D, SOLVER solver, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval >
void
PPC_KO_base<D, solver, k_imp, valid_err_ret, cv_strat, cv_err_eval>::KO_algo()
{ 
  //finding the PPCs
  auto ppcs_ret = this->PPC_retained();
    
  //saving their number if retained 
  if constexpr( k_imp == K_IMP::NO)
  {
    m_k = std::get<0>(ppcs_ret);
  }
  
  //cumulative explanatory power: if ex_solver, is coherent. If gep_solver, is only the relative magnitude of the retained eigenvalues
  m_explanatory_power.resize(m_k);
  std::partial_sum(std::get<1>(ppcs_ret).begin(),std::get<1>(ppcs_ret).end(),m_explanatory_power.begin());        
  std::for_each(m_explanatory_power.begin(),m_explanatory_power.end(),[this](auto &el){el=el/m_tot_exp_pow;});
  
  //Weights (b_i): if gep, their for free. If not, cross-covariance has to be applie
  if constexpr(solver == SOLVER::ex_solver){m_b = m_CovRegRoot*std::get<2>(ppcs_ret);}  else{m_b = std::get<2>(ppcs_ret);}
  
  //Directions (a_i)
  m_a = m_CrossCov*m_b;

  //autoregressive operator estimate
  m_rho = m_a*(m_b.transpose());

}



/*!
* @brief Performs one-step ahead prediction of the fts. The mean function is added
* @return the array of the prediction
*/
template< class D, SOLVER solver, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval >
KO_Traits::StoringArray
PPC_KO_base<D, solver, k_imp, valid_err_ret, cv_strat, cv_err_eval>::prediction()
const 
{
  //Applying the estimated autoregressive operator and adding the mean function
  return (m_rho*m_X.col(m_n-1)).array() + m_means;
}



/*!
* @brief Computes the scores of the PPCs, defined as scalar product between the direction and the fts at the last instant
* @return a vector containing the score of each PPC
*/
template< class D, SOLVER solver, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval >
std::vector<double>
PPC_KO_base<D, solver, k_imp, valid_err_ret, cv_strat, cv_err_eval>::scores()
const
{ 
  std::vector<double> scores;
  scores.reserve(m_k);
  
  for(std::size_t i = 0; i < m_k; ++i)
  {
    scores.emplace_back((m_X.col(m_n-1)).dot((m_a.col(i))));
  }
  
  return scores;
}



/*!
* @brief Computes the standard deviation of the scores of directions and weights
* @return a vector (of size equal to the number of PPCs) containing arrays with two elements (standard deviation of direction and weight score of the PPC)
* @details - Scores of directions are computed as the scalar product of the direction and the fts in the instants between 2 and n. 
*          - Scores of weights are computed as the scalar product of the weight and the fts in the instants between 1 and n-1     
*/
template< class D, SOLVER solver, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval >
std::vector<std::array<double,2>>
PPC_KO_base<D, solver, k_imp, valid_err_ret, cv_strat, cv_err_eval>::sd_scores_dir_wei()
const
{
  std::vector<std::array<double,2>> standard_dev;
  standard_dev.reserve(m_k);
  
  //the number of computed scalar products is n - 1 for both directions and weights
  std::size_t n = m_X.cols() - 1;
  
  // for each one of the PPC
  for(std::size_t comp = 0; comp < m_k; ++comp){
    
    std::vector<double> scores_dir;
    scores_dir.reserve(n);
    std::vector<double> scores_wei;
    scores_wei.reserve(n);
    
    //computing the scores of directions (it has to be for the next value) and weights (current value)
    for(std::size_t i = 0; i < n; ++i){
      scores_dir.emplace_back(m_X.col(i+1).dot(m_a.col(comp)));
      scores_wei.emplace_back(m_X.col(i).dot(m_b.col(comp)));
    }
    
    //compute the standard deviation
    // mean
    double mean_dir = std::accumulate(scores_dir.begin(), scores_dir.end(), 0.0)/n;
    double mean_wei = std::accumulate(scores_wei.begin(), scores_wei.end(), 0.0)/n;
    
    // variance 
    double variance_dir = std::transform_reduce(scores_dir.begin(), 
                                                scores_dir.end(), 
                                                0.0,
                                                std::plus{},
                                                [mean_dir](auto el) {
                                                return std::pow(el-mean_dir,2);})/n;
    
    double variance_wei = std::transform_reduce(scores_wei.begin(), 
                                                scores_wei.end(), 
                                                0.0,
                                                std::plus{},
                                                [mean_wei](auto el) {
                                                return std::pow(el-mean_wei,2);})/n;
    
    // standard deviation
    standard_dev.emplace_back(std::array<double, 2>{std::sqrt(variance_dir),std::sqrt(variance_wei)});
    
    scores_dir.clear();
    scores_wei.clear();
  }
  
  return standard_dev;
}