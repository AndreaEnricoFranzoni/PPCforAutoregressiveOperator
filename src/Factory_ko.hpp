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

#ifndef KO_FACTORY_HPP
#define KO_FACTORY_HPP

#include <string>
#include <memory>
#include <stdexcept>
#include <utility>
#include <iostream>

#include "traits_ko.hpp"
#include "parameters_wrapper.hpp"
#include "PPC_KO_wrapper.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif


/*!
* @file Factory_ko.hpp
* @brief Factory for generating the PPCKO solver
* @author Andrea Enrico Franzoni
*/


/*!
* @class KO_Factory
* @brief Generating the PPCKO solver runtime according to an input string
* @tparam solver if algorithm solved inverting the regularized covariance or avoiding it through gep (not possible if retaining the number of PPCs with explanatory power criterion)
* @tparam k_imp if k is imposed or has to be found through explanatory power criterion
* @tparam valid_err_ret if validation error are stored
* @tparam cv_strat strategy for splitting training/validation sets
* @tparam err_eval how to evaluate the loss between prediction on validation set and validation set
*/
template< SOLVER solver, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval > 
class KO_Factory
{	
public:
  //! Static method that takes a string as identifier and builds a pointer to the right object for the cross-validation requested
  /*!
  * @brief Generating the PPCKO solver runtime according to an input string. It raises an error if a wrong one is passed
  * @param id input string:
  *           - 'NoCV': no cross-validation is performed. Regularization is imposed. Number of PPCs can be imposed or retrieved through explanatory power criterion;
  *           - 'CV_alpha': cross-validation for regularization parameter. Number of PPCs can be imposed or retrieved through explanatory power criterion;
  *           - 'CV_k': cross-validation for number of retained PPCs. Regularization parameter is imposed;
  *           - 'CV': cross-validation for both regularization parameter and number of retained PPCs.
  * @param X matrix containing the fts
  * @param alpha regularization parameter
  * @param k number of retained PPCs:
  *          - 'k' = 0: if possible (no cv on ok), number of retained PPCs is selected through explanatory power criterion;
  *          - 'k' > 0: if possible (no cv on ok), number of retained PPCs is imposed.
  * @param threshold_ppc requested explanatory power from the PPCs. Used only for selecting 'k' through explanatory power criterion
  * @param alphas input space for regularization parameter
  * @param k_s input space for the number of retained PPCs
  * @param toll the cv on the number of retained PPCs continues only if between two parameters, that are checked in increasing order, the absolute difference between two validation errors is bigger than tolerance*trace(covariance). If not, stops and look for k only between the tested ones 
  * @param min_size_ts smallest training set size (number of time instants)
  * @param max_size_ts biggest training set size (number of time instants)
  * @param num_threads number of threads for OMP
  * @return a unique pointer to a PPC_KO_wrapper<solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>
  * @details the method is static. std::unique_ptr provides a safer architecture for factory purposes
  */
  static 
  std::unique_ptr<PPC_KO_wrapper<solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>> 
    KO_solver(const std::string &id,
              KO_Traits::StoringMatrix && X,
              double alpha,
              int k,
              double threshold_ppc,
              const std::vector<double>& alphas,
              const std::vector<int>& k_s,
              double toll,
              int min_size_ts,
              int max_size_ts,
              int num_threads)
    {
      
#ifdef _OPENMP
      if(num_threads==1)
      {
        std::cout << "Running parallel version with " << num_threads << " thread" << std::endl;
      }
      else
      {
        std::cout << "Running parallel version with " << num_threads << " threads" << std::endl;
      }
#else
      std::cout << "Running serial version" << std::endl;
#endif
      
      if (id == CV_algo::CV1)   //No CV:          if (k=0): explanatory power criterion
      {
        return k==0 ? std::make_unique<PPC_KO_wrapper_no_cv<solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>>(std::move(X),alpha,threshold_ppc,num_threads) : std::make_unique<PPC_KO_wrapper_no_cv<solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>>(std::move(X),alpha,k,num_threads);
      }
      
      if (id == CV_algo::CV2)   //CV on alpha:    if (k=0): explanatory power criterion
      {
        return k==0 ? std::make_unique<PPC_KO_wrapper_cv_alpha<solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>>(std::move(X),alphas,threshold_ppc,min_size_ts,max_size_ts,num_threads) : std::make_unique<PPC_KO_wrapper_cv_alpha<solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>>(std::move(X),alphas,k,min_size_ts,max_size_ts,num_threads);
      }
      
      if (id == CV_algo::CV3)   //CV on k
      {
        return std::make_unique<PPC_KO_wrapper_cv_k<solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>>(std::move(X),alpha,k_s,toll,min_size_ts,max_size_ts,num_threads);
      }
      
      if (id == CV_algo::CV4)   //CV on both alpha and k
      {
        return std::make_unique<PPC_KO_wrapper_cv_alpha_k<solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>>(std::move(X),alphas,k_s,toll,min_size_ts,max_size_ts,num_threads);
      }
      
      else
      {
        std::string error_message = "Wrong input string";
        throw std::invalid_argument(error_message);
      }
    }
};

#endif //KO_FACTORY_HPP