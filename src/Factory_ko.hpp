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



template< SOLVER solver, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval > 
class KO_Factory
{	
public:
  //! Static method that takes a string as identifier and builds a pointer to the right object for the cross-validation requested
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
      
      if (id == CV_algo::CV1)   //No CV
      {
        return k==0 ? std::make_unique<PPC_KO_wrapper_no_cv<solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>>(std::move(X),alpha,threshold_ppc,num_threads) : std::make_unique<PPC_KO_wrapper_no_cv<solver,k_imp,valid_err_ret,cv_strat,cv_err_eval>>(std::move(X),alpha,k,num_threads);
      }
      
      if (id == CV_algo::CV2)   //CV on alpha
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
