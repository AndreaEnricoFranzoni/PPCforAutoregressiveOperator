#ifndef KO_FACTORY_HPP
#define KO_FACTORY_HPP

#include <string>
#include <memory>
#include <stdexcept>
#include "utility"

#include "traits_ko.hpp"
#include "PPC_KO_wrapper.hpp"
#include "mesh.hpp"




template< DOM_DIM dom_dim, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval > 
class KO_Factory
{	
public:
  //! Static method that takes a string as identifier and builds a pointer to the right object for the cross-validation requested
  static 
  std::unique_ptr<PPC_KO_wrapper<dom_dim,k_imp,valid_err_ret,cv_strat,cv_err_eval>> 
    KO_solver(const std::string &id,
              KO_Traits::StoringMatrix && X,
              Geometry::Mesh1D && grid_func_eval,
              double alpha,
              int k,
              double threshold_ppc,
              const std::vector<double>& alphas,
              const std::vector<int>& k_s,
              double toll)
    {
      
      if (id == "NoCV")
      {
        if(k==0){return std::make_unique<PPC_KO_wrapper_no_cv<dom_dim,k_imp,valid_err_ret,cv_strat,cv_err_eval>>(std::move(X),std::move(grid_func_eval),alpha,threshold_ppc);}
        else    {return std::make_unique<PPC_KO_wrapper_no_cv<dom_dim,k_imp,valid_err_ret,cv_strat,cv_err_eval>>(std::move(X),std::move(grid_func_eval),alpha,k);}
      }
      
      if (id == "CV_alpha")
      {
        if(k==0){return std::make_unique<PPC_KO_wrapper_cv_alpha<dom_dim,k_imp,valid_err_ret,cv_strat,cv_err_eval>>(std::move(X),std::move(grid_func_eval),alphas,threshold_ppc);}
        else    {return std::make_unique<PPC_KO_wrapper_cv_alpha<dom_dim,k_imp,valid_err_ret,cv_strat,cv_err_eval>>(std::move(X),std::move(grid_func_eval),alphas,k);}
      }
      
      if (id == "CV_k")
      {
        return std::make_unique<PPC_KO_wrapper_cv_k<dom_dim,k_imp,valid_err_ret,cv_strat,cv_err_eval>>(std::move(X),std::move(grid_func_eval),alpha,k_s,toll);
      }
      
      if (id == "CV")
      {
        return std::make_unique<PPC_KO_wrapper_cv_alpha_k<dom_dim,k_imp,valid_err_ret,cv_strat,cv_err_eval>>(std::move(X),std::move(grid_func_eval),alphas,k_s,toll);
      }
      
      else
      {
        std::string error_message = "Wrong input string";
        throw std::invalid_argument(error_message);
      }
    }
};

#endif //KO_FACTORY_HPP