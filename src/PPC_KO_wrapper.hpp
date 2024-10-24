#ifndef PPC_KO_WRAPPER_HPP
#define PPC_KO_WRAPPER_HPP

#include <vector>
#include <tuple>
#include "utility"

#include "traits_ko.hpp"
#include "PPC_KO_include.hpp"
#include "mesh.hpp"



//base
template< DOM_DIM dom_dim, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval > 
class PPC_KO_wrapper
{
private:
  KO_Traits::StoringMatrix m_data;      //data
  Geometry::Mesh1D m_grid_func_data;    //grid where the there is the evaluation of the functional data
  results_t<valid_err_ret> m_results;                  //storing the results
  
public:
  template<typename STOR_OBJ,typename GRID>
  PPC_KO_wrapper(STOR_OBJ&& data, GRID&& grid_func_data)
    :  m_data{std::forward<STOR_OBJ>(data)}, m_grid_func_data{std::forward<GRID>(grid_func_data)}    {}
  
  /*!
   * Virtual destructor
   */
  virtual ~PPC_KO_wrapper() = default;
  
  /*!
   * Virtual method to call the correct version of KO at runtime
   */
  virtual void call_ko() = 0; 
  
  /*!
   * Getter for m_data
   */
  inline KO_Traits::StoringMatrix data() const {return m_data;};
  
  /*!
   * Getter for m_grid_func_data
   */
  inline Geometry::Mesh1D grid_func_data() const {return m_grid_func_data;};
  
  /*!
   * Getter for m_results
   */
  inline results_t<valid_err_ret> results() const {return m_results;};
  
  /*!
   * Setter for m_results
   */
  inline results_t<valid_err_ret> & results() {return m_results;};
};



//no cv
template< DOM_DIM dom_dim, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval > 
class PPC_KO_wrapper_no_cv  : public PPC_KO_wrapper<dom_dim,k_imp,valid_err_ret,cv_strat,cv_err_eval>
{
private:
  double m_alpha;                       //reg parameter
  int m_k;                              //number of PPCs
  double m_threshold_ppc;               //threshold for retaining the number of PPCs

public:
  //constructor if the k is already known (k_imp = YES)
  template<typename STOR_OBJ,typename GRID>
  PPC_KO_wrapper_no_cv(STOR_OBJ&& data, GRID&& grid_func_data, double alpha, int k)
    : PPC_KO_wrapper<dom_dim,k_imp,valid_err_ret,cv_strat,cv_err_eval>(std::move(data),std::move(grid_func_data)), m_alpha(alpha), m_k(k)    {}
  
  //constructor for finding the number of PPCs (k_imp = NO)
  template<typename STOR_OBJ,typename GRID>
  PPC_KO_wrapper_no_cv(STOR_OBJ&& data, GRID&& grid_func_data, double alpha, double threshold_ppc)
    : PPC_KO_wrapper<dom_dim,k_imp,valid_err_ret,cv_strat,cv_err_eval>(std::move(data),std::move(grid_func_data)), m_alpha(alpha), m_threshold_ppc(threshold_ppc)  {}
  
  void call_ko() override;
};



//cv alpha
template< DOM_DIM dom_dim, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval > 
class PPC_KO_wrapper_cv_alpha  : public PPC_KO_wrapper<dom_dim,k_imp,valid_err_ret,cv_strat,cv_err_eval>
{
private:
  std::vector<double> m_alphas;         //CV on alpha
  int m_k;
  double m_threshold_ppc;

public:
  //k given
  template<typename STOR_OBJ,typename GRID>
  PPC_KO_wrapper_cv_alpha(STOR_OBJ&& data, GRID&& grid_func_data, const std::vector<double> & alphas, int k)
    : PPC_KO_wrapper<dom_dim,k_imp,valid_err_ret,cv_strat,cv_err_eval>(std::move(data),std::move(grid_func_data)), m_alphas(alphas), m_k(k) {}
  
  //k to be found
  template<typename STOR_OBJ,typename GRID>
  PPC_KO_wrapper_cv_alpha(STOR_OBJ&& data, GRID&& grid_func_data, const std::vector<double> & alphas, double threshold_ppc)
    : PPC_KO_wrapper<dom_dim,k_imp,valid_err_ret,cv_strat,cv_err_eval>(std::move(data),std::move(grid_func_data)), m_alphas(alphas), m_threshold_ppc(threshold_ppc) {}
  
  void call_ko() override;
};



//cv k
template< DOM_DIM dom_dim, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval > 
class PPC_KO_wrapper_cv_k  : public PPC_KO_wrapper<dom_dim,k_imp,valid_err_ret,cv_strat,cv_err_eval>
{
private:
  double m_alpha;
  std::vector<int> m_k_s;
  double m_toll;
  
public:
  template<typename STOR_OBJ,typename GRID>
  PPC_KO_wrapper_cv_k(STOR_OBJ&& data, GRID&& grid_func_data, double alpha, const std::vector<int> & k_s, double toll)
    : PPC_KO_wrapper<dom_dim,k_imp,valid_err_ret,cv_strat,cv_err_eval>(std::move(data),std::move(grid_func_data)), m_alpha(alpha), m_k_s(k_s), m_toll(toll) {}
  
  void call_ko() override;
};



//cv alpha-k
template< DOM_DIM dom_dim, K_IMP k_imp, VALID_ERR_RET valid_err_ret, CV_STRAT cv_strat, CV_ERR_EVAL cv_err_eval > 
class PPC_KO_wrapper_cv_alpha_k  : public PPC_KO_wrapper<dom_dim,k_imp,valid_err_ret,cv_strat,cv_err_eval>
{
private:
  std::vector<double> m_alphas;
  std::vector<int> m_k_s;
  double m_toll;
  
public:
  template<typename STOR_OBJ,typename GRID>
  PPC_KO_wrapper_cv_alpha_k(STOR_OBJ&& data, GRID&& grid_func_data, const std::vector<double> &alphas, const std::vector<int> &k_s, double toll)
    : PPC_KO_wrapper<dom_dim,k_imp,valid_err_ret,cv_strat,cv_err_eval>(std::move(data),std::move(grid_func_data)), m_alphas(alphas), m_k_s(k_s), m_toll(toll) {}
  
  void call_ko() override;
};


#include "PPC_KO_wrapper_imp.hpp"

#endif  //PPC_KO_WRAPPER_HPP