#ifndef CV_CRTP_PPC_HPP
#define CV_CRTP_PPC_HPP

#include <iostream>
#include <functional>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <concepts>
#include <memory>
#include <utility>
#include <type_traits>

#include "traits_ko.hpp"
#include "strategy_cv.hpp"
#include "cv_eval_valid_err.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif



//to do tag dispatching for correct way of evaluating the error between prediction on training set and validation
template <CV_ERR_EVAL err_eval>
using ERR_EVAL_T = std::integral_constant<CV_ERR_EVAL, err_eval>;


//function to make prediction time
//using pred_func_t       = std::function<KO_Traits::StoringVector(KO_Traits::StoringMatrix,double,double,int)>;
//if k not imp: for using KO I need data, alpha and the threshold
//if k imp: for using KO I need data, alpha and k
using pred_func_k_yes_t = std::function<KO_Traits::StoringVector(KO_Traits::StoringMatrix,double,int,int)>;
using pred_func_k_no_t  = std::function<KO_Traits::StoringVector(KO_Traits::StoringMatrix,double,double,int)>;

template <K_IMP k_imp>
using pred_func_t = typename std::conditional<k_imp, pred_func_k_yes_t,pred_func_k_no_t>::type;


template< class D, CV_STRAT cv_strat, CV_ERR_EVAL err_eval, K_IMP k_imp, VALID_ERR_RET valid_err_ret >
class CV_base
{
  
private:
  //data in which you do CV
  KO_Traits::StoringMatrix m_Data;
  
  //strategy for CV
  cv_strategy<cv_strat> m_strategy;

  //how to evaluate the error between prediction on training set and validation set
  double err_valid_set_eval(const KO_Traits::StoringVector &pred, const KO_Traits::StoringVector &valid, int number_threads, ERR_EVAL_T<CV_ERR_EVAL::MSE>) const;
  
  //number of threads for OMP
  int m_number_threads;
  
public:
  
  //constructor
  template<typename STOR_OBJ,typename STRATEGY>
  CV_base(STOR_OBJ&& Data, STRATEGY && strategy, int number_threads)
    : m_Data{std::forward<STOR_OBJ>(Data)}, m_strategy{std::forward<STRATEGY>(strategy)}, m_number_threads(number_threads)  {}
  
  //getter for data
  inline KO_Traits::StoringMatrix Data() const {return m_Data;}
  
  //getter for the strategy
  inline cv_strategy<cv_strat> strategy() const {return m_strategy;}
  
  //getter for the number of threads
  inline int number_threads() const {return m_number_threads;}
  
  //error evaluation
  double err_valid_set_eval(const KO_Traits::StoringVector &pred, const KO_Traits::StoringVector &valid, int number_threads) const { return err_valid_set_eval(pred,valid,number_threads,ERR_EVAL_T<err_eval>{});};
  
};


#include "CV_err_valid_eval.hpp"

#endif  //CV_CRTP_PPC_HPP
