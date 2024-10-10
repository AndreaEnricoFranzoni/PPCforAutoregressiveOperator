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

#include "KO_Traits.hpp"
#include "default_parameters.hpp"
#include "error_function.hpp"





//to do tag dispatching for correct way of creating train and vaidation set
template <DEF_PARAMS_PPC::cv_strat_type cv_strat>
using CV_STRAT_T = std::integral_constant<DEF_PARAMS_PPC::cv_strat_type, cv_strat>;


//to do tag dispatching for correct way of evaluating the error between prediction and validation
template <DEF_PARAMS_PPC::cv_err_eval_type err_eval>
using ERR_EVAL_T = std::integral_constant<DEF_PARAMS_PPC::cv_err_eval_type, err_eval>;




//strategy for cv type
//using cv_strategy_t     = std::map<std::vector<int>,std::vector<int>>;  
using cv_strategy_t     = std::vector<std::pair<std::vector<int>,std::vector<int>>>;
//train and valid set type
using train_valid_set_t = std::pair<KO_Traits::StoringMatrix,KO_Traits::StoringMatrix>;
//function to make prediction time
using pred_func_t       = std::function<KO_Traits::StoringVector(KO_Traits::StoringMatrix,double,double,int)>;


namespace CV_PPC
{

template<class D, DEF_PARAMS_PPC::cv_strat_type cv_strat, DEF_PARAMS_PPC::cv_err_eval_type err_eval>
class CV_KO_PPC_CRTP
{

private:
  //data in which you do CV
  KO_Traits::StoringMatrix m_Data;
  
  //strategy for doing CV
  //for each element of the map: first element is strategy (cols) for deciding the train set, second element is strategy for validations set
  cv_strategy_t train_validation_set_strategy(const KO_Traits::StoringMatrix &dataset, CV_STRAT_T<DEF_PARAMS_PPC::cv_strat_type::AUGMENTING_WINDOW>) const;
  //first element is the training set, the second one is the validation set
  train_valid_set_t train_validation_set(const std::pair<std::vector<int>,std::vector<int>> &strat, CV_STRAT_T<DEF_PARAMS_PPC::cv_strat_type::AUGMENTING_WINDOW>) const;
  
  //how to evaluate the error between prediction and validation
  double err_valid_set_eval(const KO_Traits::StoringVector &pred, const KO_Traits::StoringVector &valid, ERR_EVAL_T<DEF_PARAMS_PPC::cv_err_eval_type::MSE>) const;
  
  
public:
  
  //constructor
  CV_KO_PPC_CRTP(KO_Traits::StoringMatrix&& Data)
    : m_Data{std::forward<KO_Traits::StoringMatrix>(Data)}  {}
  
  //getter for data
  inline KO_Traits::StoringMatrix Data() const {return m_Data;}
  
  //how to define validation and training set
  cv_strategy_t train_validation_set_strategy(const KO_Traits::StoringMatrix &dataset) const { return train_validation_set_strategy(dataset, CV_STRAT_T<cv_strat>{});};
  train_valid_set_t train_validation_set(const std::pair<std::vector<int>,std::vector<int>> &strat) const { return train_validation_set(strat, CV_STRAT_T<cv_strat>{});};
  
  //error evaluation
  double err_valid_set_eval(const KO_Traits::StoringVector &pred, const KO_Traits::StoringVector &valid) const { return err_valid_set_eval(pred,valid,ERR_EVAL_T<err_eval>{});};
  
};
  
}   //end namespace CV_PPC

#include "CV_CRTP_train_valid_strat_imp.hpp"
#include "CV_CRTP_err_valid_eval.hpp"

#endif  //CV_CRTP_PPC_HPP