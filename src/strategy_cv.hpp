#ifndef STRAT_CVPPC_HPP
#define STRAT_CVPPC_HPP

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


//to do tag dispatching for correct way of creating train and validation set
template <CV_STRAT cv_strat>
using CV_STRAT_T = std::integral_constant<CV_STRAT, cv_strat>;



//strategy for cv type (a vector that contains pairs: first element col indices of train set, second of validation set)
using cv_strategy_t     = std::vector<std::pair<std::vector<int>,std::vector<int>>>;
//a single cv iter type (a pair that identifies a specific training and valdiation set)
using iter_cv_t         = std::pair<std::vector<int>,std::vector<int>>;
//train and valid sets
using train_valid_set_t = std::pair<KO_Traits::StoringMatrix,KO_Traits::StoringMatrix>;



template<CV_STRAT cv_strat>
class cv_strategy
{
private:
  //creating the strategy: method that is called in the constructor
  void train_validation_set_strategy(int min_dim_ts, int max_dim_ts, CV_STRAT_T<CV_STRAT::AUGMENTING_WINDOW>);
  //for each element: the elements of the first identify the training set, the second's ones the validation set
  cv_strategy_t m_strategy;
  //for a specific pair of train and validation sets, returns them
  train_valid_set_t train_validation_set(const KO_Traits::StoringMatrix &data, const iter_cv_t &strat, CV_STRAT_T<CV_STRAT::AUGMENTING_WINDOW>) const;
  
public:
  //constructor
  cv_strategy(int min_dim_ts, int max_dim_ts) { this->train_validation_set_strategy(min_dim_ts,max_dim_ts);}
  
  //getter for the strategy
  inline cv_strategy_t strategy() const {return m_strategy;}
  
  ////creating the strategy: method that is called in the constructor: tag dispatching to how to do it
  void train_validation_set_strategy(int n) { return train_validation_set_strategy(n, CV_STRAT_T<cv_strat>{});};
  //for a specific pair of train and validation sets, returns them
  train_valid_set_t train_validation_set(const KO_Traits::StoringMatrix &data, const iter_cv_t &strat) const { return train_validation_set(data, strat, CV_STRAT_T<cv_strat>{});};
  
};


#include "strategy_cv_imp.hpp"

#endif  //STRAT_CVPPC_HPP